//
// Created by David Heddle on 6/1/20.
//

#include "magfielddraw.h"
#include "svg.h"
#include "mapcolor.h"
#include "magfieldutil.h"

static int del = 2;

static int marginLeft =  50;
static int marginRight =  90;
static int marginTop =  50;
static int marginBottom =  50;
static int imageWidth;
static int imageHeight;
static int width;
static int height;
static int xmin;
static int xmax;
static int ymin;
static int ymax;
static svg* psvg;
static double coilThickness = 12; //cm

#define ROOT3OVER2 0.8660254037844386468

/**
 * Check whether the given point is withing the coils. This allows
 * us to exclude weird map values.
 * @param x the x coordinate.
 * @param y the y coordinate.
 * @return true if the given point is deemed to be (approximately) in the coils.
 */
static bool inCoils(double x, double y) {
    if (fabs(x) < coilThickness) {
        return true;
    }

    double dist = (x/2 - ROOT3OVER2*y);
    if (fabs(dist) < coilThickness) {
        return true;
    }

    dist = (x/2 + ROOT3OVER2*y);
    if (fabs(dist) < coilThickness) {
        return true;
    }

    return false;
}

/**
 * Common methof to draw a gradient on the right
 * @param psvg the svg object.
 * @param colorMap the color map.
 * @param width the width of the picture.
 * @param marginTop the top margin in pixels.
 */
static void gradient(svg* psvg, ColorMapPtr colorMap, int width) {

    char *label = (char *) malloc(255);
    int x = width - 75;
    int y = marginTop + 40;
    int nc = colorMap->numColors;
    int gw = 20;
    int gh = 4;
    for (int i = 0; i < colorMap->numColors; i++) {
        svgRectangle(psvg, gw, gh, x, y, colorMap->colors[i], "none", 0, 0, 0);
        bool lab = ((i == 0)  || (i == nc/4) || (i == nc/2) || (i == 3*nc/4));

        if (lab) {
            sprintf(label, " %-3.1f kG", colorMap->values[i]);
            svgText(psvg, x+gw+4, y+5, "times", 11, "black", "none", label);
        }
        y+=gh;
    }

    sprintf(label, " %-3.1f kG", colorMap->values[nc+1]);
    svgText(psvg, x+gw+4, y+5, "times", 11, "black", "none", label);

    svgRectangle(psvg, gw, nc*gh, width-75, marginTop+40, "none", "black", 1, 0, 0);
    free (label);
}

/**
 * Add the title, axis labels, and gradient.
 * @param psvg the svg object.
 * @param colorMap the colorMap being used.
 * @param title the plot title.
 * @param xlabel the x axis label.
 * @param ylabel  the y axis label.
 */
static void adornments(svg *psvg, ColorMapPtr colorMap, char *title, char *xlabel, char *ylabel) {
    //border
    svgRectangle(psvg, imageWidth+1, imageHeight+1, marginLeft, marginTop, "none", "black", 1, 0, 0);

    //title and axes labels
    svgText(psvg, marginLeft-10, 25, "times", 14, "black", "none", title);
    svgText(psvg, marginLeft + imageWidth/2, height-15, "times", 14, "black", "none", xlabel);
    svgRotatedText(psvg, 20, marginTop + imageHeight/2, "times", 14, "black", "none", -90, ylabel);

    //the gradient
    gradient(psvg, colorMap, width);
}

/**
 * Prepare to make a plot.
 * @param path the path to the svg file.
 * @param minX the minimum horizontal value of the plot.
 * @param maxX the maximum horizontal value of the plot.
 * @param minY the minimum vertical value of the plot.
 * @param maxY the maximum vertical value of the plot.
 */
static void init(char *path, int minX, int maxX, int minY, int maxY) {
    xmin = minX;
    xmax = maxX;
    ymin = minY;
    ymax = maxY;
    imageWidth = xmax - xmin;
    imageHeight = ymax - ymin;
    width = imageWidth + marginLeft + marginRight;
    height = imageHeight + marginTop + marginBottom;
    psvg = svgStart(path, width, height);
    svgFill(psvg, "#f0f0f0");
    fprintf(stdout, "\nStarting svg image creation for: [%s]", path);
}

/**
 * Create an SVG image of the fields at a fixed value of z.
 * @param path the path to the svg file.
 * @param title the title of the plot.
 * @param z the fixed value of z in cm.
 * @param fieldPtr torus the torus field (can be NULL).
 * @param fieldPtr solenoid the solenoid field (can be NULL).
 */
void createSVGImageFixedZ(char *path, char *title, double z, MagneticFieldPtr torus, MagneticFieldPtr solenoid) {

    ColorMapPtr colorMap;
    double rhoMax;

    colorMap = defaultColorMap();
     if (z > 99) {
         rhoMax = 360;
     }
     else {
         rhoMax = 240;
     }
    init(path, -rhoMax, rhoMax, -rhoMax, rhoMax);

    FieldValuePtr fieldValuePtr = (FieldValuePtr) malloc(sizeof (FieldValue));

    int y = ymin+del;
    while (y < ymax+del) {
        if ((y % 50) == 0) {
            fprintf(stdout, ".");
        }
        int yPic = y - ymin + marginTop; //x is vertical
        int x = xmin;
        while (x < xmax) {
            int xPic = x - xmin + marginLeft;
            bool inCoil = inCoils((double)x, (double) y);

            if (inCoil) {
                svgRectangle(psvg, del, del, xPic, yPic, "#555555", "none", 0, 0, 0);
            }
            else {
                getCompositeFieldValue(fieldValuePtr, x, y, z, torus, solenoid);
                double magnitude = fieldMagnitude(fieldValuePtr);
                char *color = getColor(colorMap, magnitude);
                svgRectangle(psvg, del, del, xPic, yPic, color, "none", 0, 0, 0);
            }

            x += del;
        }
        y += del;
    }

    free(fieldValuePtr);

    //labels
    char *label = (char *) malloc(255);

    int x = marginLeft - 6;
    y = ymax;
    while (y >= ymin) {
        int yPic = ymax - y + marginTop;

        sprintf(label, "%d", y);
        svgRotatedText(psvg, x, yPic+8, "times", 12, "black", "none", -90, label);
        svgLine(psvg, "#cccccc", 1, marginLeft, yPic, marginLeft+imageWidth, yPic);
        y -= 120;
    }


    x = xmin;
    y = marginTop + imageHeight + 20;
    while (x <= xmax) {

        int xPic = x - xmin + marginLeft;
        sprintf(label, "%d", x);
        svgText(psvg, xPic-12, y, "times", 12, "black", "none", label);
        svgLine(psvg, "#cccccc", 1, xPic, marginTop, xPic, marginTop+imageHeight);

        x += 120;
    }

    adornments(psvg, colorMap, title, "x (cm)", "y (cm)");

    free(label);
    svgEnd(psvg);

    fprintf(stdout, "done.\n");
}

/**
 * Create an SVG image of difference of two fields at a fixed value of phi. It plots
 * the magnitude of the vector difference between two fields. Presently the only
 * use is likely the difference between the symmetric and full torus maps.
 * @param path the path to the svg file.
 * @param title the title of the plot.
 * @param z the fixed value of z in cm.
 * @param fieldPtr field1 the first field.
 * @param fieldPtr fiel2 the second field.
 */
void createSVGImageFixedZDiff(char *path, char *title, double z, MagneticFieldPtr field1, MagneticFieldPtr field2) {
    double rhoMax;

    if (z > 99) {
        rhoMax = 300;
    }
    else {
        rhoMax = 240;
    }
    init(path, -rhoMax, rhoMax, -rhoMax, rhoMax);

    FieldValuePtr fieldValue1Ptr = (FieldValuePtr) malloc(sizeof (FieldValue));
    FieldValuePtr fieldValue2Ptr = (FieldValuePtr) malloc(sizeof (FieldValue));
    FieldValuePtr fieldValueDiffPtr = (FieldValuePtr) malloc(sizeof (FieldValue));

    //first get the max difference
    double maxDiff= 0;
    int y = ymin+del;
    while (y < ymax+del) {
        int x = xmin;
        while (x < xmax) {

            bool inCoil = inCoils((double)x, (double) y);

            if (!inCoil) {
                getFieldValue(fieldValue1Ptr, x, y, z, field1);
                getFieldValue(fieldValue2Ptr, x, y, z, field2);

                fieldValueDiffPtr->b1 = fieldValue2Ptr->b1 - fieldValue1Ptr->b1;
                fieldValueDiffPtr->b2 = fieldValue2Ptr->b2 - fieldValue1Ptr->b2;
                fieldValueDiffPtr->b3 = fieldValue2Ptr->b3 - fieldValue1Ptr->b3;

                double magnitude = fieldMagnitude(fieldValueDiffPtr);
                if (magnitude > maxDiff) {
                    maxDiff = magnitude;
                }
            }

            x += del;
        }
        y += del;
    }

    //now plot
    ColorMapPtr colorMap = getColorMap(1.01*maxDiff);
    y = ymin+del;
    while (y < ymax+del) {
        if ((y % 50) == 0) {
            fprintf(stdout, ".");
        }
        int yPic = y - ymin + marginTop; //x is vertical
        int x = xmin;
        while (x < xmax) {

            int xPic = x - xmin + marginLeft;

            bool inCoil = inCoils((double)x, (double) y);

            if (inCoil) {
                svgRectangle(psvg, del, del, xPic, yPic, "#555555", "none", 0, 0, 0);
            }
            else {
                getFieldValue(fieldValue1Ptr, x, y, z, field1);
                getFieldValue(fieldValue2Ptr, x, y, z, field2);

                fieldValueDiffPtr->b1 = fieldValue2Ptr->b1 - fieldValue1Ptr->b1;
                fieldValueDiffPtr->b2 = fieldValue2Ptr->b2 - fieldValue1Ptr->b2;
                fieldValueDiffPtr->b3 = fieldValue2Ptr->b3 - fieldValue1Ptr->b3;

                double magnitude = fieldMagnitude(fieldValueDiffPtr);

                char *color = getColor(colorMap, magnitude);
                svgRectangle(psvg, del, del, xPic, yPic, color, "none", 0, 0, 0);
            }
            x += del;
        }
        y += del;
    }

    free(fieldValue1Ptr);
    free(fieldValue2Ptr);
    free(fieldValueDiffPtr);


    //axes labels
    char *label = (char *) malloc(255);

    int x = marginLeft - 6;
    y = ymax;
    while (y >= ymin) {
        int yPic = ymax - y + marginTop;

        sprintf(label, "%d", y);
        svgRotatedText(psvg, x, yPic+8, "times", 12, "black", "none", -90, label);
        svgLine(psvg, "#cccccc", 1, marginLeft, yPic, marginLeft+imageWidth, yPic);
        y -= 50;
    }

    x = xmin;
    y = marginTop + imageHeight + 20;
    while (x <= xmax) {

        int xPic = x - xmin + marginLeft;
        sprintf(label, "%d", x);
        svgText(psvg, xPic-12, y, "times", 12, "black", "none", label);
        svgLine(psvg, "#cccccc", 1, xPic, marginTop, xPic, marginTop+imageHeight);

        x += 120;
    }

    adornments(psvg, colorMap, title, "x (cm)", "y (cm)");



    free(label);
    svgEnd(psvg);

    fprintf(stdout, "done.\n");
}

/**
 * Create an SVG image of the fields at a fixed value of phi.
 * @param path the path to the svg file.
 * @param title the title of the plot.
 * @param phi the fixed value of phi in degrees. For the canonical
 * sector 1 midplane, use phi = 0;
 * @param fieldPtr torus the torus field (can be NULL).
 * @param fieldPtr solenoid the solenoid field (can be NULL).
 */
void createSVGImageFixedPhi(char *path, char * title, double phi, MagneticFieldPtr torus, MagneticFieldPtr solenoid) {

    ColorMapPtr colorMap = defaultColorMap();

    double phiRad = toRadians(phi);
    double cosPhi = cos(phiRad);
    double sinPhi = sin(phiRad);

    int zmin = -100;
    int zmax = 500;
    int rmin = 0;
    int rmax = 360;

    init(path, zmin, zmax, rmin, rmax);

    FieldValuePtr fieldValuePtr = (FieldValuePtr) malloc(sizeof (FieldValue));

    int rho = rmin + del;
    while (rho < rmax + del) {
        if ((rho % 50) == 0) {
            fprintf(stdout, ".");
        }
        int rhoPic = marginTop + imageHeight - rho; //rho is vertical


        int z = zmin;
        while (z < zmax) {

            int zPic = z - zmin + marginLeft;

            getCompositeFieldValue(fieldValuePtr, rho*cosPhi, rho*sinPhi, z, torus, solenoid);
            double magnitude = fieldMagnitude(fieldValuePtr);

            char *color = getColor(colorMap, magnitude);
            svgRectangle(psvg, del, del, zPic, rhoPic, color, "none", 0, 0, 0);

            z += del;
        }
        rho += del;
    }

    free(fieldValuePtr);


    char *label = (char *) malloc(255);

    int z = marginLeft - 6;
    rho = rmin;
    while (rho <= rmax) {
        int xPic = rho - rmin + marginTop;

        sprintf(label, "%d", rmax - rho);
        svgRotatedText(psvg, z, xPic+8, "times", 12, "black", "none", -90, label);
        svgLine(psvg, "#cccccc", 1, marginLeft, xPic, marginLeft+imageWidth, xPic);
        rho += 60;
    }


    z = zmin;
    rho = marginTop + imageHeight + 20;
    while (z <= zmax) {

        int zPic = z - zmin + marginLeft;
        sprintf(label, "%d", z);
        svgText(psvg, zPic-12, rho, "times", 12, "black", "none", label);
        svgLine(psvg, "#cccccc", 1, zPic, marginTop, zPic, marginTop+imageHeight);

        z += 100;
    }

    adornments(psvg, colorMap, title, "z (cm)", "rho (cm)");

    free(label);
    svgEnd(psvg);

    fprintf(stdout, "done.\n");
}


/**
 * Create an SVG image of difference of two fields at a fixed value of phi. It plots
 * the magnitude of the vector difference between two fields. Presently the only
 * use is likely the difference between the symmetric and full torus maps.
 * @param path the path to the svg file.
 * @param title the title of the plot.
 * @param phi the fixed value of phi in degrees.
 * @param fieldPtr field1 the first field.
 * @param fieldPtr fiel2 the second field.
 */
void createSVGImageFixedPhiDiff(char *path, char *title, double phi, MagneticFieldPtr field1, MagneticFieldPtr field2) {

    double phiRad = toRadians(phi);
    double cosPhi = cos(phiRad);
    double sinPhi = sin(phiRad);

    int zmin = field1->zGridPtr->minVal;
    int zmax = field1->zGridPtr->maxVal;;
    int rmin = 0;
    int rmax = 360;

    init(path, zmin, zmax, rmin, rmax);

    FieldValuePtr fieldValue1Ptr = (FieldValuePtr) malloc(sizeof (FieldValue));
    FieldValuePtr fieldValue2Ptr = (FieldValuePtr) malloc(sizeof (FieldValue));
    FieldValuePtr fieldValueDiffPtr = (FieldValuePtr) malloc(sizeof (FieldValue));

    //first get the max difference
    double maxDiff= 0;

    int rho = rmin + del;
    while (rho < rmax + del) {
        int z = zmin;
        while (z < zmax) {

            getFieldValue(fieldValue1Ptr, rho*cosPhi, rho*sinPhi, z, field1);
            getFieldValue(fieldValue2Ptr, rho*cosPhi, rho*sinPhi, z, field2);

            fieldValueDiffPtr->b1 = fieldValue2Ptr->b1 - fieldValue1Ptr->b1;
            fieldValueDiffPtr->b2 = fieldValue2Ptr->b2 - fieldValue1Ptr->b2;
            fieldValueDiffPtr->b3 = fieldValue2Ptr->b3 - fieldValue1Ptr->b3;

            double magnitude = fieldMagnitude(fieldValueDiffPtr);
            if (magnitude > maxDiff) {
                maxDiff = magnitude;
            }

            z += del;
        }
        rho += del;
    }


    ColorMapPtr colorMap = getColorMap(1.01*maxDiff);
    //now plot
    rho = rmin + del;
    while (rho < rmax + del) {
        if ((rho % 50) == 0) {
            fprintf(stdout, ".");
        }
        int rhoPic = marginTop + imageHeight - rho; //rho is vertical


        int z = zmin;
        while (z < zmax) {

            int zPic = z - zmin + marginLeft;

            getFieldValue(fieldValue1Ptr, rho*cosPhi, rho*sinPhi, z, field1);
            getFieldValue(fieldValue2Ptr, rho*cosPhi, rho*sinPhi, z, field2);

            fieldValueDiffPtr->b1 = fieldValue2Ptr->b1 - fieldValue1Ptr->b1;
            fieldValueDiffPtr->b2 = fieldValue2Ptr->b2 - fieldValue1Ptr->b2;
            fieldValueDiffPtr->b3 = fieldValue2Ptr->b3 - fieldValue1Ptr->b3;

            double magnitude = fieldMagnitude(fieldValueDiffPtr);

            char *color = getColor(colorMap, magnitude);
            svgRectangle(psvg, del, del, zPic, rhoPic, color, "none", 0, 0, 0);

            z += del;
        }
        rho += del;
    }

    free(fieldValue1Ptr);
    free(fieldValue2Ptr);
    free(fieldValueDiffPtr);

    char *label = (char *) malloc(255);

    int z = marginLeft - 6;
    rho = rmin;
    while (rho <= rmax) {
        int xPic = rho - rmin + marginTop;

        sprintf(label, "%d", rmax - rho);
        svgRotatedText(psvg, z, xPic+8, "times", 12, "black", "none", -90, label);
        svgLine(psvg, "#cccccc", 1, marginLeft, xPic, marginLeft+imageWidth, xPic);
        rho += 60;
    }


    z = zmin;
    rho = marginTop + imageHeight + 20;
    while (z <= zmax) {

        int zPic = z - zmin + marginLeft;
        sprintf(label, "%d", z);
        svgText(psvg, zPic-12, rho, "times", 12, "black", "none", label);
        svgLine(psvg, "#cccccc", 1, zPic, marginTop, zPic, marginTop+imageHeight);

        z += 100;
    }

    adornments(psvg, colorMap, title, "z (cm)", "rho (cm)");

    free(label);
    svgEnd(psvg);

    fprintf(stdout, "done.\n");
}



