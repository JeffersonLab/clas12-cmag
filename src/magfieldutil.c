//
//  magfieldutil.c
//  cMag
//
//  Created by David Heddle on 5/22/20.
//

#include "magfield.h"
#include "magfieldutil.h"
#include "munittest.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//used for comparing real numbers
double const TINY = 1.0e-10;

//used for degrees <--> radians;
const double PIOVER180 = M_PI/180.;

//some strings for prints
const char *csLabels[] = { "cylindrical", "Cartesian" };
const char *lengthUnitLabels[] = { "cm", "m" };
const char *angleUnitLabels[] = { "degrees", "radians" };
const char *fieldUnitLabels[] = { "kG", "G", "T" };

//local prototypes
static void freeGrid(GridPtr gridPtr);
static void resetCell(Cell3D *);

/**
 * Convert an angle from radians to degrees.
 * @param angRad  the angle in radians.
 * @return the angle in degrees.
 */
double toDegrees(double angRad) {
    return  angRad/PIOVER180;
}

/**
 * Convert an angle from degrees to radians.
 * @param angDeg the angle in degrees.
 * @return the angle in radians.
 */
double toRadians(double angDeg) {
    return angDeg * PIOVER180;
}

/**
 * The usual test to see if two floating point numbers are close
 * enough to be considered equal. Test accuracy depends on the
 * global const TINY, set to 1.0e-10.
 * @param v1 one value.
 * @param v2 another value.
 * @return true if the values are close enough to be considered equal.
 */
bool sameNumber(double v1, double v2) {
    if (v1 == v2) {
        return true;
    }

    double del = fabs(v2 - v1);
    v1 = fabs(v1);
    v2 = fabs(v2);
    double vmax = (v2 > v1) ? v2 : v1;
    return del/vmax < TINY;
}

/**
 * Converts 2D Cartesian coordinates to polar. This is used because the two coordinate systems we use
 * are Cartesian and cylindrical, whose 3D transformations are equivalent to 2D Cartesian to polar. Note
 * the azimuthal angle output is in degrees, not radians.
 * @param x the x component.
 * @param y the y component.
 * @param phi will hold the angle, in degrees, in the range [0, 360).
 * @param rho the longitudinal component.
 */
void cartesianToCylindrical(const double x, const double y, double *phi, double *rho) {
    *phi = atan2(y, x);
    *phi = toDegrees(*phi);
    normalizeAngle(phi);
    *rho = sqrt(x*x + y*y);
}

/**
 * Converts polar coordinates to 2D Cartesian. This is used because the two coordinate systems we use
 * are Cartesian and cylindrical, whose 3D transformations are equivalent to 2D polar to Cartesian. Note
 * the azimuthal angle input is in degrees, not radians.
 * @param x will hold the x component.
 * @param y will hold the y component.
 * @param phi the azimuthal angle, in degrees.
 * @param rho the longitudinal component.
 */
void cylindricalToCartesian(double *x, double *y, const double phi, const double rho) {

    double dphi = toRadians(phi);
    *x = rho*cos(dphi);
    *y = rho*sin(dphi);
}

/**
 * This will normalize an angle in degrees. We use for normaliztion that
 * the angle should be in the range [0, 360).
 * @param angDeg the angle in degrees. It will be normalized.
 */
void normalizeAngle(double *angDeg) {
    while(*angDeg < 0) {
        *angDeg += 360;
    }
    while(*angDeg >= 360) {
        *angDeg -= 360;
    }
}

/**
 * Get the magnitude of a field value.
 * @param fvPtr a pointer to a field value
 * @return return: the magnitude of a field value.
 */
double fieldMagnitude(FieldValue *fvPtr) {
    double b1 = fvPtr->b1;
    double b2 = fvPtr->b2;
    double b3 = fvPtr->b3;
    return sqrt(b1 * b1 + b2 * b2 + b3 * b3);
}
/**
 * Print a summary of the map for diagnostics and debugging.
 * @param fieldPtr the pointer to the map.
 * @param stream a file stream, e.g. stdout.
 */
void printFieldSummary(MagneticFieldPtr fieldPtr, FILE *stream) {
    FieldMapHeaderPtr headerPtr = fieldPtr->headerPtr;

    fprintf(stream, "\n========================================\n");
    fprintf(stream, "%s: [%s]\n", (fieldPtr->type == TORUS) ? "TORUS" : "SOLENOID", fieldPtr->path);
    fprintf(stream, "Created: %s\n", fieldPtr->creationDate);
    fprintf(stream, "Symmetric: %s\n", fieldPtr->symmetric ? "true" : "false");
    fprintf(stream, "scale factor: %-6.2f\n", fieldPtr->scale);

    //print the grid info for the three coordinate grids
    fprintf(stream, "%s\n", gridStr(fieldPtr->phiGridPtr));
    fprintf(stream, "%s\n", gridStr(fieldPtr->rhoGridPtr));
    fprintf(stream, "%s\n", gridStr(fieldPtr->zGridPtr));

    fprintf(stream, "num field values: %d\n", fieldPtr->numValues);
    fprintf(stream, "grid cs: %s\n", csLabels[headerPtr->gridCS]);
    fprintf(stream, "field cs: %s\n", csLabels[headerPtr->fieldCS]);
    fprintf(stream, "length unit: %s\n",
            lengthUnitLabels[headerPtr->lengthUnits]);
    fprintf(stream, "angular unit: %s\n",
            angleUnitLabels[headerPtr->angleUnits]);
    fprintf(stream, "field unit: %s\n", fieldUnitLabels[headerPtr->fieldUnits]);

    //now the metrics
    fprintf(stream, "max field at index: %d\n",
            fieldPtr->metricsPtr->maxFieldIndex);


    fprintf(stream, "max field magnitude: %-10.6f %s\n",
            fieldPtr->metricsPtr->maxFieldMagnitude, fieldUnits(fieldPtr));

    FieldValuePtr fieldValPtr= getFieldAtIndex(fieldPtr, fieldPtr->metricsPtr->maxFieldIndex);
    fprintf(stdout, "max field vector");
    printFieldValue(fieldValPtr, stdout);

    //get the location of the max field
    int phiIndex, rhoIndex, zIndex;

    invertCompositeIndex(fieldPtr, fieldPtr->metricsPtr->maxFieldIndex, &phiIndex, &rhoIndex, &zIndex);
    double phi = fieldPtr->phiGridPtr->values[phiIndex];
    double rho = fieldPtr->rhoGridPtr->values[rhoIndex];
    double z = fieldPtr->zGridPtr->values[zIndex];
    fprintf(stdout, "max field location (phi, rho, z) = (%-6.2f, %-6.2f, %-6.2f)\n", phi, rho, z);
    fprintf(stdout, "avg field magnitude: %-10.6f %s", fieldPtr->metricsPtr->avgFieldMagnitude, fieldUnits(fieldPtr));
}

 /**
  * Get the field units of the magnetic field
  * @param fieldPtr a pointer to the field.
  * @return a string representing the field units, e.g. "kG".
  */
const char *fieldUnits(MagneticFieldPtr fieldPtr) {
    unsigned int funits = fieldPtr->headerPtr->fieldUnits;
    return fieldUnitLabels[funits];
}

 /**
  * Get the length units of the magnetic field
  * @param fieldPtr a pointer to the field.
  * @return a string representing the length units, e.g. "cm".
  */
const char *lengthUnits(MagneticFieldPtr fieldPtr) {
    int lunits = fieldPtr->headerPtr->lengthUnits;
    return lengthUnitLabels[lunits];
}

 /**
  * Print the components and magnitude of field value.
  * @param fvPtr a pointer to the field value
  * @param stream a file stream, e.g. stdout.
  */
void printFieldValue(FieldValuePtr fvPtr, FILE *stream) {
    fprintf(stream, "(%-9.5f, %-9.5f, %-9.5f), magnitude: %12.5f\n", fvPtr->b1,
            fvPtr->b2, fvPtr->b3, fieldMagnitude(fvPtr));
}

 /**
  * Allocate a field map with no content.
  * @return a pointer to an empty field map structure.
  */
MagneticFieldPtr createFieldMap() {
     MagneticFieldPtr fieldPtr = (MagneticFieldPtr) malloc(sizeof(MagneticField));
     fieldPtr->metricsPtr = (FieldMetricsPtr) malloc(sizeof(FieldMetrics));
     fieldPtr->scale = 1;
     fieldPtr->shiftX = 0;
     fieldPtr->shiftY = 0;
     fieldPtr->shiftZ = 0;

     return fieldPtr;
}


/**
 * Free the memory associated with a field map.
 * @param fieldPtr a pointer to the field
 */
void freeFieldMap(MagneticFieldPtr fieldPtr) {
    free(fieldPtr->metricsPtr);
    freeGrid(fieldPtr->phiGridPtr);
    freeGrid(fieldPtr->rhoGridPtr);
    freeGrid(fieldPtr->zGridPtr);

    if (fieldPtr->type == TORUS) {
        freeCell3D(fieldPtr->cell3DPtr);
    }
    else {
        freeCell2D(fieldPtr->cell2DPtr);
    }
    free(fieldPtr);
}

/**
 * Free the memory associated with a coordinate grid.
 * @param gridPtr the poiner to deallocate.
 */
static void freeGrid(GridPtr gridPtr) {
    free(gridPtr->name);
    free(gridPtr->values);
    free(gridPtr);
}

/**
 * Copy a string and create the pointer
 * @param dest on input a pointer to an unallocated string.
 * On output the string will be allocated and contain a copy of src.
 * @param src the string to be copied.
 */
void stringCopy(char **dest, const char *src) {
    unsigned long len = strlen(src);
    *dest = (char*) malloc(len + 1);
    strncpy(*dest, src, len);
}

/**
 * Obtain a random int in an inclusive range[minVal, maxVal]. Used for testing.
 * @param minVal the minimum value
 * @param maxVal the maximum Value;
 * @return
 */
int randomInt(int minVal, int maxVal) {
    int del = maxVal - minVal;
    int randVal = rand();

    return minVal + (randVal % (del + 1) );
}

/**
 * Obtain a random double in the range[minVal, maxVal]. Used for testing.
 * @param minVal the minimum value
 * @param maxVal the maximum Value;
 * @return
 */
double randomDouble(double minVal, double maxVal) {
    double range = maxVal - minVal;
    return minVal + range*drand48();
}

/**
 * A unit test for the conversions
 * @return an error message if the test fails, or NULL if it passes.
 */
char *conversionUnitTest() {
    srand48(time(0)); //seed random

    int num = 1000000;

    double x, y, phi, rho, tx, ty;

    for (int i = 0; i < 10000; i++) {
        x = randomDouble(-100, 600);
        y = randomDouble(-100, 600);

        cartesianToCylindrical(x, y, &phi, &rho);
        cylindricalToCartesian(&tx, &ty, phi, rho);

        bool result = sameNumber(x, tx) && sameNumber(y, ty);

        if (!result) {
            fprintf(stdout, "Conversions did not invert x: [%-6.3f to %-6.3f] y: [%-6.3f to %-6.3f] \n", x, tx, y, ty);
        }

        mu_assert("Conversions did not invert", result);
    }

    fprintf(stdout, "\nPASSED conversionUnitTest\n");
    return NULL;
}

/**
 * A unit test for the random number generator
 * @return an error message if the test fails, or NULL if it passes.
 */
char *randomUnitTest() {
    srand48(time(0)); //seed random

    int minVal = 0;
    int maxVal = 301;

    int count = 100000;
    bool result;

    for (int i = 0; i < count; i++) {
        int val = randomInt(minVal, maxVal);

        //result should be true if we pass
        bool result = (val >= minVal) && (val <= maxVal);

        if (!result) {
            fprintf(stdout, "OUT OF RANGE: val = %d range: [%d, %d] \n", val, minVal, maxVal);
        }

        mu_assert("Random number generated out of range.", result);
    }

    fprintf(stdout, "\nPASSED randomUnitTest\n");
    return NULL;
}

/**
 * Reset the cell to "contains nothing" values.
 * @param cell the cell to reset.
 */
static void resetCell(Cell3D *cell) {
    cell->phiMax = -INFINITY;
    cell->phiMin = INFINITY;
    cell->rhoMax = -INFINITY;
    cell->rhoMin = INFINITY;
    cell->zMax = -INFINITY;
    cell->zMin = INFINITY;
}

/**
 * Get the hex color string from color components
 * @param colorStr must be at last 8 characters
 * @param r the red component [0..255]
 * @param g the green component [0..255]
 * @param b the blue component [0..255]
 */
void colorToHex(char * colorStr, int r, int g, int b) {
    sprintf(colorStr, "#%02x%02x%02x", r, g, b);
}
