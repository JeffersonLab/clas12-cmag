//
// Created by David Heddle on 5/30/20.
// Used only to make the svg picture of the map
//

#include "mapcolor.h"
#include "magfieldutil.h"
#include <stdlib.h>

/**
 * Get a color from a color map.
 * @param cmapPtr a valid pointer to a color map.
 * @param value the value to convert into a color.
 * @return the color in "#rrggbb" format.
 */
char *getColor(ColorMapPtr cmapPtr, double value) {

    //the values are in descending order

    if (value > cmapPtr->values[0]) {
        return "black";
    }
    if (value <= cmapPtr->values[cmapPtr->numColors]) {
        return "white";
    }

    int index;
    for (index = 0; index < cmapPtr->numColors; index++) {
        if (cmapPtr->values[index+1] < value) {
            break;
        }
    }
    return cmapPtr->colors[index];
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

/**
 * Get a color map for a given maximum value. Unlike the default
 * map, here the values are generated linearly.
 * @return a color map corresponding to a range [0..maxVal].
 */
ColorMapPtr getColorMap(double maxVal) {

    ColorMapPtr cmapPtr = (ColorMapPtr) malloc(sizeof(ColorMap));

    int colorlen = 73;
    cmapPtr->numColors = colorlen;
    stringCopy(&(cmapPtr->tooSmallColor), "#ffffff");
    stringCopy(&(cmapPtr->tooBigColor), "#000000");

    int r, g, b;

    b = 0;
    r = 255;
    int i = 0;
    for (int j = 0; j <19; j++) {
        int del = (int)((j*255.)/19);
        g = del;
        cmapPtr->colors[i] = (char *) malloc(8);
        colorToHex(cmapPtr->colors[i], r, g, b);
        i++;
    }

    for (int j = 0; j <36; j++) {
        int del = (int)((j*255.)/36);
        r = 255 - del;
        g = 255 - del/2;
        b = del;
        cmapPtr->colors[i] = (char *) malloc(8);
        colorToHex(cmapPtr->colors[i], r, g, b);
        i++;
    }


    g = 255;
    b = 255;
    for (int j = 0; j <18; j++) {
        int del = (int)((j*255.)/18);
        r = del;
        g = 128 + del/2;
        cmapPtr->colors[i] = (char *) malloc(8);
        colorToHex(cmapPtr->colors[i], r, g, b);
        i++;
    }

    double del = maxVal/colorlen;
    //linear
    for (int i = 0; i <= colorlen; i++) {
        cmapPtr->values[i] = maxVal*(1. - ((double)i)/colorlen);
    }

    return cmapPtr;
}

/**
 * Get the default color map optimized for displaying torus and solenoid
 * @return the default color map.
 */
ColorMapPtr defaultColorMap() {

    ColorMapPtr cmapPtr = (ColorMapPtr) malloc(sizeof(ColorMap));

    int colorlen = 73;
    cmapPtr->numColors = colorlen;
    stringCopy(&(cmapPtr->tooSmallColor), "#ffffff");
    stringCopy(&(cmapPtr->tooBigColor), "#000000");

    char* colors[] = {
            "#7f007f","#8e0070","#9c0062","#aa0054","#b80046","#c70038",
            "#d5002a","#e3001c","#f1000e","#ff0000","#ff0f00","#ff1d00",
            "#ff2b00","#ff3900","#ff4800","#ff5600","#ff6400","#ff7200",
            "#ff8000","#ff8f00","#ff9d00","#ffab00","#ffb900","#ffc700",
            "#ffd500","#ffe300","#fff100","#ffff00","#f1f306","#e3e80c",
            "#d6dd11","#c8d217","#bac71d","#adbc22","#9fb128","#91a62e",
            "#849b33","#799c32","#6e9d30","#649e30","#599f2f","#4f9f2e",
            "#44a02d","#3aa12c","#2fa22b","#25a22a","#35a93f","#44b054",
            "#54b769","#63be7e","#73c593","#82cca8","#92d3bd","#a1dad2",
            "#b0e0e6","#b9e4e5","#c2e7e5","#cbeae5","#d4eee4","#d7f1e4",
            "#dcf1e4","#d9f1e4","#dbf1e4","#ddf1e4","#dff1e4","#f1f1e4",
            "#f3f1e4","#f5f1e4","#f7f1e4","#f9f1e4","#fbf1e4","#fdf1e4",
            "#fef4f2"};



    double values[] = {
            66.00000, 61.63079, 57.55083, 53.74095, 50.18330, 46.86116,
            43.75894, 40.86209, 38.15702, 35.63102, 33.27224, 31.06962,
            29.01280, 27.09215, 25.29865, 23.62388, 22.05997, 20.59960,
            19.23591, 17.96249, 16.77337, 15.66297, 14.62608, 13.65783,
            12.75368, 11.90939, 11.12098, 10.38477, 9.69730, 9.05534,
            8.45588, 7.89610, 7.37337, 6.88526, 6.42945, 6.00382,
            5.60637, 5.23523, 4.88866, 4.56503, 4.26282, 3.98062,
            3.71711, 3.47103, 3.24125, 3.02668, 2.82631, 2.63921,
            2.46450, 2.30135, 2.14900, 2.00674, 1.87389, 1.74984,
            1.63400, 1.52583, 1.42482, 1.33050, 1.24242, 1.16017,
            1.08337, 1.01165, 0.94468, 0.88214, 0.82375, 0.76921,
            0.71829, 0.67074, 0.62634, 0.58488, 0.54616, 0.51000,
            0.47624, 0.1};

    for (int i = 0; i < colorlen; i++) {
        cmapPtr->colors[i] = colors[i];
    }

    for (int i = 0; i <= colorlen; i++) {
        cmapPtr->values[i] = values[i];
    }

    return cmapPtr;
}



