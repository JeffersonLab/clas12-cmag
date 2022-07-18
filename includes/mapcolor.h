//
// Created by David Heddle on 5/30/20.
//

#ifndef CMAG_MAPCOLOR_H
#define CMAG_MAPCOLOR_H

typedef struct colormap *ColorMapPtr;

typedef struct colormap {
    int numColors; //number of colors in map
    char *tooSmallColor; //for values < min value
    char *tooBigColor; //for values > max value

    double values[74]; //n+1 descending ordered values
    char *colors[73]; //n colors of map
} ColorMap;

//external prototypes
extern char *getColor(ColorMapPtr, double);

extern ColorMapPtr getColorMap(double);
extern ColorMapPtr defaultColorMap();
#endif //CMAG_MAPCOLOR_H
