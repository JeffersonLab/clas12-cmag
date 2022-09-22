//
//  maggrid.h
//  cMag
//
//  Created by David Heddle on 5/22/20.
//  Copyright © 2020 David Heddle. All rights reserved.
//

#ifndef maggrid_h
#define maggrid_h

#include <stdio.h>

typedef struct grid *GridPtr;

/**
 * Holds the uniformly spaced grid values for a coordinate.
 */
typedef struct grid {
        char *name;    //the name of the coordinate, e.g., "phi".
        double minVal;  //the min value of the coordinate
        double maxVal;  //the max value of the coordinate
        unsigned int numPoints; //the number of vals including the ends
        double delta; // (max-min)/(n-1)
        double *values; //values of the coordinates
} Grid;

//external prototypes
extern GridPtr createGrid(const char*, double, double, unsigned int);
extern char *gridStr(GridPtr);
extern void printGrid(GridPtr, FILE *);
extern double valueAtIndex(GridPtr, int);
extern char *gridUnitTest(void);
extern int getIndex(const GridPtr, const double);

#endif /* maggrid_h */
