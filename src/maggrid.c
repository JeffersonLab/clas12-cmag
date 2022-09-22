//
//  maggrid.c
//  @file
//
//  Created by David Heddle on 5/22/20.
//

#include "maggrid.h"
#include "magfieldutil.h"
#include "munittest.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>


/**
 * Create a uniform (equally spaced) coordinate coordinate grid.
 * @param the name of the coordinate, e.g. "phi".
 * @param minVal the minimum value of the grid.
 * @param maxVal the maximum value of the grid.
 * @param num the number of points on the grid, including the ends.
 * @return a pointer to the coordinate grid.
 */
GridPtr createGrid(const char *name, const double minVal, const double maxVal,
        const unsigned int num) {
    GridPtr gridPtr;

    gridPtr = (GridPtr) malloc(sizeof(Grid));

    stringCopy(&(gridPtr->name), name);
    gridPtr->minVal = minVal;
    gridPtr->maxVal = maxVal;
    gridPtr->numPoints = num;

    if (num < 2) { //happens for solenoid q1 only
        gridPtr->delta = INFINITY;
        gridPtr->values = (double *) malloc(sizeof(double));
        gridPtr->values[0] = 0;
    } else {
        gridPtr->delta = (maxVal - minVal) / (num - 1);
        gridPtr->values = (double *) malloc(num * sizeof(double));

        gridPtr->values[0] = minVal;
        gridPtr->values[num - 1] = maxVal;

        for (int i = 1; i < (num - 1); i++) {
            gridPtr->values[i] = minVal + i * gridPtr->delta;
        }
    }

    return gridPtr;
}

/**
 * Get a string representation of the grid.
 * @param gridPtr the pointer to the coordinate grid.
 * @return  a string representation of the grid.
 */
char *gridStr(GridPtr gridPtr) {
    char *str = (char*) malloc(128);

    sprintf(str, "%3s min: %6.1f  max: %6.1f  Np: %4d  delta: %6.1f",
            gridPtr->name, gridPtr->minVal, gridPtr->maxVal, gridPtr->numPoints, gridPtr->delta);
    return str;
}

/**
 * Print the grid values
 * @param gridPtr the pointer to the coordinate grid.
 * @param stream the stream to print to, e.g., stdout
 */
void printGrid(GridPtr gridPtr, FILE * stream) {
    fprintf(stream, "\nValues for grid [%s]\n", gridPtr->name);
}

/**
 * Get the index of a value.
 * @param gridPtr the pointer to the coordinate grid.
 * @param val the value to index.
 * @return the index, [0..N-2] where, or -1 if out of bounds. The value
 * should be bounded by values[index] and values[index+1].
 */
int getIndex(const GridPtr gridPtr, const double val) {

    int index;

    if (gridPtr->numPoints < 2) { //solenoid phi (q1) grid
        return 0;
    }

    if ((val < gridPtr->minVal) || (val > gridPtr->maxVal)) {
        index = -1;
    }

    double fract = (val - gridPtr->minVal) / gridPtr->delta;
    index = (int) (fract);

    //guard against weird possibility of landing right on the max grid value
    if (index > (gridPtr->numPoints - 2)) {
        index = gridPtr->numPoints - 2;
    }

    return index;
}

/**
 * Get the value of the grid at a given index
 * @param gridPtr the pointer to the grid
 * @param index the index
 * @return the value of the grid at the given index, or NAN if
 * the index is out of range
 */
double valueAtIndex(GridPtr gridPtr, int index) {
    if ((index < 0) || (index >= gridPtr->numPoints)) {
        return NAN;
    }

    return gridPtr->values[index];
}

/**
 * A unit test for the coordinate grid code.
 * @return an error message if the test fails, or NULL if it passes.
 */
char *gridUnitTest() {
    //for the test, we will create a CLAS like grid, and
    //generate a bunch of random points, and make sure
    //they always give the correct index.

    double minVal = -300;
    double maxVal = 300;
    unsigned int numPoints = 1201;

    GridPtr gridPtr = createGrid("TestGrid", minVal, maxVal, numPoints);

    int numTestPoints = 100000;
    srand48(time(0)); //seed random
    double range = maxVal - minVal;

    for (int i = 0; i < numTestPoints; i++) {


        double val = randomDouble(minVal, maxVal);
        int index = getIndex(gridPtr, val);

        //result should be true if we pass
        bool result = (index >= 0) && (index <= (gridPtr->numPoints - 2));

        if (!result) {
            fprintf(stdout, "BAD INDEX: Val = %f Index = %d \n", val, index);
        }

        mu_assert("Bad index", result);
    }
    fprintf(stdout, "\nPASSED gridUnitTest\n");

    return NULL;
}
