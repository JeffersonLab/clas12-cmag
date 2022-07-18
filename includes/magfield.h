//
//  magfield.h
//  cMag
//
//  Created by David Heddle on 5/22/20.
//

#ifndef magfield_h
#define magfield_h


#include "maggrid.h"
#include <stdio.h>
#include <stdbool.h>

#define IS_INDEXABLE(arg) (sizeof(arg[0]))
#define IS_ARRAY(arg) (IS_INDEXABLE(arg) && (((void *) &arg) == ((void *) arg)))
#define ARRAYSIZE(arr) (IS_ARRAY(arr) ? (sizeof(arr) / sizeof(arr[0])) : 0)


#define FMDEBUG 1

#define ROOT3OVER2 0.8660254037844386468

// debug print
#define debugPrint(fmt, ...) \
  do { if (FMDEBUG) fprintf(stdout, fmt, __VA_ARGS__); } while (0)

//magic word used to test if byte swapping is required
#define MAGICWORD 0xced

//pointers to structures defined below
typedef struct fieldmapheader *FieldMapHeaderPtr;
typedef struct magneticfield *MagneticFieldPtr;
typedef struct fieldmetrics *FieldMetricsPtr;
typedef struct fieldvalue *FieldValuePtr;
typedef struct cell3d *Cell3DPtr;
typedef struct cell2d *Cell2DPtr;

//some strings for prints
extern const char *csLabels[];
extern const char *lengthUnitLabels[];
extern const char *angleUnitLabels[];
extern const char *fieldUnitLabels[];

//used for unit testing
extern MagneticFieldPtr testFieldPtr;

//the header of the binary files
typedef struct fieldmapheader {
    unsigned int magicWord;   // if not 0xced = 3309 then byteswap
    unsigned int gridCS;      // 0 for cylindrical, 1 for Cartesian
    unsigned int fieldCS;     // 0 for cylindrical, 1 for Cartesian
    unsigned int lengthUnits; // 0 for cm, 1 for meters
    unsigned int angleUnits;  // 0 for degrees, 1 for radians
    unsigned int fieldUnits;  // 0 for kG, 1 for G, 2 for T
    float q1min;  // min value of slowest varying coordinate
    float q1max;  // max value of slowest varying coordinate
    unsigned int nq1; // numColors equally spaced in q1 direction including ends
    float q2min;  // min value of medium varying coordinate
    float q2max;  // max value of medium varying coordinate
    unsigned int nq2; // numColors equally spaced in q2 direction including ends
    float q3min;  // min value of fastest varying coordinate
    float q3max;  // max value of fastest varying coordinate
    unsigned int nq3; // numColors equally spaced in q3 direction including ends
    int cdHigh; //high word of unix creation date of map
    int cdLow; //low word of unix creation date of map
    unsigned int reserved3; //reserved
    unsigned int reserved4; //reserved
    unsigned int reserved5; //reserved
} FieldMapHeader;

//holds a single field value
typedef struct fieldvalue {
    float b1; //1st component of field
    float b2; //2nd component of field
    float b3; //3rd component of field
} FieldValue;

//holds some metrics of the field
typedef struct fieldmetrics {
    unsigned int maxFieldIndex; //the index where we find the max field value
    double maxFieldMagnitude; //the value of the max field index
    double avgFieldMagnitude; //the average field magnitude

} FieldMetrics;

//3d cell is used by torus
typedef struct cell3d {

    MagneticFieldPtr fieldPtr; //the owner of the cell

    //cell boundaries
    double phiMin;  //minimum value of phi for the cell
    double phiMax;  //maximum value of phi for the cell

    double rhoMin;  //minimum value of rho for the cell
    double rhoMax;  //maximum value of rho for the cell

    double zMin;  //minimum value of z for the cell
    double zMax;  //maximum value of z for the cell

    int phiIndex; //phi index of cell
    int rhoIndex; //rho index of cell
    int zIndex;   //z index of cell

    double phiNorm; //cached value to speed up evaluation
    double rhoNorm; //cached value to speed up evaluation
    double zNorm; //cached value to speed up evaluation

    // space used by the cell.
    double f[3];
    double g[3];
    double a[8];

    FieldValuePtr b[2][2][2]; //field at 4 corners of cell
} Cell3D;

//2d cell is used by solenoid
typedef struct cell2d {

    MagneticFieldPtr fieldPtr; //the owner of the cell

    //cell boundaries
    double rhoMin;  //minimum value of rho for the cell
    double rhoMax;  //maximum value of rho for the cell

    double zMin;  //minimum value of z for the cell
    double zMax;  //maximum value of z for the cell

    int rhoIndex; //rho index of cell
    int zIndex;   //z index of cell

    double rhoNorm; //cached value to speed up evaluation
    double zNorm; //cached value to speed up evaluation

    FieldValuePtr b[2][2]; //field at 4 corners of cell

} Cell2D;

typedef enum {TORUS, SOLENOID} FieldType;
typedef enum {INTERPOLATION, NEAREST_NEIGHBOR} Algorithm;

//holds the entire field map
typedef struct magneticfield {
    FieldMapHeaderPtr headerPtr; //pointer to the header data

    char *path; //the path to the file
    char *name; //a descriptive name
    bool symmetric; //is this a symmetric grid (solenoid always is)
    FieldType  type;
    char *creationDate;  //date the map was created
    unsigned int numValues;  //total number of field values

    //a grid pointer for each coordinate
    GridPtr phiGridPtr;
    GridPtr rhoGridPtr;
    GridPtr zGridPtr;

    FieldMetricsPtr metricsPtr; //some field metrics

    Cell3DPtr cell3DPtr;  //the cell, or probe for torus
    Cell2DPtr cell2DPtr;  //the cell, or probe for solenoid

    double scale; //scale factor of the field

    double shiftX; //misplacement shift in the x direction (cm)
    double shiftY; //misplacement shift in the y direction (cm)
    double shiftZ; //misplacement shift in the z direction (cm)

    //some auxiliary data to cache
    unsigned int N23; // for faster indexing

    //use 1D array which will require manual indexing
    FieldValue *fieldValues;
} MagneticField;

// external function prototypes
extern int getCompositeIndex(MagneticFieldPtr, int, int, int);
extern void invertCompositeIndex(MagneticFieldPtr fieldPtr, int index, int *phiIndex, int *rhoIndex, int *zIndex);
extern char *compositeIndexUnitTest();
extern char *containsUnitTest();
extern char *nearestNeighborUnitTest();
extern FieldValuePtr getFieldAtIndex(MagneticFieldPtr, int );
extern void getFieldValue(FieldValuePtr, double, double, double, MagneticFieldPtr);
extern void getCompositeFieldValue(FieldValuePtr, double, double, double, MagneticFieldPtr, MagneticFieldPtr);
extern void setAlgorithm(Algorithm);
extern Algorithm getAlgorithm();
bool containsCartesian(MagneticFieldPtr, double, double, double);
bool containsCylindrical(MagneticFieldPtr, double, double);
extern void resetCell3D(Cell3DPtr, double, double, double);
extern void resetCell2D(Cell2DPtr, double, double);
extern void getCoordinateIndices(MagneticFieldPtr, double, double, double,
                          int *, int *, int *);


#endif /* magfield_h */
