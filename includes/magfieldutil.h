//
//  magfieldutil.h
//  cMag
//
//  Created by David Heddle on 5/22/20.
//  Copyright © 2020 David Heddle. All rights reserved.
//

#ifndef magfieldutil_h
#define magfieldutil_h

#include "magfield.h"

//used for comparing real numbers
extern const double TINY;


//external prototypes
extern void stringCopy(char **, const char *);
extern const char *fieldUnits(MagneticFieldPtr);
extern const char *lengthUnits(MagneticFieldPtr);
extern double fieldMagnitude(FieldValue *);
extern void printFieldSummary(MagneticFieldPtr, FILE *);
extern void printFieldValue(FieldValue *, FILE *);
extern void printFieldValueFull(FieldValue *, char *, FILE *);
extern MagneticFieldPtr createFieldMap(void);
extern void freeFieldMap(MagneticFieldPtr);
extern int randomInt(int, int);
extern double randomDouble(double, double);
extern char *randomUnitTest();
bool sameNumber(double, double);
void normalizeAngle(double *);
extern double min(double, double);
extern double max(double, double);

extern double toDegrees(double);
extern double toRadians(double);
extern void cartesianToCylindrical(const double, const double, double *, double *);
extern void cylindricalToCartesian(double *, double *, const double, const double);
extern char *conversionUnitTest();
int binarySearch(double *, int, int, double);
int descBinarySearch(double *, int, int, double);
extern char *binarySearchUnitTest();
extern int sign(double);
extern void sortArray(double *, int);
extern double relativePhi(double);
extern int getSector(double);

#endif /* magfieldutil_h */
