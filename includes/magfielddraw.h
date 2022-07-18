//
// Created by David Heddle on 6/1/20.
//

#ifndef CMAG_MAGFIELDDRAW_H
#define CMAG_MAGFIELDDRAW_H

#include "magfield.h"

//external function prototypes
extern void createSVGImageFixedPhiDiff(char *, char *, double, MagneticFieldPtr, MagneticFieldPtr);
extern void createSVGImageFixedZDiff(char *, char *, double, MagneticFieldPtr, MagneticFieldPtr);
extern void createSVGImageFixedPhi(char *, char *, double, MagneticFieldPtr, MagneticFieldPtr);
extern void createSVGImageFixedZ(char *, char *, double, MagneticFieldPtr, MagneticFieldPtr);

#endif //CMAG_MAGFIELDDRAW_H
