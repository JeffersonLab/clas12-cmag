//
// Created by David Heddle on 5/31/20.
//

#include "magfieldio.h"
#include "magfieldutil.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <arpa/inet.h>

//do we have to swap bytes?
//since the fields were produced by Java which uses
//new format (BigEndian) we probably will have to swap.
static bool swapBytes = false;

//local prototypes
static FieldMapHeaderPtr readMapHeader(FILE *);
static MagneticFieldPtr readField(const char *);
static long getFileSize(FILE*);
static void swap32(char*, int);
static char* getCreationDate(MagneticFieldPtr);
static void computeFieldMetrics(MagneticFieldPtr);

/**
 * Initialize the torus field.
 * @param torusPath a path to a torus field map. If you want to use environment variables, pass NULL
 * in this parameter, in which case the code make two attempts two attempts at finding the field. The first
 * will be to try the a path specified by the COAT_MAGFIELD_TORUSMAP environment variable. If that
 * fails, it will then check TORUSMAP.
 * @return a valid field pointer on success, NULL on failure.
 */
MagneticFieldPtr initializeTorus(const char *torusPath) {
    if (torusPath == NULL) {
        torusPath = getenv("COAT_MAGFIELD_TORUSMAP");
        if (torusPath == NULL) {
            torusPath = getenv("TORUSMAP");
        }
    }

    if (torusPath == NULL) {
        fprintf(stderr, "\ncMag ERROR null torus path even after trying environment variables.\n");
        return NULL;
    }
    return readField(torusPath);
}

/**
 * Initialize the solenoid field.
 * @param solenoidPath a path to a solenoid field map. If you want to use environment variables, pass NULL
 * in this parameter, in which case the code make two attempts two attempts at finding the field. The first
 * will be to try the a path specified by the COAT_MAGFIELD_SOLENOIDMAP environment variable. If that
 * fails, it will then check SOLENOIDMAP.
 * @return a valid field pointer on success, NULL on failure.
 */
MagneticFieldPtr initializeSolenoid(const char *solenoidPath) {
    if (solenoidPath == NULL) {
        solenoidPath = getenv("COAT_MAGFIELD_SOLENOIDMAP");
        if (solenoidPath == NULL) {
            solenoidPath = getenv("SOLENOIDMAP");
        }
    }

    if (solenoidPath == NULL) {
        fprintf(stderr, "\ncMag ERROR null solenoid path even after trying environment variables.\n");
        return NULL;
    }
    return readField(solenoidPath);
}


/**
 * Read a binary field map at the given location.
 * @param path the full path to a field map file.
 * @return a valid field pointer on success, NULL on failure.
 */
static MagneticFieldPtr readField(const char *path) {

    debugPrint("\nAttempting to read field map from [%s]\n", path);
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        fprintf(stderr, "\ncMag ERROR could not read field map file: [%s]\n", path);
        return NULL;
    }

    //get the header
    FieldMapHeaderPtr headerPtr = readMapHeader(file);
    if (headerPtr == NULL) {
        fclose(file);
        fprintf(stderr, "\ncMag ERROR could not read field map header from: [%s]\n", path);
        return NULL;
    }

    MagneticFieldPtr fieldPtr = createFieldMap();

    //copy the path and name
    stringCopy(&(fieldPtr->path), path);

    fieldPtr->headerPtr = headerPtr;
    fieldPtr->numValues = headerPtr->nq1 * headerPtr->nq2 * headerPtr->nq3;
    fieldPtr->creationDate = getCreationDate(fieldPtr);

    //malloc the data array
    fieldPtr->fieldValues = malloc(fieldPtr->numValues * sizeof(FieldValue));

    //did we have enough memory?
    if (fieldPtr->fieldValues == NULL) {
        fprintf(stderr, "\ncMag ERROR out of memory when allocating space for field map.\n");
        fclose(file);
        return NULL;
    }

    //now we can read the field. Reading the header should have left
    //the file pointer positioned at the right spot.
    fread(fieldPtr->fieldValues, sizeof(FieldValue), fieldPtr->numValues, file);
    fclose(file);

    //swap?
    if (swapBytes) {
        for (int i = 0; i < fieldPtr->numValues; i++) {
            swap32((char*) (fieldPtr->fieldValues + i), 3);
        }
    }

    //create the coordinate grids
    //CLAS fields always have cylindrical grids
    //with q1 = phi, q2 = rho and q3 = z
    fieldPtr->phiGridPtr = createGrid("phi", headerPtr->q1min,
                                      headerPtr->q1max, headerPtr->nq1);
    fieldPtr->rhoGridPtr = createGrid("rho", headerPtr->q2min,
                                      headerPtr->q2max, headerPtr->nq2);
    fieldPtr->zGridPtr = createGrid("z", headerPtr->q3min,
                                    headerPtr->q3max, headerPtr->nq3);

    //this is useful to cache for indexing purposes
    fieldPtr->N23 = headerPtr->nq2 * headerPtr->nq3;

    //solenoid files have nq1 (nPhi) = 1 and are symmetric (no phi dependence apart from rotation)
    //torus symmetric if phi max = 30.
    fieldPtr->symmetric = false;

    if (headerPtr->nq1 < 2) {
        fieldPtr->type = SOLENOID;
        fieldPtr->symmetric = true;
        fieldPtr->cell3DPtr = NULL;
        createCell2D(fieldPtr);
    }
    else {
        fieldPtr->type = TORUS;
        if ((headerPtr->q1max - headerPtr->q1min) < 31) {
            fieldPtr->symmetric = true;
        }
        createCell3D(fieldPtr);
        fieldPtr->cell2DPtr = NULL;
    }


    //compute some metrics
    computeFieldMetrics(fieldPtr);

    printFieldSummary(fieldPtr, stdout);
    return fieldPtr;
}

/**
 * Create a 3D cell, which is used by the torus.
 * Note that nothing is
 * returned, the field's 2D cell pointer is made to point at the cell,
 * and the cell is given a reference to the field.
 * @param fieldPtr a pointer to the solenoid field.
 */
void createCell3D(MagneticFieldPtr fieldPtr) {
    Cell3DPtr cell3DPtr = (Cell3DPtr) malloc(sizeof(Cell3D));
    cell3DPtr->phiMin = INFINITY;
    cell3DPtr->phiMax = -INFINITY;
    cell3DPtr->rhoMin = INFINITY;
    cell3DPtr->rhoMax = -INFINITY;
    cell3DPtr->zMin = INFINITY;
    cell3DPtr->zMax = -INFINITY;
    cell3DPtr->fieldPtr = fieldPtr;
    fieldPtr->cell3DPtr = cell3DPtr;
}

/**
 * Create a 2D cell, which is used by the solenoid, since the lack
 * of phi dependence renders the solenoidal field effectively 2D.
 * Note that nothing is
 * returned, the field's 2D cell pointer is made to point at the cell,
 * and the cell is given a reference to the field.
 * @param fieldPtr a pointer to the solenoid field.
 */
void createCell2D(MagneticFieldPtr fieldPtr) {
    Cell2DPtr cell2DPtr = (Cell2DPtr) malloc(sizeof(Cell2D));
    cell2DPtr->rhoMin = INFINITY;
    cell2DPtr->rhoMax = -INFINITY;
    cell2DPtr->zMin = INFINITY;
    cell2DPtr->zMax = -INFINITY;
    cell2DPtr->fieldPtr = fieldPtr;
    fieldPtr->cell2DPtr = cell2DPtr;
}

/**
 * Free the memory associated with a 3D cell.
 * @param cell3DPtr a pointer to the cell.
 */
void freeCell3D(Cell3DPtr cell3DPtr) {
    free(cell3DPtr);
}

/**
 * Free the memory associated with a 2D cell.
 * @param cell2DPtr a pointer to the cell.
 */
void freeCell2D(Cell2DPtr cell2DPtr) {
    free(cell2DPtr);
}

/**
 * Read the 80 byte field map header.
 * @param fd the file descriptor.
 * @return a valid pointer to a field map header, or NULL upon failure.
 */
static FieldMapHeaderPtr readMapHeader(FILE *fd) {

    //create space for the header
    FieldMapHeaderPtr headerPtr = (FieldMapHeaderPtr) malloc(
            sizeof(FieldMapHeader));

    //get the actual file size in bytes
    long actualFileSize = getFileSize(fd);
    debugPrint("Actual file size: %ld bytes\n", actualFileSize);

    //get the magic word and see if byteswap required
    fread(&(headerPtr->magicWord), sizeof(unsigned int), 1, fd);
    swapBytes = (headerPtr->magicWord != MAGICWORD);

    debugPrint("byteswap required: %s\n", swapBytes ? "yes" : "no");

    if (swapBytes) {
        headerPtr->magicWord = htonl(headerPtr->magicWord);
    }

    //if still not a match, it is fatal
    if (headerPtr->magicWord != 0xced) {
        fprintf(stderr,
                "\ncMag ERROR Magic word doesn't match, even after a byte swap.\n");
        free(headerPtr);
        return NULL;
    }

    //now that we know if we have to swap, let's read the header all at once
    rewind(fd);
    fread(headerPtr, sizeof(FieldMapHeader), 1, fd);

    if (swapBytes) {
        swap32((char*) headerPtr, sizeof(FieldMapHeader) / 4);
    }

    //if debugging, print header
    debugPrint("\nHEADER:\nmagic word: \"%03x\"\n", headerPtr->magicWord);
    debugPrint("grid CS: %s\n", csLabels[headerPtr->gridCS]);
    debugPrint("field CS: %s\n", csLabels[headerPtr->fieldCS]);

    debugPrint("length units: %s\n", lengthUnitLabels[headerPtr->lengthUnits]);
    debugPrint("angular units: %s\n", angleUnitLabels[headerPtr->angleUnits]);
    debugPrint("field units: %s\n", fieldUnitLabels[headerPtr->fieldUnits]);
    debugPrint("q1Min: %-5.2f\n", headerPtr->q1min);
    debugPrint("q1Max: %-5.2f\n", headerPtr->q1max);
    debugPrint("NumQ2: %d\n", headerPtr->nq1);
    debugPrint("q2Min: %-5.2f\n", headerPtr->q2min);
    debugPrint("q2Max: %-5.2f\n", headerPtr->q2max);
    debugPrint("NumQ2: %d\n", headerPtr->nq2);
    debugPrint("q3Min: %-5.2f\n", headerPtr->q3min);
    debugPrint("q3Max: %-5.2f\n", headerPtr->q3max);
    debugPrint("NumQ3: %d\n", headerPtr->nq3);

    //get the number of field values and the computed file size
    int numFieldValues = headerPtr->nq1 * headerPtr->nq2 * headerPtr->nq3;
    long computedFileSize = sizeof(FieldMapHeader) + 4 * 3 * numFieldValues;

    debugPrint("Computed file size: %ld bytes\n", computedFileSize);
    if (actualFileSize != computedFileSize) {
        fprintf(stderr,
                "\ncMag ERROR computed file size and actual file size do not match.\n");
        free(headerPtr);
        return NULL;
    }

    return headerPtr;
}

/**
 * Compute some diagnostic metrics for this field.
 * @param fieldPtr the field map pointer.
 */
static void computeFieldMetrics(MagneticFieldPtr fieldPtr) {

    FieldMetricsPtr metrics = fieldPtr->metricsPtr;
    metrics->maxFieldIndex = 0;
    metrics->maxFieldMagnitude = 0;
    metrics->avgFieldMagnitude = 0;

    for (unsigned int i = 0; i < fieldPtr->numValues; i++) {
        FieldValuePtr fieldValuePtr = fieldPtr->fieldValues+i;

        double magnitude = fieldMagnitude(fieldValuePtr);
        if (magnitude > metrics->maxFieldMagnitude) {
            metrics->maxFieldMagnitude = magnitude;
            metrics->maxFieldIndex = i;
        }

        metrics->avgFieldMagnitude += magnitude;
    }

    metrics->avgFieldMagnitude /= fieldPtr->numValues;
}


/**
 * Get the creation date of a field map.
 * @param fieldPtr a pointer to the field map.
 * @return A string representation of the date and the the field map
 * was created from the engineering data.
 */
static char *getCreationDate(MagneticFieldPtr fieldPtr) {

    int high = fieldPtr->headerPtr->cdHigh;
    int low = fieldPtr->headerPtr->cdLow;

    //the divide by 1000 is because Java creation time (which was used) is in nS
    long dlow = low & 0x00000000ffffffffL;
    time_t utime = (((long) high << 32) | (dlow & 0xffffffffL)) / 1000;

    return ctime(&utime);

}

/**
 * Get the file size in bytes.
 * @param fd the file descriptor.
 * @return the file size in bytes.
 */
static long getFileSize(FILE *fd) {
    long size;
    fseek(fd, 0L, SEEK_END);
    size = ftell(fd);
    rewind(fd);
    return size;
}

/**
 * Byte swap a block of 32-bit entities
 * @param ptr a pointer to the block.
 * @param num32 the number of entities.
 */
static void swap32(char *ptr, int num32) {
    int i, j, k;
    char temp[4];

    for (i = 0, j = 0; i < num32; i++, j += 4) {
        for (k = 0; k < 4; k++) {
            temp[k] = ptr[j + k];
        }
        for (k = 0; k < 4; k++) {
            ptr[j + k] = temp[3 - k];
        }
    }
}
