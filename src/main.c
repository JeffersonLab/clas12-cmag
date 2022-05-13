//
//  main.c
//  cMag
//
//  Created by David Heddle on 5/22/20.
//  The main program runs the unit tests.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "magfield.h"
#include "magfieldio.h"
#include "munittest.h"
#include "magfieldutil.h"
#include "magfielddraw.h"

//the three fields we'll try to initialize
static MagneticFieldPtr symmetricTorus;
static MagneticFieldPtr fullTorus;
static MagneticFieldPtr solenoid;

//the one being used now
static MagneticFieldPtr currentField;

/**
 * Run the unit tests
 * return: NULL or the first failed test message
 */
static char *allTests() {
    fprintf(stdout, "\n\n***** Unit tests ****** \n");
    mu_run_test(gridUnitTest);
    mu_run_test(randomUnitTest);
    mu_run_test(conversionUnitTest);
    mu_run_test(binarySearchUnitTest);

    fprintf(stdout, "\n  [SYMMETRIC TORUS]");
    testFieldPtr = symmetricTorus;
    mu_run_test(compositeIndexUnitTest);
    mu_run_test(containsUnitTest);
    mu_run_test(nearestNeighborUnitTest);

    fprintf(stdout, "\n  [FULL  TORUS]");
    testFieldPtr = fullTorus;
    mu_run_test(compositeIndexUnitTest);
    mu_run_test(containsUnitTest);
    mu_run_test(nearestNeighborUnitTest);

    testFieldPtr = solenoid;
    fprintf(stdout, "\n  [SOLENOID]");
    mu_run_test(compositeIndexUnitTest);
    mu_run_test(containsUnitTest);
    mu_run_test(nearestNeighborUnitTest);

    fprintf(stdout, "\n ***** End of unit tests ******\n");
    return NULL;
}

/**
 * Generate some svg images of the maps
 */
static void svgImages() {

    //create some pictures
    char *title = (char *) malloc(255);
    double z = 375;
    double phi = 0;
    char *home = getenv("HOME");
    char *picPath = (char*) malloc(255);

    sprintf(title, "Symmetric Torus at phi = %-4.1f degrees", phi);
    sprintf(picPath, "%s/%s", home, "magfield.svg");
    createSVGImageFixedPhi(picPath, title, phi, symmetricTorus, solenoid);

    sprintf(title, "Symmetric Torus at z = %-4.1f cm", z);
    sprintf(picPath, "%s/%s", home, "symmFixedZ.svg");
    createSVGImageFixedZ(picPath, title, z, symmetricTorus, solenoid);


    sprintf(picPath, "%s/%s", home, "magfieldDiffFixedPhi.svg");
    sprintf(title, "Magnitude of difference between full and symmetric torus maps at phi = %-4.1f degrees", phi);
    createSVGImageFixedPhiDiff(picPath, title, phi, fullTorus, symmetricTorus);

    sprintf(picPath, "%s/%s", home, "magfieldDiffFixedZ.svg");
    sprintf(title, "Magnitude of difference between full and symmetric torus maps at z = %-4.1f cm", z);
    createSVGImageFixedZDiff(picPath, title, z, fullTorus, symmetricTorus);

    free(picPath);
    free(title);

}

/**
 * Get an print a field value
 * @param fieldValuePtr holds the value
 * @param x x coordinate in cm
 * @param y y coordinate in cm
 * @param z z coordinate in cm
 * @param field which field to use
 */
static void printFV(FieldValuePtr fieldValuePtr, double x, double y, double z, MagneticFieldPtr field) {

    getFieldValue(fieldValuePtr, x, y, z, field);

    printf("(%7.2f, %7.2f, %7.2f) cm  ", x, y, z);
   // printFieldValue(fieldValuePtr, stdout);
    printFieldValueFull(fieldValuePtr, "kG", stdout);
}

/**
 * Test some special field values
 */
static void testSpecialValues() {
    FieldValuePtr fieldValuePtr = (FieldValuePtr) malloc (sizeof(FieldValue));

    printFV(fieldValuePtr, -t30, 0, 380.0, currentField);
    printFV(fieldValuePtr, -46.5, 0, 380, currentField);
    printFV(fieldValuePtr, 32.8, 32.8, 380, currentField);
    printFV(fieldValuePtr, -32.8, -32.8, 380, currentField);
    printFV(fieldValuePtr, 23.2, 40.2, 380, currentField);



}

/**
 * Give the users some choices
 */
static void userCommands() {

    currentField = symmetricTorus;


    //for selections
    char *choice = (char *) malloc(80);
    bool done = false;

    //unit test results
    char *testResult;

    while (!done) {
        printf("\nOptions:");
        printf("\n\tq\tquit and exit program");
        printf("\n\tr\trun all tests");
        printf("\n\ts\tgenerate svg images");
        printf("\n\tt\ttest special values");

        printf("\n> ");
        scanf("%s", choice);

        printf("Selected [%c]\n", choice[0]);

        switch (choice[0]) {

            case 'q':
                exit(0);

            case 'r':
                //run the unit tests
                testResult = allTests();

                if (testResult != NULL) {
                    fprintf(stdout, "Unit test failed: [%s]\n", testResult);
                } else {
                    fprintf(stdout, "\nProgram ran successfully.\n");
                }
                break;

            case 's':
                svgImages();
                break;

            case 't':
                testSpecialValues();
                break;

            default:
                printf("Unrecognized command [%c]\n", choice[0]);
                break;
        }
    }
    free(choice);
}

/**
 * The main method of the test application.
 * @param argc the number of arguments
 * @param argv the command line argument. Only one is processed, the path
 * to the directory containing the magnetic fields. If that argument is missing,
 * it will look in $(HOME)/magfield.
 * @return 0 on successful completion, 1 if any error occurred.
 */
int main(int argc, const char * argv[]) {
    
    //paths to field maps
    char *solenoidPath = (char*) malloc(255);
    char *torusSymmetricPath = (char*) malloc(255);
    char *torusFullPath = (char*) malloc(255);

    //for testing will look for mag fields in either the first
    // command line argument, if provided, or in $HOME/magfield.
    //Here we build the paths
    const char *dataDir;

    if (argc > 1) {
        dataDir = argv[1];
        fprintf(stdout, "Using command line data directory: [%s]", dataDir);
    }
    else {
        dataDir = strcat(getenv("HOME"), "/magfield");
    }

    //build the full paths
    sprintf(solenoidPath, "%s/Symm_solenoid_r601_phi1_z1201_13June2018.dat", dataDir);
    sprintf(torusSymmetricPath, "%s/Symm_torus_r2501_phi16_z251_24Apr2018.dat", dataDir);
    sprintf(torusFullPath, "%s/Full_torus_r251_phi181_z251_03March2020.dat", dataDir);

    fprintf(stdout, "\nTesting the cMag library\n");

    //try to read the symmetric torus
    symmetricTorus = initializeTorus(torusSymmetricPath);
    if (symmetricTorus == NULL) {
        fprintf(stderr, "\ncMag ERROR failed to read symmetric torus map from [%s]\n", torusSymmetricPath);
        return 1;
    }

    //try to read the full torus
    fullTorus = initializeTorus(torusFullPath);
    if (fullTorus == NULL) {
        fprintf(stderr, "\ncMag ERROR failed to read full torus map from [%s]\n", torusFullPath);
        return 1;
    }


    //try to read the solenoid
    solenoid = initializeSolenoid(solenoidPath);
    if (solenoid == NULL) {
        fprintf(stderr, "\ncMag ERROR failed to read solenoid map from [%s]\n",
                solenoidPath);
        return 1;
    }

    //process user commands
    userCommands();

    //free the fieldmaps
    freeFieldMap(symmetricTorus);
    freeFieldMap(solenoid);
    freeFieldMap(fullTorus);

}




