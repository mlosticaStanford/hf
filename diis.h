//*****************************************************************//
// functions running DIIS
//*****************************************************************//

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "matOps.h"
#include "vector.h"
#include "basis.h"
#include "matrixElts.h"



//*****************************************************************//
// compute error matrix
///////////////////////////////

void errVec(int blen, double e[], double F[], double P[], 
  double S[], double Xt[], double X[], double workA[], double workB[]);

//*****************************************************************//



//*****************************************************************//
// initialize bmatrix
///////////////////////////////


void initB(int numPrev, double B[]);


//*****************************************************************//



//*****************************************************************//
// initialize array of error or Fock matrices
///////////////////////////////


void initVecArray(int blen, int diisNum, double *vecArr[]);


//*****************************************************************//



//*****************************************************************//
// update an element of vector array
///////////////////////////////


void upVecArray(int blen, int diisNum, int iter, 
  double vecArr[], double newVec[]);


//*****************************************************************//



//*****************************************************************//
// make element of bmatrix
///////////////////////////////

double bElt(double ei[], double ej[], int blen, double work[]);

//*****************************************************************//



//*****************************************************************//
// update bmatrix
///////////////////////////////

void upB(int blen, int iter, int numPrev, double *errVecs[], 
  double work[], double B[]);

//*****************************************************************//



//*****************************************************************//
// extrapolate Fock matrix
///////////////////////////////

double newFock(int blen, int diisNum, double coeffs[], 
  double fockArr[], double newF[]);

//*****************************************************************//



//*****************************************************************//
// run DIIS
///////////////////////////////

void runDIIS(int blen, int diisNum, int iter,
  double *fockArr[], double *errVecs[],
  double B[], double coeffs[], 
  double F[], double upErr[], 
  double P[], double S[], double Xt[], double X[],
  double workA[], double workB[], double bCopy[], double errNorm);

//*****************************************************************//
