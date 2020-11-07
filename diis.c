//*****************************************************************//
// functions for running DIIS
//*****************************************************************//

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "diis.h"
#include "matOps.h"
#include "vector.h"
#include "basis.h"
#include "matrixElts.h"



//*****************************************************************//
// compute error matrix
///////////////////////////////

void errVec(int blen, double e[], double F[], double P[], 
  double S[], double X[], double workA[], double workB[]) {

  // compute XtFPS, store it in workB
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, X, blen, F, blen,
                                        0.0, workA, blen);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, workA, blen, P, blen,
                                        0.0, workB, blen);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, workB, blen, S, blen,
                                        0.0, workA, blen);

  // compute SPFX, store it in e
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, S, blen, P, blen,
                                        0.0, workB, blen);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, workB, blen, F, blen,
                                        0.0, e, blen);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, e, blen, X, blen,
                                        0.0, workB, blen);

  // compute FPS - SPF = [comm], store in e
  addMats(blen, blen, false, workA, workB, 1.0, -1.0, e);


  //printf("current err Vec\n");
  //printMat(blen,blen, false,e);
  //printf("\n\n");

} 

//*****************************************************************//



//*****************************************************************//
// is max value in err vec below threshold?
///////////////////////////////

int errTol(int blen, double e[], double tol) {
  
  double maxVal = fabs(e[0]);
  int i;

  for (i = 0; i < blen*blen; i++) {
    if (fabs(e[i]) > maxVal) {
      maxVal = fabs(e[i]);
    }
  }

  if (maxVal < tol) {
    return 1;
  }
  else {
    return 0;
  }
} 

//*****************************************************************//



//*****************************************************************//
// initialize bmatrix
///////////////////////////////


void initB(int numPrev, double B[]) {
  // numPrev, number of error vectors used
  // B, initialized b matrix. Note that it has dim (numPrev + 1)**2
  
  // form B iteratively
  int i,j,k,l;
  for (i = 0; i < numPrev + 1; i++) {
    for (j = 0; j <= i; j++) {

      // 1D index
      k = i + (j*(numPrev +1 ));
      l = (i*(numPrev + 1)) + j;

      // work in the area of matrix that's not
      // the bottom row
      // or the rightmost column
      if (i < numPrev && j < numPrev) {
        B[k] = 0.0; 
        
        // assign off-diag elts
        if (k != l) {
          B[l] = 0.0; 
        }
      } 

      // work in bottom row/rightmost column
      else if (i == numPrev && j != i) {
        B[k] = 1.0; 
        B[l] = 1.0; 

      // bottom right elt
      } 
      else {
        B[k] = 0.0;
      }
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// initialize array of error or Fock matrices
///////////////////////////////


void initVecArray(int blen, int diisNum, double *vecArr[]) {
  // blen, length of basis
  // diisNum, number of previous Fock matrices to use
  // vecArr, array of matrices

  // indices
  int i, arrIndex;

  for (arrIndex = 0; arrIndex < diisNum; arrIndex++) {

    // allocate memory for matrix
    vecArr[arrIndex] = (double *)malloc(sizeof(double)*blen*blen);

    // initialize matrix as zeros
    for (i = 0; i < blen*blen; i++) {
      vecArr[arrIndex][i] = 0.0;
    }
  }
}


//*****************************************************************//



//*****************************************************************//
// update an element of vector array
///////////////////////////////


void upVecArray(int blen, int diisNum, int iter, 
  double *vecArr[], double newVec[]) {
  // blen, length of basis
  // diisNum, number of previous Fock matrices to use
  // iter, iteration of SCF
  // vecArr, array of matrices
  // newVec, new vector to be added to the array

  // indices
  int arrIndex;

  // only update a single element of the array
  arrIndex = iter % diisNum;

  // copy the new vector into array
  copyMat(blen, blen, newVec, vecArr[arrIndex]);

  int i;
  for (i = 0; i < diisNum; i++) {
    printf("vector %i is: \n", i);
    printMat(blen, blen, false, vecArr[i]);
  }
}


//*****************************************************************//



//*****************************************************************//
// make element of bmatrix
///////////////////////////////

double bElt(double ei[], double ej[], int blen, double work[]) {
  // ei,ej, respectively ith and jth error vectors
  // blen, length of basis (dim of error vecs)
  // work, work matrix
  
  return dotMat(blen, false, ei, ej, work);
}

//*****************************************************************//



//*****************************************************************//
// update bmatrix
///////////////////////////////

void upB(int blen, int iter, int numPrev, double *errVecs[], 
  double work[], double B[]) {
  // blen, length of basis
  // iter, index of iteration
  // numPrev, number of error vectors used
  // errVecs, array of error matrices
  // work, work matrix
  // B, initialized b matrix. Note that it has dim (numPrev + 1)**2

  // iterate over each error vector
  int i,j,k,l,bIndex;

  // maxErr to get index of worst error
  // epsilon for regularization
  double maxErr, eps = 0.01;

  // only update using the new vector
  bIndex = iter % numPrev;

  // update using the worst error 
  if (iter >= numPrev * 2) {
    maxErr = B[0];
    bIndex = 0;
    for (i = 1; i < numPrev; i++) {
      k = i + (i*(numPrev +1 ));
      if (B[k] > maxErr) {
        maxErr = B[k];
        bIndex = i;
      }
    }
  }
  else {
    bIndex = iter % numPrev;
  }

  printf("bIndex for upB = %i\n", bIndex);
  

  for (i = 0; i < numPrev + 1; i++) {
    for (j = 0; j <= i; j++) {

      // 1D index
      k = i + (j*(numPrev +1 ));
      l = (i*(numPrev + 1)) + j;

      // work in the area of matrix that's not
      // the bottom row
      // or the rightmost column
      if (i < numPrev && j < numPrev && (i == bIndex || j == bIndex)) {
        B[k] = bElt(errVecs[i], errVecs[j], blen, work);

        
        // assign off-diag elts
        if (k != l) {
          B[l] = B[k]; 
        }
        // regularization
        else {
          //B[k] += eps*exp((-1.0 * B[k]) / eps);
          B[k] += eps;
        }
      } 

      // work in bottom row/rightmost column
      else if (i == numPrev && j != i) {
        B[k] = 1.0; 
        B[l] = 1.0; 

      // bottom right elt
      } 
      else if (i == numPrev && j == numPrev){
        B[k] = 0.0;
        //printf("made B[%i] = 0.0\n", k, B[k]);
      }
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// extrapolate Fock matrix
///////////////////////////////

double newFock(int blen, int diisNum, double coeffs[], 
  double *fockArr[], double newF[]) {
  // blen, length of basis
  // diisNum, number of previous Fock matrices to use
  // coeffs, coefficients of old matrices, from solving B
  // fockArr, array of Fock matrices
  // newF, extrapolated Fock matrix

  // extrapolate Fock matrix
  // weight old Fmats with coeffs
  int i;

  // zero out new Fock matrix
  for (i = 0; i < pow(blen,2); i++) {
    newF[i] = 0.0;
  }
  for (i = 0; i < diisNum; i++) {
    printf("fockArr[%i] is \n", i);
    printMat(blen,blen,false,fockArr[i]);
    //printf("coeffs[%i] = %f \n", i, coeffs[i]);
    addMats(blen, blen, false, newF, 
            fockArr[i], 1.0, coeffs[i], newF);
  }
}

//*****************************************************************//



//*****************************************************************//
// run DIIS
///////////////////////////////

void runDIIS(int blen, int diisNum, int iter,
  double *fockArr[], double *errVecs[],
  double B[], double coeffs[], 
  double F[], double upErr[], 
  double P[], double S[], double X[],
  double workA[], double workB[], double bCopy[], double errNorm) {
  // blen, length of basis
  // diisNum, number of previous Fock matrices to use
  // iter, iteration of SCF
  // fockArr, array of Fock matrices
  // errVecs, array of error matrices
  // B, matrix formed from dotting error matrices
  // coeffs, coefficients of old matrices, from solving B
  // upF, Fock matrix used to update fockArr
  // upErr, error matrix used to update errVecs
  // newF, initialized Fock matrix to be extrapolated
  // P, density matrix
  // S, overlap matrix
  // workA, workB, work matrices
  

  // only update a single vector
  int bIndex, i, k;
  double maxErr;

  bIndex = iter % diisNum;

  //// update using the worst error 
  if (iter >= diisNum * 2) {
    maxErr = B[0];
    bIndex = 0;
    for (i = 0; i < diisNum; i++) {
      k = i + (i*(diisNum +1 ));
      if (B[k] > maxErr) {
        maxErr = B[k];
        bIndex = i;
      }
    }
  }
  else {
    bIndex = iter % diisNum;
  }

  // update fock vectors
  copyMat(blen, blen, F, fockArr[bIndex]);

  // compute new error vector
  errVec(blen, upErr, F, P, S, X, workA, workB);
  copyMat(blen, blen, upErr, errVecs[bIndex]);

  printf("Old B = \n");
  printMat(diisNum+1, diisNum + 1, false, B);

  // update B matrix
  upB(blen, iter, diisNum, errVecs, 
      workA, B);

  errNorm = dotMat(blen, false, errVecs[bIndex], errVecs[bIndex], workA);
  printf("errNorm = %f\n", errNorm);

  if (iter >= diisNum) {

    printf("current Fock matrix: \n");
    printMat(blen,blen,false,fockArr[bIndex]);

    printf("current density matrix: \n");
    printMat(blen,blen,false,P);

    printf("current error matrix: \n");
    printMat(blen,blen,false,errVecs[bIndex]);

    // copy B matrix because it will be altered by dgesv
    copyMat(diisNum+1,diisNum+1,B,bCopy);

    printf("new B = \n");
    printMat(diisNum+1,diisNum+1,false,B);

    // initialize coefficients
    // note that it's a (diisNum + 1) x 1 vector
    int i, info;
    int ipiv[diisNum + 1];
    for (i = 0; i < diisNum; i++) {
      coeffs[i] = 0.0;
    }
    coeffs[diisNum] = 1.0;

    // solve B matrix
    info = LAPACKE_dgesv( LAPACK_COL_MAJOR, diisNum + 1, 1, bCopy, 
          diisNum + 1, ipiv, coeffs, diisNum + 1);

    printf("solved coeffs matrix is \n");
    printMat(diisNum + 1, 1, false,coeffs);

    // compute the extrapolated Fock Matrix
    newFock(blen, diisNum, coeffs, fockArr, F);

    // copy b matrix back
    //copyMat(diisNum+1,diisNum+1,bCopy,B);

    printf("F after DIIS\n");
    printMat(blen,blen, false, F);
    printf("\n\n");

  }
}

//*****************************************************************//
