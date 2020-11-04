// operations for matrices
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>
#include <math.h>
#include "matOps.h"

// copy a matrix
// the copy will be the same form (row or col)
// as the input
void copyMat(int N, int M, double *a, double *copy) {
  // N, rows of matrix
  // M, cols of matrix
  // *a, pointer to 1D array for matrix
  // *copy, pointer to 1D array to store the copy of a

  // i,   1D index
  int i;

  // copy elements
  for (i = 0; i < (N*M); i++) {
    copy[i] = a[i];
  }
}

// print sqaure matrix
void printMat(int N, int M, bool rowf, double *a) {
  // N, rows of matrix
  // M, cols of matrix
  // rowf, is matrix in row form?
  // *a, pointer to 1D array for matrix

  // i,j, 2D indices
  // k,   1D index
  int i,j,k;
  bool doPrint = true;

  // is matrix too large to print?
  if (N > 10 || M > 10 || N*M > 100) {
    printf("This matrix is probably too large to print :| \n");
    doPrint = false;
  }

  // print matrix if it's not too big
  if (rowf && doPrint) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        k = (i*M) + j;
        printf("%.4f   ", a[k]);
      }
      printf("\n");
    }
  }

  else if (!rowf || doPrint) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        k = i + (j*N);
        printf("%.4f   ", a[k]);
      }
      printf("\n");
    }
  }
}

// switch form from row (col) to col (row)
void switchForm(int N, int M, bool rowf, double *a) {
  // N, rows of matrix a
  // M, cols of matrix a
  // rowf, is matrix a in row form?
  // *a, pointer to 1D array for matrix

  // initialize work array
  double *w;
  w = (double *)malloc(sizeof(double)*N*M);

  // initialize indices
  int i,j,ak,wk,k;

  // if we want to have a go from row to col
  if (rowf) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        ak = (i*M) + j;
        wk = i + (j*N);
        w[wk] = a[ak];
      }
    }
  }

  // if we want to have a go from col to row
  else {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        ak = i + (j*N);
        wk = (i*M) + j;
        w[wk] = a[ak];
      }
    }
  }

  // now reassign values of a
  for (k = 0; k < (N*M); k++) {
    a[k] = w[k];
  }

  // free memory
  free(w);
}


// transpose square matrix
void transpose(int N, int M, bool arowf, bool browf, double *a, double *b) {
  // N, rows of matrix a
  // M, cols of matrix a
  // arowf, is matrix a in row form?
  // browf, is matrix b in row form?
  // *a, pointer to 1D array for matrix
  // *b, pointer to 1D array to store result

  // i,j, 2D indices
  // ak, bk   1D index
  int i,j,ak,bk;

  // if a and b are in row form
  if (arowf && browf) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        ak = (i*M) + j;
        bk = (j*N) + i;
        b[bk] = a[ak];
      }
    }
  }

  // elif a is in col form and b is in row form
  else if (!arowf && browf) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        ak = i + (j*N);
        bk = (j*N) + i;
        b[bk] = a[ak];
      }
    }
  }

  // elif a is in row form and b is in col form
  else if (arowf && !browf) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        ak = (i*M) + j;
        bk = j + (i*M);
        b[bk] = a[ak];
      }
    }
  }

  // else a and b are in col form
  else {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        ak = i + (j*N);
        bk = j + (i*M);
        b[bk] = a[ak];
      }
    }
  }
}

// add matrices
void addMats(int N, int M, bool rowf, double *a, double *b, double alpha, double beta, double *result) {
  // N, rows of matrices
  // M, cols of matrices
  // rowf, are matrices in row form?
  // !!! all matrices must be in the same form !!!
  // *a, *b, *result, pointers to 1D arrays of matrices
  // !!! result must be initialized first !!!
  // alpha, beta, scalars for a and b respectively

  // initialize indices
  int i,j,k;

  // if matrices are in row form
  if (rowf) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        k = (i*M) + j;
        result[k] = (alpha * a[k]) + (beta * b[k]);
      }
    }
  }

  // if matrices are in col form
  else {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        k = i + (j*N);
        result[k] = (alpha * a[k]) + (beta * b[k]);
      }
    }
  }
}

// multiply matrices
void multMats(int N, int L, int M, bool rowf, double *a, double *b, double alpha, double beta, double *result) {
  // N, rows of matrix a
  // L, cols (rows) of matrix a (b)
  // M, cols of matrix b
  // rowf, are matrices in row form?
  // !!! all matrices must be in the same form !!!
  // *a, *b, *result, pointers to 1D arrays of matrices
  // !!! result must be initialized first !!!
  // alpha, beta, scalars for a and b respectively
  
  // initialize indices
  // and total to be retured
  int i,j,k,aik,bkj,cij;
  double cElt;

  // if matrices are in row form
  if (rowf) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        cElt = 0.0;
        cij = (i*M) + j;
        for (k = 0; k < L; k++) {
          aik = (i*L) + k;
          bkj = (k*M) + j;
          //printf("%f is a[aik] and %f is b[bkj] \n", a[aik], b[bkj]);
          cElt += (alpha * a[aik]) * (beta * b[bkj]); 
        }
        result[cij] = cElt;
      }
    }
  }

  // if matrices started in col form
  else {
    for (i = 0; i < N; i++) {
      for (j = 0; j < M; j++) {
        cElt = 0.0;
        cij = i + (j*N);
        for (k = 0; k < L; k++) {
          aik = i + (k*N);
          bkj = k + (j*L);
          //printf("%f is a[aik] and %f is b[bkj] \n", a[aik], b[bkj]);
          cElt += (alpha * a[aik]) * (beta * b[bkj]); 
        }
        result[cij] = cElt;
      }
    }
  }
}

// compute trace of matrix
double trace(int N, bool rowf, double *a) {
  // N, rows of matrice
  // rowf, is matrix in row form?
  // *a, pointers to 1D array of matrix

  // initialize index 
  // and result to be returned
  int i,j,k;
  double total = 0.0;

  // if row form
  if (rowf) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        
        // compound index
        k = (i*N) + j;
        
        // sum over diagonals
        if (i == j) {
        total += a[k];
        }
      }
    }
  }

  // if col form
  else {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        
        // compound index
        k = i + (j*N);
        
        // sum over diagonals
        if (i == j) {
        total += a[k];
        }
      }
    }
  }

  return total;
}

// inner product of matrices
double dotMat(int N, bool rowf, double *a, double *b, double *work) {
  // N, rows of matrices
  // rowf, are matrices in row form?
  // !!! all matrices must be in the same form !!!
  // *a, *b, *work, pointers to 1D arrays of matrices
  // !!! work must be initialized first !!!

  // initialize double for trace 
  double t;


  // if in row form
  if (rowf) {
  // multiply matrices
	  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                                          N, N, N,
                                          1.0, a, N, b, N,
                                          0.0, work, N);
    t = trace(N, true, work);
  }

  // if in col form
  else {
  // multiply matrices
	  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                                          N, N, N,
                                          1.0, a, N, b, N,
                                          0.0, work, N);
    t = trace(N, false, work); 
  }

  return t;
}

// diagoanlize a matrix a
// !!! the output eigenvector array will be in column form !!!
void diagMat(int N, int LDA, double eVecs[], double eVals[]) {
  int n = N, lda = LDA, info, lwork;
  double wkopt;
  double* work;

  // MKL stuff
  lwork = -1;
  dsyev( "Vectors", "Upper", &n, eVecs, &lda, eVals, &wkopt, &lwork, &info );
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  dsyev( "Vectors", "Upper", &n, eVecs, &lda, eVals, work, &lwork, &info );


  /* Check for convergence */
  if( info > 0 ) {
    printf( "The algorithm failed to compute eigenvalues.\n" );
    exit( 1 );
  }

  free( (void*)work );
}
