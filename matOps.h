#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>
#include <math.h>

// copy a matrix
void copyMat(int N, int M, double *a, double *copy);

// print sqaure matrix
void printMat(int N, int M, bool rowf, double *a);

// switch form from row (col) to col (row)
void switchForm(int N, int M, bool rowf, double *a);

// transpose square matrix
void transpose(int N, int M, bool arowf, bool browf, double *a, double *b);

// add matrices
void addMats(int N, int M, bool rowf, double *a, double *b, double alpha, double beta, double *result);

// multiply matrices
void multMats(int N, int K, int M, bool rowf, double *a, double *b, double alpha, double beta, double *result);

// compute trace of matrix
double trace(int N, bool rowf, double *a);

// inner product of matrices
double dotMat(int N, bool rowf, double *a, double *b, double *work);

// prototype for matrix diagonalization
void diagMat(int N, int LDA, double *eVecs, double *eVals);
