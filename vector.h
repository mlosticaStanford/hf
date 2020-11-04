// vector structure
// used for positions of nuclei
// based on code from Rabani group

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>

#define sqrt(x)  ((x) * (x))

// structure declarations

typedef struct vector_ {
  double x, y, z;
  double mag;
} vector;

// functions for vector structure in vector.c
vector retAddedVectors(vector vect1, vector vect2);
vector retSubtractedVectors(vector vect1, vector vect2);
vector scaleVector(vector vect, double scalar);
double retVectorMagnitude(vector vect);
double retDotProduct(vector vect1, vector vect2);
vector newVec(double x, double y, double z);

#endif
