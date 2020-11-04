#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include"vector.h"
#include"basis.h"


//*****************************************************************//
// parameters structure
///////////////////////////////
//
#ifndef CALCPARAMS_H 
#define CALCPARAMS_H

typedef struct calcParams_{
  
  // which basis to use
  char basis[1024];

  // uhf or rhf
  int uhf;

  // tolerance
  double tol;

  // uhf int/boolean
  // alpha and beta elecrons
  // max iterations
  // size of subspace for DIIS
  int aElecs, bElecs, maxIters, diisNum;
  
} calcParams;

#endif

double assignZ(char *s);

int getNat(char *fName);

void readMolecule(char *fName, nuc *nuclei[], int nat);

int getNorbs(char *fName);

//void readGTO(int orbLen, orb o);

void makeOrbs(int nat, int *orbCount, char *bName, nuc *nuclei[], orb *allOrbs[]);

void readParams(char *fName, calcParams *p);
