#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include<mkl_lapacke.h>
#include "matOps.h"
#include "basis.h"
#include "vector.h"
#include "matrixElts.h"
#include "makeMats.h"
#include "diis.h"

// calculate nuclear repulsion energy for a set of nuclei
double calcTotalNucE(int nat, nuc *nuclei[]);

// calculate ground state energy
// used to evaluate convergence
double getE(int blen, double P[], double Hcore[], double F[]);

// uhf energy
double getUrE(int blen, double Pa[], double Pb[], double Pt[],
  double Fa[], double Fb[], double Hcore[]); 

// run scf calculation
void scf(int ur, int nat, int blen, int maxIters, double eTol, 
  nuc *nuclei[], orb *allOrbs[], int nA, int nB, int diisNum);
