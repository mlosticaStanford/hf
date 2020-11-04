//*****************************************************************//
// functions for forming matrices
//*****************************************************************//

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "matOps.h"
#include "vector.h"
#include "basis.h"
#include "matrixElts.h"



//*****************************************************************//
// compute compound indices
///////////////////////////////

int compoundIndices(int blen, int i, int j, int k, int l);

//*****************************************************************//



//*****************************************************************//
// overlap matrix S
///////////////////////////////

void makeS(int blen, orb *allOrbs[], double S[]);

//*****************************************************************//



//*****************************************************************//
// kinetic energy matrix T
///////////////////////////////

void makeT(int blen, orb *allOrbs[], double T[]);

//*****************************************************************//



//*****************************************************************//
// nuclear potential energy matrix Vn
///////////////////////////////

void makeVn(int blen, int nats, orb *allOrbs[], double T[], nuc *nuclei[]);

//*****************************************************************//



//*****************************************************************//
// core hamiltonian Hcore
///////////////////////////////

void makeHcore(int blen, int nats, orb *allOrbs[], 
  double T[], double Vn[], double Hcore[], nuc *nuclei[]);

//*****************************************************************//



//*****************************************************************//
// yeet eri
///////////////////////////////

void makeERI(int blen, orb *allOrbs[], double eri[]);

//*****************************************************************//



//*****************************************************************//
// make x matrix for canonical orthogonalization of s
///////////////////////////////

void makeX(int blen, double S[], double U[], double sp[], 
  double X[], double Xt[]);

//*****************************************************************//



//*****************************************************************//
// make g matrix
///////////////////////////////

void makeG(int blen, double P[], double eri[], double G[]); 

//*****************************************************************//



//*****************************************************************//
// make uhf g matrix
///////////////////////////////

void makeUrG(int blen, double Pa[], double Pb[], double Pt[],
  double eri[], double Ga[], double Gb[]);

//*****************************************************************//



//*****************************************************************//
// update density matrix
///////////////////////////////

void upP(int ur, int blen, int nElec, double P[], double C[]); 

//*****************************************************************//