//*****************************************************************//
// basis strctures for s gaussian-type orbitals
//*****************************************************************//

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"



//*****************************************************************//
// primitive gto
///////////////////////////////
//
#ifndef PRIM_H
#define PRIM_H

typedef struct prim_{
  
  // ctr, where gto is ctrd
  vector ctr;

  // ex, exponent
  // co, coefficient
  // pN, primitive normalization const
  double ex, co, pN;
  
} prim;


///////////////////////////////
// prim funcs
///////////////////////////////

double normPrim(prim p);

#endif

//*****************************************************************//



//*****************************************************************//
// gto
///////////////////////////////
//
#ifndef ORB_H
#define ORB_H 

typedef struct orb_ {

  // orbLen, number of prims in orb
  int orbLen;

  // oN, orb normalization const
  double oN;

  // prims, array of primitive gtos
  prim orbPrims[];
  
} orb;


///////////////////////////////
// orb funcs
///////////////////////////////
double normOrbPrims(prim o[], int orbLen);

#endif

//*****************************************************************//



//*****************************************************************//
// nucleus strut
///////////////////////////////

#ifndef NUC_H
#define NUC_H 

typedef struct nuc_ {

  // atmSym, atomic symbol
  char atmSym[1024];

  // orbLen, number of prims in orb
  vector coords;

  // z, charge of nucleus 
  double z;

} nuc;

#endif

//*****************************************************************//

void normAllOrbs(int orbCount, orb *allOrbs[]);
