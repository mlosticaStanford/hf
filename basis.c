//*****************************************************************//
// functions for prim and orb strutures
//*****************************************************************//

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include "basis.h"
#include "vector.h"
#define pi 3.14159265358979323846



//*****************************************************************//
// prim funcs
///////////////////////////////

double normPrim(prim p){

  // p, prim instance

  return pow((pi/(2.0*p.ex)),-0.75);
}

//*****************************************************************//



//*****************************************************************//
// orb funcs
///////////////////////////////

double normOrbPrims(prim o[], int orbLen){

  // o, array of prims
  // orbLen, number of prims
  // !!! note that each prim must be normalized
  // before calling normOrbPrims!!!

  // initialize indices and total to be returned
  int p, q;
  double total = 0.0;

  // initialize vars for calculation
  vector vecDiff;
  double distSq, pre1, pre2, exponent;

  // normalize each prim
  for(p = 0; p < orbLen; p++) {
    o[p].pN = normPrim(o[p]);
  }

  // compute oN
  for( p = 0 ; p < orbLen ; p++){
    for( q = 0 ; q < orbLen ; q++){
      
      // get square distances btw prim orbs
      vecDiff = retSubtractedVectors(o[p].ctr, o[q].ctr);
      distSq = retDotProduct(vecDiff, vecDiff);

      // compute factors
      pre1 = o[p].co* o[q].co* o[p].pN * o[q].pN;
      pre2 = pow((pi/(o[p].ex + o[q].ex)),1.5);
      exponent = exp(-1.0*((o[p].ex*o[q].ex)/(o[p].ex+o[q].ex))*distSq);

      // increment total
      total += (pre1 * pre2 * exponent);
    }
  }

  return(pow(total,-0.5));
}

//*****************************************************************//


void normAllOrbs(int orbCount, orb *allOrbs[]) {
  int i;

  for (i = 0; i < orbCount; i++) {
    allOrbs[i]->oN = normOrbPrims(allOrbs[i]->orbPrims, allOrbs[i]->orbLen);
  }
}
