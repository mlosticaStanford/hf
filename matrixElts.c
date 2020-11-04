//*****************************************************************//
// functions for computing matrix elements
//*****************************************************************//

#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "basis.h"
#include "matrixElts.h"
#include "matOps.h"
#define pi 3.14159265358979323846



//*****************************************************************//
// define boys function
///////////////////////////////

double boys(double t){ 

  double pre = .5 * pow((pi/t),0.5);
  return pre * erf(pow(t,0.5));
}



//*****************************************************************//
// overlap matrix S
///////////////////////////////

double sElt(double p_oN, double q_oN, prim orbP[], prim orbQ[], int pLen, int qLen){

  // p_oN, q_oN, normalization for total gtos
  // orbP[], orbQ[], arrays of prim orbs
  // pLen, qLen, nums of prim orbs for each total gto

  // initialize indices
  // and value to be returned
  int p,q;
  double total = 0.0;

  // initialize quantities for calculation
  vector vecDiff;
  double distSq, pre1, pre2, exponent;

  // loop over all primitive orbitals
  for( p = 0 ; p < pLen ; p++){
    for( q = 0 ; q < qLen ; q++){

    // get sq distances btw prim orbs
    vecDiff = retSubtractedVectors(orbP[p].ctr, orbQ[q].ctr);
    distSq = retDotProduct(vecDiff, vecDiff);

    // compute factors
    pre1 = p_oN * q_oN * orbP[p].co * orbQ[q].co * orbP[p].pN * orbQ[q].pN;
    pre2 = pow((pi/(orbP[p].ex + orbQ[q].ex)), 1.5);
    exponent = exp(-1.0 * ((orbP[p].ex*orbQ[q].ex)/(orbP[p].ex+orbQ[q].ex))*distSq);

    // increment total
    total += pre1 * pre2 * exponent;

    }
  }
  return total;
}

//*****************************************************************//



//*****************************************************************//
// kinetic energy matrix T
///////////////////////////////

double tElt(double p_oN, double q_oN, prim orbP[], prim orbQ[], int pLen, int qLen){

  // p_oN, q_oN, normalization for total gtos
  // orbP[], orbQ[], arrays of prim orbs
  // pLen, qLen, nums of prim orbs for each total gto


  // initialize indices
  // and value to be returned
  int p,q;
  double total = 0.0;
   
  // initialize quantities for calculation
  vector vecDiff;
  double distSq, pre1, pre2, pre3, pre4, exponent;
   
  // loop over all primitive orbitals
  for( p = 0 ; p < pLen ; p++){
    for( q = 0 ; q < qLen ; q++){

    // get sq distances btw prim orbs
    vecDiff = retSubtractedVectors(orbP[p].ctr, orbQ[q].ctr);
    distSq = retDotProduct(vecDiff, vecDiff);

    // compute factors
    pre1 = p_oN * q_oN * orbP[p].co * orbQ[q].co * orbP[p].pN * orbQ[q].pN;
    pre2 = -2.0 * pow(((orbP[p].ex*orbQ[q].ex)/(orbP[p].ex + orbQ[q].ex)), 2.0) * distSq;
    pre3 = 3.0 * ((orbP[p].ex*orbQ[q].ex)/(orbP[p].ex + orbQ[q].ex));
    pre4 = pow((pi/(orbP[p].ex + orbQ[q].ex)), 1.5);
    exponent = exp(-1.0 * ((orbP[p].ex*orbQ[q].ex)/(orbP[p].ex+orbQ[q].ex))*distSq);

    // increment total
    total += pre1 * (pre2 + pre3) * pre4 * exponent;
    }
  }
  return total;
}

//*****************************************************************//



//*****************************************************************//
// nuclear potential energy matrix Vn
///////////////////////////////

double vnElt(double p_oN, double q_oN, prim orbP[], prim orbQ[], int pLen, int qLen, 
  nuc nucleus){

  // p_oN, q_oN, normalization for total gtos
  // orbP[], orbQ[], arrays of prim orbs
  // pLen, qLen, nums of prim orbs for each total gto
  // z, charge of nucleus
  // nucCoords, position of nucleus
  
  
  // initialize indices
  // and value to be returned
  int p,q;
  double total = 0.0;

  // initialize quantities for calculation
  vector vecDiff, rp, rpDist;
  double distSq, rpFac, rpSq, pre1, pre2, exponent, boysFac;
  
  // loop over all primitive orbitals
  for( p = 0 ; p < pLen ; p++){
    for( q = 0 ; q < qLen ; q++){
  
    // get sq distances btw prim orbs
    vecDiff = retSubtractedVectors(orbP[p].ctr, orbQ[q].ctr);
    distSq = retDotProduct(vecDiff, vecDiff);

    // get square difference between weighted avg of prim orbs centers
    // and nuc coords
    rpFac = 1.0/(orbP[p].ex + orbQ[q].ex);
    rp = scaleVector(retAddedVectors(
      scaleVector(orbP[p].ctr, orbP[p].ex), 
      scaleVector(orbQ[q].ctr, orbQ[q].ex)
      ), rpFac);
    rpDist = retSubtractedVectors(rp, nucleus.coords);
    rpSq = retDotProduct(rpDist, rpDist);

    // compute factors
    pre1 = p_oN * q_oN * orbP[p].co * orbQ[q].co * orbP[p].pN * orbQ[q].pN;
    pre2 = -2.0 * nucleus.z * (pi/(orbP[p].ex+orbQ[q].ex));
    exponent = exp(-1.0 * ((orbP[p].ex*orbQ[q].ex)/(orbP[p].ex+orbQ[q].ex))*distSq);

    // do orbs have the same coords?
    if (rpSq < 0.00001) {
      boysFac = 1.0;
    }
    else {
      boysFac = boys((orbP[p].ex + orbQ[q].ex)*rpSq);
    }

    // increment total val
    total += pre1 * pre2 * exponent * boysFac;
    }
  }
  return total;
}

//*****************************************************************//



//*****************************************************************//
// yeet eri
///////////////////////////////

double eriElt(double p_oN, double q_oN, double r_oN, double s_oN,
  prim orbP[], prim orbQ[], prim orbR[], prim orbS[],
  int pLen, int qLen, int rLen, int sLen) {
  
  // p_oN, q_oN, r_oN, s_oN, normalization for total gtos
  // orbP[], orbQ[], orbR[], orbS[], arrays of prim orbs
  // pLen, qLen, rLen, sLen, nums of prim orbs for each total gto

  vector distPQ, distRS, pqAvg, rsAvg, wDist;
  double pqSq, rsSq, pqAvgFac, rsAvgFac, wDistSq,
    rhoNumer, rhoDenom, rho,
    pre1, pre2, pre3, pre4,
    exponent1, exponent2, boysFac;

  // initialize indices
  // and value to be returned
  int p,q,r,s;
  double total = 0.0;

  // loop over all primitive orbitals
  for( p = 0 ; p < pLen ; p++){
    for( q = 0 ; q < qLen ; q++){
      for( r = 0 ; r < rLen ; r++){
        for( s = 0 ; s < sLen ; s++){

        // get square distance between prim orbs
        distPQ = retSubtractedVectors(orbP[p].ctr, orbQ[q].ctr);
        distRS = retSubtractedVectors(orbR[r].ctr, orbS[s].ctr);
        pqSq = retDotProduct(distPQ, distPQ);
        rsSq = retDotProduct(distRS, distRS);

        // get weighted average of prim orbs ctrs
        // and get square distance
        pqAvgFac = 1.0/(orbP[p].ex + orbQ[q].ex);
        rsAvgFac = 1.0/(orbR[r].ex + orbS[s].ex);
        pqAvg = scaleVector(retAddedVectors(
          scaleVector(orbP[p].ctr, orbP[p].ex),
          scaleVector(orbQ[q].ctr, orbQ[q].ex)
          ), pqAvgFac);
        rsAvg = scaleVector(retAddedVectors(
          scaleVector(orbR[r].ctr, orbR[r].ex),
          scaleVector(orbS[s].ctr, orbS[s].ex)
          ), rsAvgFac);
        wDist = retSubtractedVectors(pqAvg, rsAvg);
        wDistSq = retDotProduct(wDist, wDist);

        // compute rho
        rhoNumer = (orbP[p].ex + orbQ[q].ex) * (orbR[r].ex + orbS[s].ex);
        rhoDenom = orbP[p].ex + orbQ[q].ex + orbR[r].ex + orbS[s].ex;
        rho = rhoNumer / rhoDenom;

        // compute factors
        pre1 = p_oN * q_oN * r_oN * s_oN * orbP[p].pN * orbQ[q].pN * orbR[r].pN * orbS[s].pN;
        pre2 = orbP[p].co * orbQ[q].co * orbR[r].co * orbS[s].co;
        pre3 = 2.0*pow(pi,2.5);
        pre4 = 1.0/(rhoNumer * pow(rhoDenom,0.5));
        exponent1 = exp(-1.0*((orbP[p].ex*orbQ[q].ex)/(orbP[p].ex+orbQ[q].ex))*pqSq);
        exponent2 = exp(-1.0*((orbR[r].ex*orbS[s].ex)/(orbR[r].ex+orbS[s].ex))*rsSq);

        // do orbs have the same coords?
        if (wDistSq < 0.00001) {
          boysFac = 1.0;
        }
        else {
          boysFac = boysFac = boys((rhoNumer / rhoDenom)*wDistSq);
        }

        // increment total val
        total += pre1 * pre2 * pre3 * pre4 * exponent1 * exponent2 * boysFac;
        }
      }
    }
  }
  return total;
}
