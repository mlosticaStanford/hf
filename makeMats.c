//*****************************************************************//
// functions for forming matrices
//*****************************************************************//

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "matOps.h"
#include "basis.h"
#include "vector.h"
#include "matrixElts.h"
#include "makeMats.h"



//*****************************************************************//
// compute compound indices for ERI
///////////////////////////////

int compoundIndices(int blen, int i, int j, int k, int l) {
  // blen, length of basis
  // i,j,k,l indices to be compounded

  // this polynomial in blen works
  // and still uses memore that is of O(blen^4)
  // expoiting the symmetry of the tensor would be better
  return i + j*blen + k*blen*blen +l*blen*blen*blen;
}

//*****************************************************************//



//*****************************************************************//
// overlap matrix S
///////////////////////////////

void makeS(int blen, orb *allOrbs[], double S[]) {
  // blen, length of basis
  // allOrbs[], orbitals used in calculation
  // S[], initialized s matrix

  // iterate over each basis function
  int i,j,k;
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {

      // 1D index
      k = i + (j*blen);

      // compute matrix elt
      S[k] = sElt(allOrbs[i] -> oN, allOrbs[j] -> oN, 
        allOrbs[i] -> orbPrims, allOrbs[j] -> orbPrims, 
        allOrbs[i] -> orbLen, allOrbs[j] -> orbLen);
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// kinetic energy matrix T
///////////////////////////////

void makeT(int blen, orb *allOrbs[], double T[]) {
  // blen, length of basis
  // allOrbs[], orbitals used in calculation
  // T[], initialized t matrix

  // iterate over each basis function
  int i,j,k;
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {

      // 1D index
      k = i + (j*blen);
      T[k] = tElt(allOrbs[i] -> oN, allOrbs[j] -> oN, 
        allOrbs[i] -> orbPrims, allOrbs[j] -> orbPrims, 
        allOrbs[i] -> orbLen, allOrbs[j] -> orbLen);
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// nuclear potential energy matrix Vn
///////////////////////////////

void makeVn(int blen, int nats, orb *allOrbs[], double Vn[], nuc *nuclei[]) {
  // blen, length of basis
  // nats, number of nuclei
  // allOrbs[], orbitals used in calculation
  // Vn[], initialized vn matrix
  // nuclei, array of nucleus instances 

  // iterate over each basis function
  int i,j,k,a;
  double total; 
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {
      
      // 1D index
      k = i + (j*blen);

      // iterate over each nucleus
      total = 0.0;
      for (a = 0; a < nats; a++) {
        total += vnElt(allOrbs[i] -> oN, allOrbs[j] -> oN, 
          allOrbs[i] -> orbPrims, allOrbs[j] -> orbPrims, 
          allOrbs[i] -> orbLen, allOrbs[j] -> orbLen,
          *nuclei[a]);
      }

      // assign total Vn matrix elt
      Vn[k] = total;
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// core hamiltonian Hcore
///////////////////////////////

void makeHcore(int blen, int nats, orb *allOrbs[], 
  double T[], double Vn[], double Hcore[], nuc *nuclei[]) {
  // blen, length of basis
  // nats, number of nuclei
  // allOrbs[], orbitals used in calculation
  // Vn[], initialized vn matrix
  // nuclei, array of nucleus instances 

  // iterate over each basis function
  int i,j,k,a;
  double total; 
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {

      // 1D index
      k = i + (j*blen);

      // compute kinetic energy elt
      T[k] = tElt(allOrbs[i] -> oN, allOrbs[j] -> oN, 
        allOrbs[i] -> orbPrims, allOrbs[j] -> orbPrims, 
        allOrbs[i] -> orbLen, allOrbs[j] -> orbLen);

      // compute e-n attraction elt
      total = 0.0;
      for (a = 0; a < nats; a++) {
        total += vnElt(allOrbs[i] -> oN, allOrbs[j] -> oN, 
          allOrbs[i] -> orbPrims, allOrbs[j] -> orbPrims, 
          allOrbs[i] -> orbLen, allOrbs[j] -> orbLen,
          *nuclei[a]);
      }
      Vn[k] = total;

      // assign value to Hcore
      Hcore[k] = T[k] + Vn[k];
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// yeet eri
///////////////////////////////

void makeERI(int blen, orb *allOrbs[], double eri[]) {
  // blen, length of basis
  // allOrbs[], orbitals used in calculation
  // eri[], initialized eri tensor

  // iterate over each basis function
  int i,j,k,l,n;
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {
      for (k = 0; k < blen; k++) {
        for (l = 0; l < blen; l++) {

          // compound indices
          n = compoundIndices(blen, i,j,k,l);

          // hella params for eriElt lol
          eri[n] = eriElt(allOrbs[i] -> oN, allOrbs[j] -> oN, allOrbs[k] -> oN, allOrbs[l] -> oN
            , allOrbs[i] -> orbPrims, allOrbs[j] -> orbPrims, allOrbs[k] -> orbPrims, allOrbs[l] -> orbPrims
            , allOrbs[i] -> orbLen, allOrbs[j] -> orbLen, allOrbs[k] -> orbLen, allOrbs[l] -> orbLen);
        }
      }
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// make x matrix for canonical orthogonalization of s
///////////////////////////////

void makeX(int blen, double S[], double U[], double sp[], 
  double X[], double Xt[]) {
  // blen, length of basis
  // allOrbs[], orbitals used in calculation
  // S, computed overlap matrix
  // U, initialized matrix to store diagonalized U
  // sp, array for eVals of S
  // X,Xt, initialized matrices to be output in col form

  // iterate over each basis function
  int i,j,k;
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {

      // 1D index
      k = i + (j*blen);

      // note that U is stored in column major format
      // !!! and this matters lol !!!
      X[k] = U[k] / pow(sp[j],0.5);
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// make g matrix
///////////////////////////////

void makeG(int blen, double P[], double eri[], double G[]) {
  // blen, length of basis
  // P, density matrix
  // eri, eri tensor
  // G, initialized matrix
  // S, computed overlap matrix
  // U, initialized matrix to store diagonalized U

  // initialize indices
  int i,j,k,l;
  int en,em,pkl,gij;

  // vars for computation
  double gElt, eriFac;

  // iterate over each basis function
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {

      // 1D index
      gij = i + j*blen;

      // (Re)initialize matrix elt
      gElt = 0.0;

      // iterate over eri elts
      for (k = 0; k < blen; k++) {
        for (l = 0; l < blen; l++) {

          // compound indices
          en = i + j*blen + k*blen*blen + l*blen*blen*blen;
          em = i + k*blen + l*blen*blen + j*blen*blen*blen;
          pkl = k + l*blen;

          // increment matrix elt
          eriFac = eri[en] - (.5*eri[em]);
          gElt += P[pkl]*eriFac;
        }
      }

      // assign total val to G matrix
      G[gij] = gElt;
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// make uhf g matrix
///////////////////////////////

void makeUrG(int blen, double Pa[], double Pb[], double Pt[],
  double eri[], double Ga[], double Gb[]) {
  // blen, length of basis
  // Pa,Pb,Pt density matrices
  // eri, eri tensor
  // Ga,Gb initialized matrices
  // S, computed overlap matrix
  // U, initialized matrix to store diagonalized U

  // initialize indices
  int i,j,k,l;
  int en,em,pkl,gij;

  // vars for computation
  double gAelt, gBelt;

  // form total density matrix Pt
  addMats(blen, blen, false, Pa, Pb, 
    1.0, 1.0, Pt);

  // iterate over each basis function
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {

      // 1D index
      gij = i + j*blen;

      // (Re)initialize matrix elt
      gAelt = 0.0;
      gBelt = 0.0;

      // iterate over eri elts
      for (k = 0; k < blen; k++) {
        for (l = 0; l < blen; l++) {

          // compound indices
          en = i + j*blen + k*blen*blen + l*blen*blen*blen;
          em = i + k*blen + l*blen*blen + j*blen*blen*blen;
          pkl = k + l*blen;

          // increment matrix elts
          gAelt += (Pt[pkl]*eri[em]) - (Pa[pkl]*eri[en]);
          gBelt += (Pt[pkl]*eri[em]) - (Pb[pkl]*eri[en]);
        }
      }

      // assign total val to G matrices
      Ga[gij] = gAelt;
      Gb[gij] = gBelt;
    }
  }
}

//*****************************************************************//



//*****************************************************************//
// update density matrix
///////////////////////////////

void upP(int ur, int blen, int nElec, double P[], double C[]) { 
  // ur, true == uhf
  // blen, length of basis
  // nElec, number of electrons
  // P, initialized density matrix
  // C, transformed AO coeffs

  // initialize indices
  int i,j,k,kMax;
  int pij, cik, cjk;
  
  // vars for computation
  double pElt, pFac;

  // assign vals to vars according to ur
  if (ur == 0) {
    kMax = nElec / 2;
    pFac = 2.0;
  }
  else {
    kMax = nElec;
    pFac = 1.0;
  }


  // iterate over each basis function
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {

      // 1D index
      pij = i + j*blen;

      // (Re)initialize matrix elt
      pElt = 0.0;

      // iterate over half the number of electrons
      for (k = 0; k < kMax; k++) {

        // 1D indices
        cik = i + k*blen;
        cjk = j + k*blen;

        // increment matrix elt
        pElt += pFac * C[cik] * C[cjk];
      }

      // assign total val to P matrix
      P[pij] = pElt;
    }
  }
}
