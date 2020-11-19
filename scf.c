// functions for running scf
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include "matOps.h"
#include "basis.h"
#include "vector.h"
#include "matrixElts.h"
#include "makeMats.h"
#include "scf.h"
#include "diis.h"

// calculate nuclear repulsion energy for a set of nuclei
double calcTotalNucE(int nat, nuc *nuclei[]) {
  int i,j;
  double distSq, dist, e = 0.0;
  vector diff;
  
  for (i = 0; i < nat; i++) {
    for (j = 0; j < i; j++) {
      diff = retSubtractedVectors(nuclei[i] -> coords, nuclei[j] -> coords);
      distSq = retDotProduct(diff, diff);
      dist = pow(distSq, 0.5);
      e += ((nuclei[i] -> z) * (nuclei[j] -> z)) / dist;
    }
  }
  return e;
}


// calculate ground state energy
// used to evaluate convergence
double getE(int blen, double P[], double Hcore[], double F[]) {

  // initialize val to be returned
  // and indices for matrix elts
  double e = 0.0;
  int i, j, pji, hij, fij;

  // calculate energy
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {

      // 1D indices
      pji = i*blen + j;
      hij = i + j*blen;
      fij = i + j*blen;

      // increment e
      e += 0.5 * P[pji] * (Hcore[hij] + F[fij]);
    }
  }

  // return ground state electronic energy
  return e;
}


// uhf energy
double getUrE(int blen, double Pa[], double Pb[], double Pt[],
  double Fa[], double Fb[], double Hcore[]) {

  // form total density matrix Pt
  addMats(blen, blen, false, Pa, Pb, 
    1.0, 1.0, Pt);
  
  // initialize val to be returned
  // and indices for matrix elts
  double e = 0.0;
  int i, j, pji, hij, fij;

  // increment energy
  for (i = 0; i < blen; i++) {
    for (j = 0; j < blen; j++) {
      
      // 1D indices
      pji = i*blen + j;
      hij = i + j*blen;
      fij = i + j*blen;

      // increment e
      e += 0.5 * Pt[pji]*Hcore[hij];
      e += 0.5 * Pa[pji]*Fa[fij];
      e += 0.5 * Pb[pji]*Fb[fij];
    }
  }

  // return ground state electronic energy
  return e;
} 



//////////////////////////////
// run scf calculation
//////////////////////////////

void scf(int ur, int nat, int blen, int maxIters, double eTol, 
  nuc *nuclei[], orb *allOrbs[], int nA, int nB, int diisNum) {


  //////////////////////////////
  // open file for output
  //////////////////////////////
  FILE *fp;
  fp = fopen("output.txt", "w");
  fputs("Running SCF calculation\n",fp);

  //////////////////////////////
  // initialize matrices and allocate memory
  //////////////////////////////
  
  // matrices commmon to rhf and uhf
  double *S, *U, *sp, *X, *Xt, *Xtf, *T, *Vn, *Hcore, *eri, *workA, *workB; 

  S = (double *)malloc(sizeof(double)*blen*blen);
  U = (double *)malloc(sizeof(double)*blen*blen);
  sp = (double *)malloc(sizeof(double)*blen);
  X = (double *)malloc(sizeof(double)*blen*blen);
  Xt = (double *)malloc(sizeof(double)*blen*blen);
  Xtf = (double *)malloc(sizeof(double)*blen*blen);
  T = (double *)malloc(sizeof(double)*blen*blen);
  Vn = (double *)malloc(sizeof(double)*blen*blen);
  Hcore = (double *)malloc(sizeof(double)*blen*blen);
  eri = (double *)malloc(sizeof(double)*(pow(blen,4)+pow(blen,3)+blen*blen+blen));
  workA = (double *)malloc(sizeof(double)*blen*blen);
  workB = (double *)malloc(sizeof(double)*blen*blen);

  // rhf mats
  double *P, *G, *C, *Cp, *eigs, *F, *Fp;

  // uhf mats
  double *Xtfa, *Xtfb, *Pa, *Pb, *Pt, *Ga, *Gb, *Cpa, *Cpb; 
  double *Ca, *Cb, *eigsA, *eigsB, *Fa, *Fb, *Fpa, *Fpb;

  // diis stuff
  double errNorm = 1.0;
  double *err, *coeffs;
  double *bCopy;

  // rhf diis
  double *errVecs[diisNum], *fockArr[diisNum]; 
  double *B; 

  // uhf diis
  double *errVecsA[diisNum], *fockArrA[diisNum], 
    *errVecsB[diisNum], *fockArrB[diisNum]; 
  double *Ba, *Bb;


  // rhf
  if (ur == 0) {
    P = (double *)malloc(sizeof(double)*blen*blen);
    G = (double *)malloc(sizeof(double)*blen*blen);
    C = (double *)malloc(sizeof(double)*blen*blen);
    Cp = (double *)malloc(sizeof(double)*blen*blen);
    eigs = (double *)malloc(sizeof(double)*blen);
    F = (double *)malloc(sizeof(double)*blen*blen);
    Fp = (double *)malloc(sizeof(double)*blen*blen);

    // if diis
    if (diisNum > 0) {
      
      // initialize arrays
      initVecArray(blen, diisNum, errVecs);
      initVecArray(blen, diisNum, fockArr);

      // initialize B and error matrices, and coeffs vec
      B = (double *)malloc(sizeof(double) * pow(diisNum+1,2)); 
      bCopy = (double *)malloc(sizeof(double) * pow(diisNum+1,2)); 
      err = (double *)malloc(sizeof(double) * pow(blen,2)); 
      coeffs = (double *)malloc(sizeof(double) * (diisNum + 1)); 
      initB(diisNum, B);
      initB(diisNum, bCopy);
    }
  }

  // uhf
  else {
    Xtfa = (double *)malloc(sizeof(double)*blen*blen);
    Xtfb = (double *)malloc(sizeof(double)*blen*blen);
    Pa = (double *)malloc(sizeof(double)*blen*blen);
    Pb = (double *)malloc(sizeof(double)*blen*blen);
    Pt = (double *)malloc(sizeof(double)*blen*blen);
    Ga = (double *)malloc(sizeof(double)*blen*blen);
    Gb = (double *)malloc(sizeof(double)*blen*blen);
    Ca = (double *)malloc(sizeof(double)*blen*blen);
    Cb = (double *)malloc(sizeof(double)*blen*blen);
    Cpa = (double *)malloc(sizeof(double)*blen*blen);
    Cpb = (double *)malloc(sizeof(double)*blen*blen);
    eigsA = (double *)malloc(sizeof(double)*blen);
    eigsB = (double *)malloc(sizeof(double)*blen);
    Fa = (double *)malloc(sizeof(double)*blen*blen);
    Fb = (double *)malloc(sizeof(double)*blen*blen);
    Fpa = (double *)malloc(sizeof(double)*blen*blen);
    Fpb = (double *)malloc(sizeof(double)*blen*blen);

    Pa[0] = .50625 * (nA+nB);
    Pb[(blen*blen) - 1] = .49375 * (nA+nB);

    if (diisNum > 0) {

      // initialize errVec arrays
      initVecArray(blen, diisNum, errVecsA);
      initVecArray(blen, diisNum, errVecsB);

      // initialize fock arrays
      initVecArray(blen, diisNum, fockArrA);
      initVecArray(blen, diisNum, fockArrB);


      // initialize B and error matrices, and coeffs vec
      Ba = (double *)malloc(sizeof(double) * pow(diisNum+1,2)); 
      Bb = (double *)malloc(sizeof(double) * pow(diisNum+1,2)); 
      bCopy = (double *)malloc(sizeof(double) * pow(diisNum+1,2)); 
      err = (double *)malloc(sizeof(double) * pow(blen,2)); 
      coeffs = (double *)malloc(sizeof(double) * (diisNum + 1)); 
      initB(diisNum, Ba);
      initB(diisNum, Bb);
      initB(diisNum, bCopy);
    }

    //int i,j,k,l;

    //for (i = 0; i < blen; i++) {
    //  for (j = 0; j <= i; j++) {
    //    
    //    k = i + j*blen;
    //    l = i*blen + j;

    //    if (i == j && j == 0) {
    //      Pa[k] = 1.0;
    //      Pb[k] = 1.0;
    //    }

    //    else {
    //      Pa[k] = Pa[l] = 0.0;
    //      Pb[k] = Pb[l] = 0.0;
    //    }
    //  }
    //} 
  }
  

  //////////////////////////////
  // calculate nuclear repulsion energy
  //////////////////////////////
  double nucE = calcTotalNucE(nat, nuclei);


  //////////////////////////////
  // make matrices for Hcore, also make eri
  //////////////////////////////

  // compute T, Vn, Hcore
  makeHcore(blen, nat, allOrbs, T, Vn, Hcore, nuclei); 

  // compute eri
  makeERI(blen, allOrbs, eri);


  //////////////////////////////
  // overlap matrix and orthogonalization
  //////////////////////////////
  
  // make S and diagonalize it
  makeS(blen, allOrbs, S);
  copyMat(blen, blen, S, U);
  diagMat(blen, blen, U, sp);

  // make X, Xt
  makeX(blen, S, U, sp, X, Xt); 
  transpose(blen, blen, false, false, X, Xt);


  //////////////////////////////
  // setup exit conditions for iterations
  //////////////////////////////
  bool notConverged = true;
  int iters = 0;
  double oldE = 1.0, newE = 1.0, dE = 1.0;  // dummy vals


  //////////////////////////////
  // do scf iterations
  //////////////////////////////
  
  while (notConverged) {
    printf("iters = %i\n", iters);

    // rhf
    if (ur == 0) {
    
      // compute G
      makeG(blen, P, eri, G);

      // compute F
      addMats(blen, blen, false, 
        Hcore, G, 1.0, 1.0, F);

      // if diis
      if (diisNum > 0) {

          // run diis with given params
          runDIIS(blen, diisNum, iters,
            fockArr, errVecs,
            B, coeffs, 
            F, err,
            P, S, X,
            workA, workB, bCopy,
            errNorm);
      } 

      // compute Fp
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, Xt, blen, F, blen,
                                        0.0, Xtf, blen);
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, Xtf, blen, X, blen,
                                        0.0, Fp, blen);

      // diagonalize Fp
      copyMat(blen, blen, Fp, Cp);
      diagMat(blen, blen, Cp, eigs);

      // compute C
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, X, blen, Cp, blen,
                                        0.0, C, blen);

      // compute ground state electronic energy
      newE = getE(blen, P, Hcore, F);

      // update density matrix
      upP(ur, blen, nA, P, C);  

    }

    // uhf
    else {

      printf("doing uhf\n");

    
      // compute G's
      makeUrG(blen, Pa, Pb, Pt, eri, Ga, Gb);

      // compute F's
      addMats(blen, blen, false, 
        Hcore, Ga, 1.0, 1.0, Fa);
      addMats(blen, blen, false, 
        Hcore, Gb, 1.0, 1.0, Fb);

      // if diis
      if (diisNum > 0) {

        // diis for Fock A
        runDIIS(blen, diisNum, iters,
          fockArrA, errVecsA,
          Ba, coeffs, 
          Fa, err,
          Pa, S, X,
          workA, workB, bCopy,
          errNorm);

        // diis for Fock B
        runDIIS(blen, diisNum, iters,
          fockArrB, errVecsB,
          Bb, coeffs, 
          Fb, err,
          Pb, S, X,
          workA, workB, bCopy, 
          errNorm);
     } 

      // compute Fp's
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, Xt, blen, Fa, blen,
                                        0.0, Xtfa, blen);
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, Xtfa, blen, X, blen,
                                        0.0, Fpa, blen);
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, Xt, blen, Fb, blen,
                                        0.0, Xtfb, blen);
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, Xtfb, blen, X, blen,
                                        0.0, Fpb, blen);

      // diagonalize Fp's
      copyMat(blen, blen, Fpa, Cpa);
      diagMat(blen, blen, Cpa, eigsA);
      copyMat(blen, blen, Fpb, Cpb);
      diagMat(blen, blen, Cpb, eigsB);

      // compute C's
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, X, blen, Cpa, blen,
                                        0.0, Ca, blen);
	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                        blen, blen, blen,
                                        1.0, X, blen, Cpb, blen,
                                        0.0, Cb, blen);


      // compute ground state electronic energy
      newE = getUrE(blen, Pa, Pb, Pt,
        Fa, Fb, Hcore); 

      // update P's
      upP(ur, blen, nA, Pa, Ca);  
      upP(ur, blen, nB, Pb, Cb);  

    }

    // evaluate convergence criteria
    dE = fabs(newE - oldE);
    notConverged = (dE > eTol) && (iters < maxIters);
    iters++;
    oldE = newE;
    // if diis, use errVec for convergence
    //if (diisNum > 0) {
    //  notConverged = (errNorm > eTol) && (iters < maxIters);
    //}
  }

  // calculate total energy
  double totalE = newE + nucE;

  // did we converge?
  if (iters >= maxIters) {
    fprintf(fp,"Maximum iterations reached! \n"); 
    fprintf(fp,"At the last iteration, we have that \n");
    fprintf(fp,"The ground state electronic energy = %f \n", newE); 
    fprintf(fp,"Final ground state energy = %f\n", totalE);
  }
  else {
    fprintf(fp,"Convergence criteria satisfied! \n");
    fprintf(fp,"The ground state electronic energy = %f \n", newE); 
    fprintf(fp,"Total ground state energy = %f\n", totalE);
  }

  // free memory
  free(S); free(U); free(sp); free(X); free(Xt); free(Xtf);
  free(T); free(Vn); free(Hcore); free(eri); free(workA); free(workB);
  
  if (ur == 0) {
    free(P); free(G); free(C); free(Cp);
    free(eigs); free(F); free(Fp);
  }

  else {
    free(Xtfa); free(Xtfb); free(Pa); free(Pb); free(Pt); 
    free(Ga); free(Gb); free(Ca); free(Cb); 
    free(Cpa); free(Cpb); free(eigsA); free(eigsB); 
    free(Fa); free(Fb); free(Fpa); free(Fpb);
  }
  fclose(fp);
}
