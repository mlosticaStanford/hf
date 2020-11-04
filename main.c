#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"vector.h"
#include"basis.h"
#include"read.h"
#include "matrixElts.h"
#include "makeMats.h"
#include "matOps.h"
#include "scf.h"

int main() {

  // get number of atoms in calc
  int nat = getNat("input.txt");
  int *orbCount;
  calcParams *params;

  orbCount = (int *)malloc(sizeof(int));
  params = (calcParams *)malloc(sizeof(char)*1024 + 
    sizeof(int)*5 + sizeof(double));

  // read parameters
  readParams("input.txt", params);

  // make list of nuclei
  nuc *nuclei[nat];
  readMolecule("input.txt", nuclei, nat);

  // make list of orbs
  orb *orbs[nat*8];

  // make orbs for calculation
  makeOrbs(nat, orbCount, params->basis, nuclei, orbs);

  // normalize orbitals
  normAllOrbs(*orbCount, orbs);

  // just do it
  scf(params->uhf, nat, *orbCount, params->maxIters, params->tol, 
    nuclei, orbs, params->aElecs, params->bElecs, params->diisNum);
}
