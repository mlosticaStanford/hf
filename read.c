#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include"vector.h"
#include"basis.h"
#include"read.h"

double assignZ(char *s) {
  if (!strcmp("H",s)) {
    return 1.0;
  }
  else if (!strcmp("He",s)) {
    return 2.0;
  }
}

// get number of atoms in a structure
int getNat(char *fName) {

  // fName, name of file to be read
  // !!! must be "input.txt" !!!
  
  FILE *fp;
  int nat;
  char buff[1024];
  // fp, file 

  // open file
  fp = fopen(fName, "r");

  // first line is "$molecule"
  fscanf(fp, "%1023s", buff);

  // second line is nat
  fscanf(fp, "%i", &nat);

  // close file
  fclose(fp);

  // return val
  return nat;
}

// read molecule
void readMolecule(char *fName, nuc *nuclei[], int nat) {
  FILE *fp;
  int buffLen = 1024;
  char buff[buffLen], atomSymbol[buffLen];
  char val[buffLen], *b;
  int i = 0; 

  fp = fopen(fName, "r");
  rewind(fp);
  
  // read first line
  fscanf(fp, "%1023s", buff);

  // next line is number of atoms
  fscanf(fp, "%1023s", buff);


  while (i < nat) {
    
    // new nuc instance to be added to array
    nuc *a;
    a = (nuc *)malloc(sizeof(nuc) + sizeof(double)*5 + sizeof(char)*1024);

    // read atom symbol and assign it
    fscanf(fp, "%1024s", atomSymbol);
    strcpy(a->atmSym, atomSymbol);

    // read and then assign x, y, z coords
    fscanf(fp, "%1024s", buff);
    strcpy(val, buff);
    a->coords.x = strtod(val, &b);
    //
    fscanf(fp, "%1024s", buff);
    strcpy(val, buff);
    a->coords.y = strtod(val, &b);
    //
    fscanf(fp, "%1024s", buff);
    strcpy(val, buff);
    a->coords.z = strtod(val, &b);
    

    // assign charge based on atom symbol
    a->z = assignZ(atomSymbol);

    // put the nucleus into the array
    nuclei[i] = a;

    // read next nucleus
    i++;
  }

  // close file
  fclose(fp);
}

// count total number of GTOs in a basis file
int getNorbs(char *fName) {
  FILE *fp;
  int i = 0;
  int buffLen = 1024;
  char buff[buffLen];

  // open file
  fp = fopen(fName, "r");

  // "1.00" appears as many times as
  // there are total GTOs in a file
  while (fscanf(fp, " %1024s", buff) == 1) {
    if(!strcmp("1.00",buff)) {
      i++;
    }
  }

  // close file
  fclose(fp);

  return i;
}


void makeOrbs(int nat, int *orbCount, char *bName, nuc *nuclei[], orb *allOrbs[])  {
  FILE *fp;
  int i,j,k,l,nOrbs,nPrims;
  int buffLen = 1024;
  char fName[buffLen];
  char buff[buffLen];
  char val[buffLen], *b;

  // iterate over all nuclei
  for (i = 0; i < nat; i++) {

    // get filename for basis
    strcpy(fName,nuclei[i]->atmSym);
    strcat(fName,"-");
    strcat(fName,bName);
    strcat(fName,".txt");

    // get number of GTOs in the basis file
    nOrbs = getNorbs(fName);
    *orbCount += nOrbs;

    // open file
    fp = fopen(fName, "r");
  
    // skip first two lines
    fscanf(fp, "%1023s", buff);
    fscanf(fp, "%1023s", buff);
    fscanf(fp, "%1023s", buff);
    
    for (j = 0; j < nOrbs; j++) {

      // compound indices
      // this is used to index allOrbs
      k = i + j;

      // this is just S,P,SP,D...
      fscanf(fp, "%1023s", buff);

      // get number of prims
      fscanf(fp, "%1024s", buff);
      strcpy(val, buff);
      nPrims = atoi(val);

      // initialize GTO to be constructed
      orb *gto;
      gto = (orb *)malloc(sizeof(orb) + sizeof(int) + 
        sizeof(double) + sizeof(prim)*nPrims);

      // assign gto its number of prims
      gto -> orbLen = nPrims;

      // assign gto to allOrbs
      allOrbs[k] = gto;

      // this is just 1.00
      fscanf(fp, "%1023s", buff);

      // read each prim
      for (l = 0; l < nPrims; l++) {
        
        // initialize prim to be constructed
        prim *p;
        p = (prim *)malloc(sizeof(prim) + sizeof(double)*7);

        // assign ctr to prim
        p->ctr.x = nuclei[i]->coords.x;
        p->ctr.y = nuclei[i]->coords.y;
        p->ctr.z = nuclei[i]->coords.z;

        // assign exponent
        fscanf(fp, "%1024s", buff);
        strcpy(val, buff);
        p->ex = strtod(val, &b);
        
        // assign coefficient
        fscanf(fp, "%1024s", buff);
        strcpy(val, buff);
        p->co = strtod(val, &b);

        // assign prim to gto
        gto->orbPrims[l] = *p;
      }
    }
  }
}

// read params
void readParams(char *fName, calcParams *params) {
  FILE *fp;
  int buffLen = 1024;
  char buff[buffLen];
  char val[buffLen], *b;
  int i = 0; 

  fp = fopen(fName, "r");
  rewind(fp);
  
  fscanf(fp, "%1023s", buff);
  strcpy(val, buff);

  // read until $calcParams
  while (strcmp(val,"$calcParams")) {
    fscanf(fp, "%1023s", buff);
    strcpy(val, buff);
  }

  // read basis
  fscanf(fp, "%1023s", buff); fscanf(fp, "%1023s", buff); fscanf(fp, "%1023s", buff);
  strcpy(params->basis, buff);

  // read uhf boolean
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%i", &i);
  params->uhf = i;
  printf("params->uhf = %i\n", params->uhf);

  // read tolerance
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%1023s", buff);
  strcpy(val,buff);
  params->tol = strtod(val, &b);
  printf("params->tol= %.9f\n", params->tol);

  // alpha elecs
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%i", &i);
  params->aElecs = i;
  printf("params->aElecs= %i\n", params->aElecs);

  // beta elecs
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%i", &i);
  params->bElecs = i;
  printf("params->bElecs= %i\n", params->bElecs);

  // max iters
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%i", &i);
  params->maxIters= i;

  // diisNum 
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%1023s", buff);
  fscanf(fp, "%i", &i);
  params->diisNum = i;

  // close file
  fclose(fp);
}
