//*****************************************************************//
// functions for computing matrix elements
//*****************************************************************//

#include "vector.h"
#include <math.h>
#include <stddef.h>
#include "basis.h"



//*****************************************************************//
// define boys function
///////////////////////////////

double boys(double t);

//*****************************************************************//



//*****************************************************************//
// overlap matrix S
///////////////////////////////

double sElt(double p_oN, double q_oN, prim orbP[], prim orbQ[], int pLen, int qLen);

//*****************************************************************//



//*****************************************************************//
// kinetic energy matrix T
///////////////////////////////

double tElt(double p_oN, double q_oN, prim orbP[], prim orbQ[], int pLen, int qLen);

//*****************************************************************//



//*****************************************************************//
// nuclear potential energy matrix Vn
///////////////////////////////

double vnElt(double p_oN, double q_oN, prim orbP[], prim orbQ[], int pLen, int qLen,
  nuc nucleus);

//*****************************************************************//



//*****************************************************************//
// yeet eri
///////////////////////////////

double eriElt(double p_oN, double q_oN, double r_oN, double s_oN,
  prim orbP[], prim orbQ[], prim orbR[], prim orbS[],
  int pLen, int qLen, int rLen, int sLen);

//*****************************************************************//
