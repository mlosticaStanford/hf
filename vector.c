// functions for vector structure

#include "vector.h"


//****************************************************************//

vector retAddedVectors(vector vect1, vector vect2) {
  vect1.x += vect2.x; vect1.y += vect2.y; vect1.z += vect2.z;
  vect1.mag = retVectorMagnitude(vect1);

  return vect1;
}

//****************************************************************//

vector retSubtractedVectors(vector vect1, vector vect2) {

  vector vect3;

  vect3.x = vect1.x - vect2.x;
  vect3.y = vect1.y - vect2.y;
  vect3.z = vect1.z - vect2.z;

  vect3.mag = retVectorMagnitude(vect3);
           
  return vect3;
}

//****************************************************************//

vector scaleVector(vector vect, double scalar) {

  vector newVect;

  newVect.x = vect.x * scalar;
  newVect.y = vect.y * scalar;
  newVect.z = vect.z * scalar;
           
  return newVect;
}
//****************************************************************//

double retDotProduct(vector vect1, vector vect2) {
  return (vect1.x * vect2.x + vect1.y * vect2.y + vect1.z * vect2.z);
}

//****************************************************************//

double retVectorMagnitude(vector vect) {
  return sqrt(retDotProduct(vect, vect));
}

//****************************************************************//

vector newVec(double x, double y, double z) {
  vector v;

  v.x = x;
  v.y = y;
  v.z = z;
  v.mag = retVectorMagnitude(v);

  return v;
}
