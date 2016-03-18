#ifndef GeometryCalculator_h
#define GeometryCalculator_h 1

#include "TVector3.h"
class GeometryCalculator {

public:
  GeometryCalculator(){};
  double CalculateIntersection(TVector3,TVector3);
private:
  double QuadraticEquation(double, double, double);
  
};
#endif

