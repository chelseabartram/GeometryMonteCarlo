#ifndef GeometryCalculator_h
#define GeometryCalculator_h 1

#include <utility>
#include "TVector3.h"
class GeometryCalculator {

public:
  GeometryCalculator(){};
  std::pair<double,double> CalculateIntersection(TVector3,TVector3);
  
private:
  double PositiveQuadraticEquation(double, double, double);
  double NegativeQuadraticEquation(double, double, double);
  
};
#endif

