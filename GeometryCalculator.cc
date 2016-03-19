// The following code calculates the intersection between a ray and a cylinder
// The ray is written in parameterized form
// The intersection points are characterized by t_pos and t_neg, which are the solutions to the quadratic equation
// t_pos and t_neg are returned to the main code, where the actual intersection point is determined and evaluated to be either inside or outside the z bounds of the APEX cylinder
// The radius of the APEX array (r) is 24.0*cm
// Author: C. Bartram
// Date 3/18/2016

#include "GeometryCalculator.hh"
#include <cmath>
#include <iostream>
#include "TVector3.h"
#include <utility>

std::pair<double,double> GeometryCalculator::CalculateIntersection(TVector3 origin, TVector3 k_i)
{
  // Get the origin of the decay
  // The origin is basically whereever the source holder is
  // Return the vector for k1, k2, k3

  // x(t)= x0+xa*t
  // origin is x0
  // offset is xd (k1.x())


  // Equation for the APEX cylinder if x*x+y*y=r*r where r=24cm

  // Calculate a b and c
  double xd = k_i.x();
  double yd = k_i.y();
  double x0 = origin.x();
  double y0 = origin.y();
  double r = 24.0;

  double a = xd*xd+yd*yd;
  double b = 2*x0*xd+2*y0*yd;
  double c = x0*x0+y0*y0-r*r;

  std::pair <double,double> t;
  double t_pos = PositiveQuadraticEquation(a,b,c);
  double t_neg = NegativeQuadraticEquation(a,b,c);
  t = std::make_pair(t_pos,t_neg);
  
  return t;
}

double GeometryCalculator::PositiveQuadraticEquation(double a, double b, double c){

  double t_pos = (-b+sqrt(b*b-4*a*c))/(2*a);

  return t_pos;
}

double GeometryCalculator::NegativeQuadraticEquation(double a, double b, double c){

  double t_neg = (-b-sqrt(b*b-4*a*c))/(2*a);

  return t_neg;
}

