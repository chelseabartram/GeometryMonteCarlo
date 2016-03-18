#include "GeometryCalculator.hh"
#include <cmath>
#include <iostream>
#include "TVector3.h"

double GeometryCalculator::CalculateIntersection(TVector3 origin, TVector3 k_i)
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
  double r = 2;

  double a = xd*xd+yd*yd;
  double b = 2*x0*xd+2*y0*yd;
  double c = x0*x0+y0*y0-r*r;


  double t = QuadraticEquation(a,b,c);

  return t;
}

double GeometryCalculator::QuadraticEquation(double a, double b, double c){

  double t_pos = (-b+sqrt(b*b-4*a*c))/(2*a);
  double t_neg = (-b-sqrt(b*b-4*a*c))/(2*a);

  return t_pos;
}