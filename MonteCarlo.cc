// The following code simulates the reconstruction of gamma ray events soley with geometry
// No physics processes are simulated--those can be simulated with the corresponding Geant code
// Units are in cm

// Author: C. Bartram
// Date: 3/18/2016


#include "GeometryCalculator.hh"
#include "MonteCarlo.hh"
#include "PsGenerator.hh"
#include <iostream>
#include <utility>
#include <math.h>
#define _USE_MATH_DEFINES

int main() {

  // Should be user-defined stuff
  //====================================================================
  // Define the decay vertex (location of the source holder--user input)
  TVector3 decayVertex(0,0,0);

  // Define a rotation of the magnetic field (this can be due to the magnet being rotated
  // or an inhomogeneity in the field (the former is more likely)
  // Units are in radians
  double rotationAboutX = 5*M_PI/180.;
  double rotationAboutY = 0*M_PI/180.;
  double rotationAboutZ = 0*M_PI/180.;

  int maximumEventNum = 100000000;
  //====================================================================
  
  // Define the output file
  TFile f1("test.root","RECREATE");
  // Define the asymmetry histogram
  TH1D *hAP1 = new TH1D("hAP1", "sin(2 theta) * cos(phi)", 200, -1, 1);

  // Reyco's nifty PS Generator
  PsGammaDecayGenerator *gen = new PsGammaDecayGenerator();
  // Initialize PS Generator
  gen->Initialize();
  // Define gamma rays k1, k2, k3
  TVector3 k1, k2, k3;
  int m;

  // GeometryCalculator class calculates intersection points of gamma rays with cylinder (APEX)
  GeometryCalculator *geometry = new GeometryCalculator();

  // Initialization of variables
  MonteCarlo reconstruction;
  double k1_x_intercept, k1_y_intercept, k1_z_intercept;
  double k1_z_intercept_pos, k1_z_intercept_neg;
  double k2_x_intercept, k2_y_intercept, k2_z_intercept;
  double k2_z_intercept_pos, k2_z_intercept_neg;
  std::pair<double,double> t1, t2;
  TVector3 k1APEX(0,0,0);
  TVector3 k2APEX(0,0,0);
  double theta, psi, phi;
  double xi1, xi2, sinxi1, cosxi1, cosxi2;
  TVector3 normal;
  double k1_mag, k2_mag;
  double sindelta;
  TVector3 sprojection;
  bool k1InRange, k2InRange;
  double t;
  bool t1first, t1second;
  bool t2first, t2second;
  double randomNumber;
  TRandom* ran = new TRandom();
  ran->SetSeed(0);

  // APEX Bar Length (length of cylindrical array: 55 cm)
  double cylinderLength = 55.0;
  
  for (int i = 0; i < maximumEventNum ; i++ ) {

    t1first = false;
    t1second = false;
    t2first = false;
    t2second = false;
    
    if ( i%100000 == 0 ) {    
    std::cout << i << std::endl;
    }

    // Need the appropriate ratios of z quantum numbers for generating o-PS
    // m = 0 happens 1/3 of the time, m = +1 happens 1/3, m = -1 happens 1/3
    // Assumes no B Field
    randomNumber = ran->Uniform(0,1);
    if ( randomNumber < 0.33 ){
      m = 0;
      gen->GenerateoPs(&k1,&k2,&k3,m);
    }
    if ( randomNumber > 0.66 ){
      m = 1;
      gen->GenerateoPs(&k1,&k2,&k3,m);
    }
    else {
      m = -1;
      gen->GenerateoPs(&k1,&k2,&k3,m);    
    }

    k1.RotateX(rotationAboutX);
    k1.RotateY(rotationAboutY);
    k1.RotateZ(rotationAboutZ);
    k2.RotateX(rotationAboutX);
    k2.RotateY(rotationAboutY);
    k2.RotateZ(rotationAboutZ);

    // Calculates the intersection of k1 with APEX
    // Returns pair of solutions (t_pos, t_neg)
    t1 = geometry->CalculateIntersection(decayVertex,k1);

    // These are the intersection points for k1
    k1_z_intercept_pos = decayVertex.z()+k1.z()*t1.first;
    k1_z_intercept_neg = decayVertex.z()+k1.z()*t1.second;

    // Check to see if both are within the bounds of the cylinder (APEX)
    //=================================================================
    if ( std::abs(k1_z_intercept_pos) < cylinderLength/2. )
      t1first = true;

    if ( std::abs(k1_z_intercept_neg) < cylinderLength/2. )
      t1second = true;
    //=================================================================

    // If only t1_pos is within bounds, use that one
    if ( t1first==true && t1second==false ){
      k1_x_intercept = decayVertex.x()+k1.x()*t1.first;
      k1_y_intercept = decayVertex.y()+k1.y()*t1.first;
      k1_z_intercept = decayVertex.z()+k1.z()*t1.first;
      k1APEX.SetXYZ(k1_x_intercept, k1_y_intercept, k1_z_intercept);
      k1InRange = true;
    }
    
    // If only t1_neg is within bounds, use that one
    if ( t1first==false && t1second==true ){
      k1_x_intercept = decayVertex.x()+k1.x()*t1.second;
      k1_y_intercept = decayVertex.y()+k1.y()*t1.second;
      k1_z_intercept = decayVertex.z()+k1.z()*t1.second;
      k1APEX.SetXYZ(k1_x_intercept, k1_y_intercept, k1_z_intercept);
      k1InRange = true;
    }

    // If both are in bounds, use the one with minimum value for t1
    if ( t1first==true && t1second==true ) {
      // Use the one with minimum t
      if ( t1.first > t1.second ) {
	k1_x_intercept = decayVertex.x()+k1.x()*t1.first;
	k1_y_intercept = decayVertex.y()+k1.y()*t1.first;
	k1_z_intercept = decayVertex.z()+k1.z()*t1.first;
	k1APEX.SetXYZ(k1_x_intercept, k1_y_intercept, k1_z_intercept);
	k1InRange = true;
      }
      else{
	k1_x_intercept = decayVertex.x()+k1.x()*t1.second;
	k1_y_intercept = decayVertex.y()+k1.y()*t1.second;
	k1_z_intercept = decayVertex.z()+k1.z()*t1.second;
	k1APEX.SetXYZ(k1_x_intercept, k1_y_intercept, k1_z_intercept);
	k1InRange = true;
      }
    }

    // If neither is within bounds, throw it away
    //    if ( t1first==false && t1second==false ) {
      // This is just here because I was curious how many events I was throwing away.
      // Might want to histogram this as diagnostic
    //    }

    //=============================================
    // Do all of that over again, this time for k2
    //=============================================

    // Calculates the intersection of k1 with APEX
    // Returns pair of solutions (t_pos, t_neg)
    t2 = geometry->CalculateIntersection(decayVertex,k2);

    // These are the intersection points for k2
    k2_z_intercept_pos = decayVertex.z()+k2.z()*t2.first;
    k2_z_intercept_neg = decayVertex.z()+k2.z()*t2.second;

    // Check to see if both are within the bounds of the cylinder (APEX)
    //=================================================================
    if ( std::abs(k2_z_intercept_pos) < cylinderLength/2. )
      t2first = true;

    if ( std::abs(k2_z_intercept_neg) < cylinderLength/2. )
      t2second = true;
    //=================================================================

    // If only t2_pos is within bounds, use that one
    if ( t2first==true && t2second==false ){
      k2_x_intercept = decayVertex.x()+k2.x()*t2.first;
      k2_y_intercept = decayVertex.y()+k2.y()*t2.first;
      k2_z_intercept = decayVertex.z()+k2.z()*t2.first;
      k2APEX.SetXYZ(k2_x_intercept, k2_y_intercept, k2_z_intercept);
      k2InRange = true;
    }

    // If only t2_neg is within bounds, use that one
    if ( t2first==false && t2second==true ){
      k2_x_intercept = decayVertex.x()+k2.x()*t1.second;
      k2_y_intercept = decayVertex.y()+k2.y()*t1.second;
      k2_z_intercept = decayVertex.z()+k2.z()*t1.second;
      k2APEX.SetXYZ(k2_x_intercept, k2_y_intercept, k2_z_intercept);
      k2InRange = true;
    }

    // If both are in bounds, use the one with minimum value for t2
    if ( t2first==true && t2second==true ) {
      // Use the one with minimum t
      if ( t2.first > t2.second ) {
	k2_x_intercept = decayVertex.x()+k2.x()*t2.first;
	k2_y_intercept = decayVertex.y()+k2.y()*t2.first;
	k2_z_intercept = decayVertex.z()+k2.z()*t2.first;
	k2APEX.SetXYZ(k2_x_intercept, k2_y_intercept, k2_z_intercept);
	k2InRange = true;
      }
      else{
	k2_x_intercept = decayVertex.x()+k2.x()*t2.second;
	k2_y_intercept = decayVertex.y()+k2.y()*t2.second;
	k2_z_intercept = decayVertex.z()+k2.z()*t2.second;
	k2APEX.SetXYZ(k2_x_intercept, k2_y_intercept, k2_z_intercept);
	k2InRange = true;
      }
    }
    
    // If neither is within bounds, throw it away
    //    if ( t2first==false && t2second==false ) {
      // Again, this could be used as diagnostic...
      //std::cout << "throwaway event" << std::endl;
    //    }

    // If we have set both the k1APEX and k2APEX vectors, do the reconstruction
    if ( k1InRange && k2InRange ) 
      reconstruction.DoReconstruction(k1APEX,k2APEX,hAP1);

  }

  std::cout << "Run over. Writing to file now." << std::endl;
  hAP1->Write();
  f1.Close();
  
  delete gen;
  delete geometry;
  delete ran;
  

}

void MonteCarlo::DoReconstruction(TVector3 k1APEX, TVector3 k2APEX, TH1D* histPtr){

  double theta, psi, phi;
  double xi1, xi2, sinxi1, cosxi1, cosxi2;
  TVector3 normal;
  double k1_mag, k2_mag;
  double sindelta;
  TVector3 sprojection;

  normal = k1APEX.Cross(k2APEX);
  theta = TMath::ACos(normal.CosTheta());
  psi = k1APEX.Angle(k2APEX);
  k1_mag = k1APEX.Mag(); k2_mag = k2APEX.Mag(); 
  sinxi1 = k1APEX.Z()/k1_mag;
  cosxi1 = TMath::Sqrt(1-k1APEX.Z()*k1APEX.Z()/k1_mag/k1_mag);
  xi1 = TMath::ACos(cosxi1);
  cosxi2 = TMath::Sqrt(1-k2APEX.Z()*k2APEX.Z()/k2_mag/k2_mag);
  xi2 = TMath::ACos(cosxi2);
  sindelta = (k1APEX.X()*k2APEX.Y()-k1APEX.Y()*k2APEX.X())/cosxi1/cosxi2/k1_mag/k2_mag;
  
  sprojection = normal;
  sprojection.Rotate(TMath::Pi()*0.5, normal.Cross(TVector3(0,0,1))); // projection of z-axis (S) onto decay plane.
  phi = k1APEX.Angle(sprojection);
  histPtr->Fill(TMath::Cos(phi) * TMath::Sin(2.*theta));
  histPtr->SetTitle("cos(phi)*sin(2theta)");
  
}

