#include "GeometryCalculator.hh"
#include "PsGenerator.hh"
#include <iostream>
// C. Bartram
// 03/18/2016

int main() {

  TFile f1("test.root","RECREATE");
  TH1D *hAP1 = new TH1D("hAP1", "sin(2 theta) * cos(phi)", 200, -1, 1);
  
  PsGammaDecayGenerator *gen = new PsGammaDecayGenerator();
  // Generate gamma rays
  gen->Initialize();
  TVector3 k1, k2, k3;
  int m = 0;

  // Calculate the intersection
  GeometryCalculator *geometry = new GeometryCalculator();
  TVector3 decayVertex(0,0,0);

  double k1_x_intercept, k1_y_intercept, k1_z_intercept;
  double k2_x_intercept, k2_y_intercept, k2_z_intercept;
  double t1, t2;
  TVector3 k1APEX(0,0,0);
  TVector3 k2APEX(0,0,0);
  double theta, psi, phi;
  double xi1, xi2, sinxi1, cosxi1, cosxi2;
  TVector3 normal;
  double k1_mag, k2_mag;
  double sindelta;
  TVector3 sprojection;
  bool k1InRange, k2InRange;
  
  double cylinderLength = 55;
  for (int i = 0; i < 1000000 ; i++ ) {

    if ( i%1000 == 0 ) {    
    std::cout << i << std::endl;
    }
    
    gen->GenerateoPs(&k1,&k2,&k3,m);  
    t1 = geometry->CalculateIntersection(decayVertex,k1);
    t2 = geometry->CalculateIntersection(decayVertex,k2);
    // These are the intersection points for k1
    k1_x_intercept = decayVertex.x()+k1.x()*t1;
    k1_y_intercept = decayVertex.y()+k1.y()*t1;
    k1_z_intercept = decayVertex.z()+k1.z()*t1;

    if ( std::abs(k1_z_intercept) < cylinderLength/2. ) {    
      k1APEX.SetXYZ(k1_x_intercept, k1_y_intercept, k1_z_intercept);
      k1InRange = true;
    }
    else {
      std::cout << "throwaway event" << std::endl;
      std::cout << cylinderLength/2. << std::endl;
      std::cout << k1_z_intercept << std::endl;
      k1InRange = false;
    }
    // These are the intersection points for k2
    k2_x_intercept = decayVertex.x()+k2.x()*t2;
    k2_y_intercept = decayVertex.y()+k2.y()*t2;
    k2_z_intercept = decayVertex.z()+k2.z()*t2;

    if ( std::abs(k2_z_intercept) < cylinderLength/2. ) {     
      k2APEX.SetXYZ(k2_x_intercept, k2_y_intercept, k2_z_intercept);
      k2InRange = true;
    }
    else {
      std::cout << "throwaway event" << std::endl;
      std::cout << cylinderLength/2. << std::endl;
      std::cout << k2_z_intercept << std::endl;
      k2InRange = false
    }

    // Do the reconstruction
    if ( k1InRange && k2InRange ) {
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
      hAP1->Fill(TMath::Cos(phi) * TMath::Sin(2.*theta));
      hAP1->SetTitle("cos(phi)*sin(2theta)");
    }



  }
  // Returns a value t
  //  Calculate x, y, z on the cylinder
  std::cout << "wheeeee" << std::endl;

  hAP1->Write();
  f1.Close();
  
  delete gen;
  delete geometry;
  

}
