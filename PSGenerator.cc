// 
// The following class PsDecayGenerator generates momenta vectors for decays of oPs.
// It handles the following types of decays:
//
// 1. o-Ps Angular distribution for CP-violating decays.
// 2. o-Ps decays using Ore-Powell Energy distribution and Bernreuther angular distributions, as predicted by QED
//
// All units are in terms of electron rest masses., ie. energy=1 corresponds to one electron rest mass.
// The spin/magnetic field axis is in the +z direction in the coordinate system of the returned vectors.
// 
// Author: R. Henning
// Date: 3/5/2015
// 

#include <iostream>
#include <TF2.h>
#include <TF1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TROOT.h>
#include "PSGenerator.hh"
 
Double_t PsGammaDecayGenerator::oPsAngularDistribution_m_1(Double_t *theta, Double_t *p)
{
// Distribution of angle between normal of decay plane and spin axis for m=+-1. 

	Double_t c = TMath::Cos(*theta); 
	return 0.5*(3-c*c);
}

// ---------------------------------------------------------------------------------- //
 		
Double_t PsGammaDecayGenerator::oPsAngularDistribution_m_0(Double_t *theta, Double_t *p)
{
// Distribution of angle between normal of decay plane and spin axis for m=0. 

	Double_t c = TMath::Cos(*theta); 
	return 1.0+c*c;
}

Double_t PsGammaDecayGenerator::testAngularDistribution_fn(Double_t *theta, Double_t *p)
{
// Distribution of angle between normal of decay plane and spin axis for m=0. 

  Double_t c = 1/(TMath::Cos(2*(*theta)-1)); 
	return c;
}
		
// ---------------------------------------------------------------------------------- //
 
Double_t PsGammaDecayGenerator::oPsEnergyDistribution_OP(Double_t *k, Double_t *p)
{
// Energy distribution of single gamma-rays emitted during o-Ps decay assuming SM physics
// k is units of electron masses.
// From Ore, Powell, Phys. Rev. 75, 1696 (1949)

	Double_t a = 1 - *k;
	Double_t b = 2 - *k;
	return 2*(*k*a/(b*b) - 2*a*a/(b*b)*TMath::Log(a) + b/(*k) + 2*a/((*k)*(*k))*TMath::Log(a)); 
}

// ---------------------------------------------------------------------------------- //

Double_t PsGammaDecayGenerator::oPsEnergyDistribution_PS(Double_t *k, Double_t *p)
{
// Energy distribution of single gamma-rays emitted during o-Ps decay from phase space only
// k is units of electron masses.
// From Ore, Powell, Phys. Rev. 75, 1696 (1949)

	Double_t kk = (*k)*(*k);
	return kk*(6-6*(*k)+kk); 
}

// ---------------------------------------------------------------------------------- //
 
void PsGammaDecayGenerator::Initialize()
{
  //	std::cout << "Initializing Ps decay gamma-ray generator."<< std::endl;

	testAngularDistribution_m_0 = new TF1("testAngularDistribution_m_0", this, &PsGammaDecayGenerator::testAngularDistribution_fn, 0, TMath::Pi(), 0);
	foPsAngularDistribution_m_1 = new TF1("oPsAngularDistribution_m_1", this, &PsGammaDecayGenerator::oPsAngularDistribution_m_1, 0, TMath::Pi(), 0);
	foPsAngularDistribution_m_0 = new TF1("oPsAngularDistribution_m_0", this, &PsGammaDecayGenerator::oPsAngularDistribution_m_0, 0, TMath::Pi(), 0);
	foPsAngularDistribution_CP_theta_phi = new TF2("oPsAngularDistribution_CP_theta_phi", "1+sin(2*y)*cos(x)", -TMath::Pi(), TMath::Pi(), 0, TMath::Pi());
	foPsEnergyDistribution_OP = new TF1("oPsEnergyDistribution_OP", this, &PsGammaDecayGenerator::oPsEnergyDistribution_OP, 0.001, 0.999, 0); //Doesn't like 0,1
	foPsEnergyDistribution_PS = new TF1("oPsEnergyDistribution_PS", this, &PsGammaDecayGenerator::oPsEnergyDistribution_PS, 0., 1., 0);

	foPsAngularDistribution_CP_theta_phi->SetNpx(200);
	foPsAngularDistribution_CP_theta_phi->SetNpy(100);

	// The following Write() statements should be removed in production code. 
	//	foPsAngularDistribution_m_0->Write();
	//	foPsAngularDistribution_m_1->Write();
	//	foPsAngularDistribution_CP_theta_phi->Write();
	//	foPsEnergyDistribution_OP->Write();
	//	foPsEnergyDistribution_PS->Write();
	
    delete gRandom;
	gRandom = new TRandom3(0); // Use UUID as seed. 
	//    std::cout << "TRandom3 seed: " << gRandom->GetSeed() << std::endl;
	
	std::cout << "Done initializing Ps decay gamma-ray generator."<< std::endl;
}

// ---------------------------------------------------------------------------------- //

Int_t PsGammaDecayGenerator::GenerateoPs(TVector3 *k1, TVector3 *k2, TVector3 *k3, Int_t m)
{

  Double_t E1, E2;
  do {
    E1 = foPsEnergyDistribution_OP->GetRandom(0.0001, 0.9999);
    E2 = foPsEnergyDistribution_OP->GetRandom(0.0001, 0.9999);
  } while (E1 + E2 < 1.0);
  DecayPlaneKinematics(k1, k2, k3, E1, E2 );

	
	Double_t angle = (m == 0) ? foPsAngularDistribution_m_0->GetRandom() : foPsAngularDistribution_m_1->GetRandom();
	//	Double_t angle = testAngularDistribution_m_0->GetRandom();
	//	Double_t nu=gRandom->Rndm();
	//	angle=TMath::ACos(2*nu-1);
	k1->RotateX(angle);
	k2->RotateX(angle);
	k3->RotateX(angle);
	
	angle = gRandom->Rndm()*TMath::Pi()*2;
	k1->RotateZ(angle);
	k2->RotateZ(angle);
	k3->RotateZ(angle);
	
	return 1;
}


// ---------------------------------------------------------------------------------- //
// ---------------------------------------------------------------------------------- //

Int_t PsGammaDecayGenerator::GenerateoPs_CP(TVector3 *k1, TVector3 *k2, TVector3 *k3)
{
    
	// The energy distribution of the photons is model-dependent, hence we don't know what
	// it will be. For now we just pick them using phase
	// space calculations from Ore-Powell (ie. the matrix element is a constant).
	//
	
	DecayPlaneKinematics(k1, k2, k3, foPsEnergyDistribution_PS->GetRandom(),
                         foPsEnergyDistribution_PS->GetRandom());

	TVector3 zAxis(0,0,1);
	
	// Find theta and phi. Psi is determined by kinematics. 
	Double_t theta, phi;
	foPsAngularDistribution_CP_theta_phi->GetRandom2(phi, theta);
	
	// Find rotation axis to rotate decay plane normal vector by. 
	TVector3 rotationAxis = k1->Cross(zAxis);
	
	// Rotate vectors on decay plane to proper phi value.
	k1->RotateZ(phi);
	k2->RotateZ(phi);
	k3->RotateZ(phi);
	
	// Rotate entire decay plane. Now the angle between the projection of the z-axis on the decay plane
	// and k1 is phi. 
	// k1 x k2 at this point may be pointing in the -z direction. Correct for this. 
	
	if(k1->Cross(*k2) * zAxis < 0) 
		theta = TMath::Pi() - theta;
	
	k1->Rotate(theta, rotationAxis);
	k2->Rotate(theta, rotationAxis);
	k3->Rotate(theta, rotationAxis);
	
	return 1;
}


// ---------------------------------------------------------------------------------- //

void PsGammaDecayGenerator::DecayPlaneKinematics(TVector3 *k1, TVector3 *k2, TVector3 *k3, Double_t x1, Double_t x2)
{
// Generates 3 random momentum vectors on the z-plane. The vectors are kinetically constrained so that they sum to zero.
// Two of the vectors will have magnitude x1 and x2, with the third being 2 - x1 - x2
// k1 and k2 will have the largest and second largest momenta, respectively. 
	
	Double_t k1_mag, k2_mag, k3_mag;
	Double_t x3 = 2.0 - x1 - x2;
	
// Sort so that k1_mag > k2_mag > k3_mag
	
	if(x1>x2)
		if(x3>x1) {
			k1_mag = x3;
			k2_mag = x1;
			k3_mag = x2;
		} else if(x2 > x3) {
			k1_mag = x1;
			k2_mag = x2;
			k3_mag = x3;
		} else {
			k1_mag = x1;
			k2_mag = x3;
			k3_mag = x2;
			}
	else if(x3 > x2) {
			k1_mag = x3;
			k2_mag = x2;
			k3_mag = x1;
		} else if(x1 > x3) {
			k1_mag = x2;
			k2_mag = x1;
			k3_mag = x3;
		} else {
			k1_mag = x2;
			k2_mag = x3;
			k3_mag = x1;
			}

// Create three vectors on z-plane using kinematic constraints. k3 points along the -x axis for simplification.
// Since we have ranked the vectors, the equations below will generate the directions with a certain handedness, ie. k1, k2 and k3
// will always be in the same rotational order when viewed from, say the +z axis. We break this bias by 
// randomizing the sign of k1_y.
    
	Double_t k1_x = -0.5*((k2_mag*k2_mag - k1_mag*k1_mag)/k3_mag - k3_mag);
	if(k1_x>k1_mag){	  
	  /*	  std::cout << "BAD" << std::endl;
	  std::cout << "k1_x: " << k1_x << std::endl;
	  std::cout << "k1_mag: " << k1_mag << std::endl;
	  std::cout << "k2_mag: " << k2_mag << std::endl;
	  std::cout << "k3_mag: " << k3_mag << std::endl;*/
	}
	else {
	  /*	  std::cout << "GOOD" << std::endl;
	  std::cout << "k1_x: " << k1_x << std::endl;
	  std::cout << "k1_mag: " << k1_mag << std::endl;
	  std::cout << "k2_mag: " << k2_mag << std::endl;
	  std::cout << "k3_mag: " << k3_mag << std::endl;*/
	}

	Double_t k1_y = TMath::Sqrt(k1_mag*k1_mag - k1_x*k1_x) * ((gRandom->Rndm() > 0.5) ? -1 : 1);
	/*	if ((k1_mag*k1_mag - k1_x*k1_x)<0) {
	  std::cout << "OH NO" << std::endl;
	  std::cout << "VAL: " << TMath::Sqrt(k1_mag*k1_mag - k1_x*k1_x) * ((gRandom->Rndm() > 0.5) ? -1 : 1) << std::endl;
	}
	else {
	  std::cout << "CONFUSED" << std::endl;
	  std::cout << "VAL: " << TMath::Sqrt(k1_mag*k1_mag - k1_x*k1_x) * ((gRandom->Rndm() > 0.5) ? -1 : 1) << std::endl;
	  }*/
	//	std::cout << "VAL: " << TMath::Sqrt(k1_mag*k1_mag - k1_x*k1_x) * ((gRandom->Rndm() > 0.5) ? -1 : 1) << std::endl;
    
	k1->SetXYZ(k1_x, k1_y, 0.);
	//	std::cout << "X COMP: " <<  k1->Y() << std::endl;
	//	std::cout << "Y COMP: " <<  k1->Y() << std::endl;
	//	std::cout << "k1Mag: " << k1->Mag() << std::endl;
	k2->SetXYZ(k3_mag - k1_x, -k1_y, 0.);
	k3->SetXYZ(-k3_mag, 0., 0.);

	
// Rotate by random angle around z-axis. k3 no longer points in -x direction. 
	Double_t angle = gRandom->Rndm()*TMath::Pi()*2;
	k1->RotateZ(angle);
	k2->RotateZ(angle);
	k3->RotateZ(angle);
}

// ---------------------------------------------------------------------------------- //
// ---------------------------------------------------------------------------------- //

		
