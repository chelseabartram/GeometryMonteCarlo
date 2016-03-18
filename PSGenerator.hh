#ifndef PSGenerator_h
#define PSGenerator_h 1
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

class PsGammaDecayGenerator {
	
	public:
		PsGammaDecayGenerator(){;}
		void Initialize();
		
		Double_t oPsAngularDistribution_m_1(Double_t *theta, Double_t *p); // See descriptions below
		Double_t oPsAngularDistribution_m_0(Double_t *theta, Double_t *p);
		Double_t testAngularDistribution_fn(Double_t *theta, Double_t *p);
		Double_t oPsEnergyDistribution_OP(Double_t *k, Double_t *p);
		Double_t oPsEnergyDistribution_PS(Double_t *k, Double_t *p);

		// Generate random momentum vectors for the given decay modes
		Int_t GenerateoPs(TVector3 *k1, TVector3 *k2, TVector3 *k3, Int_t m);
		Int_t GenerateoPs_CP(TVector3 *k1, TVector3 *k2, TVector3 *k3);
		
	private:
	
		void DecayPlaneKinematics(TVector3 *k1, TVector3 *k2, TVector3 *k3, Double_t E1, Double_t E2);
		
	// This set of functions will be defined during initialization. They are not defined when the 
	// user requests the vertices to limit CPU waste on dynamic memory allocation.
	//

		TF1 *foPsAngularDistribution_m_0;			// Angular distribution of o-Ps decay from QED for m=0.
		TF1 *foPsAngularDistribution_m_1;			// Angular distribution of o-Ps decay from QED for m=+-1.
		TF2 *foPsAngularDistribution_CP_theta_phi;	// Angular distribution for CP-violating decays. Using convention
													// from Yamazaki, PRL 104, 083401 (2010)

		TF1 *testAngularDistribution_m_0;			// Angular distribution of o-Ps decay from QED for m=0.
		TF1 *foPsEnergyDistribution_OP;				// Energy distribution for single gammas in o-Ps decay per Ore-Powell QED
		TF1 *foPsEnergyDistribution_PS;				// Energy distribution for single gammas in o-Ps decay from phase space only.
};
#endif
