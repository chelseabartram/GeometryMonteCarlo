{
#include <sstream>
#include <fstream>
#include <string>
  gROOT->SetBatch(kTRUE);
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TFile f("first50.root");
  f->cd(0);
  TIter nextkey(f.GetListOfKeys());
  TKey *key;
  //  std::string namestring;
  std::string classstring;
  ostringstream ss("");
  int nHist=0;
  char namestring[16];
  char keyname[16];
  std::string keynamestring;
  std::string originalnamestring;
  std::string endTag;
  std::string endTag1;
  std::string endTag2;
  std::string endTagCut1;
  std::string endTagCut2;
  std::string endTagCut3;
  while ((key=(TKey*)nextkey())) {
    TObject *obj = key->ReadObj();
    printf("at key:%s, object
class:%s\n",key->GetName(),obj->ClassName());
    //    sprintf( keyname, key->GetName);
    keynamestring.assign(key->GetName(),0,7);
    std::cout << keynamestring << std::endl;
    if ( keynamestring.compare("cosXi_r") == 0 ) {
      std::cout << "WHAT THE HELL" << std::endl;
      c1->SetLogy();
    }
    else if (keynamestring.compare("hitTime") == 0 ) {
      c1->SetLogy();
    }
    else{
      c1->SetLogy(0);
    }
    originalnamestring.assign(key->GetName());
    std::cout << "WHAT THE HELL: " << originalnamestring << std::endl;
    //    endTag.assign(originalnamestring);//originalnamestring,originalnamestring.end()-6,3);

    endTag="";
    endTag1="";
    endTag2="";
    endTagCut1="";
    endTagCut2="";
    endTagCut3="";
    
    double stringlength = originalnamestring.length();
    if(stringlength>6)
      endTag.assign(originalnamestring,stringlength-6,6);
    std::cout << "ENDTAG: " << endTag << std::endl;

    if(stringlength>3)
      endTag1.assign(originalnamestring,stringlength-3,3);
    else{
      endTag1="";
    }

    if(stringlength>7)
      endTag2.assign(originalnamestring,stringlength-7,7);

    if(stringlength>12){
      endTagCut1.assign(originalnamestring,stringlength-12,12);
    }
    else{
      endTagCut1="";
    }

    
    std::cout << "ENDTAG: " << endTag << std::endl;
    std::cout << "ENDTAGCUT1: " << endTagCut1 << std::endl;
    //    namestring=obj->GetName();
    //    if(namestring=="cosXi_2")
    //c1->SetLogy();
    classstring=obj->ClassName();
    if(classstring.compare("TH2D")==0){


      if (endTag=="Vertex") {
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kGray);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTag1=="_wS") {
	std::cout << "WITH SOURCE HOLDER" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kGreen);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTag2=="APEX_RR") {
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kRed);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_RR_Cut1"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kOrange+10);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_RR_Cut2"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kOrange-3);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_RR_Cut3"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kRed-3);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTag2=="APEX_YM"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kBlue);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_YM_Cut1"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kBlue-7);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_YM_Cut2"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kCyan);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_YM_Cut3"){
	std::cout << "IT WORKS" << std::endl;
	std::cout << namestring << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kAzure-7);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }      

      else{
	std::cout << "ARG" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetTitleTextColor(kBlack);
	obj->Draw("hist");
      }
      
      obj->Draw("colz");

    }
    else{
      obj->Draw("hist");
      //      obj->GetXaxis();
      //      obj->GetXaxis()->SetBinLabel(5,"psi00_3gamma");
      if ( keynamestring.compare("PsSubType") == 0 ) {
	//	obj->GetXaxis()->SetBinLabel(5,"psi00_3gamma");
	/*      key->GetXaxis()->SetBinLabel(4,"psi10_2gamma");
		key->GetXaxis()->SetBinLabel(3,"psi00_2gamma");
		key->GetXaxis()->SetBinLabel(2,"oPs_unperturbed");
		key->GetXaxis()->SetBinLabel(1,"pPs_unperturbed");
		key->GetXaxis()->SetBinLabel(6,"psi10_3gamma");*/
      }

      if (endTag=="Vertex") {
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kGray);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTag1=="_wS") {
	std::cout << "WITH SOURCE HOLDER" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kGreen);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTag2=="APEX_RR") {
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kRed);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_RR_Cut1"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kOrange+10);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_RR_Cut2"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kOrange-3);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_RR_Cut3"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kRed-3);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTag2=="APEX_YM"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kBlue);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_YM_Cut1"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kBlue-3);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_YM_Cut2"){
	std::cout << "IT WORKS" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kCyan);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }
      else if (endTagCut1=="APEX_YM_Cut3"){
	std::cout << "IT WORKS" << std::endl;
	std::cout << namestring << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kAzure-7);
	gStyle->SetTitleTextColor(kWhite);
	obj->Draw("hist");
      }      

      else{
	std::cout << "ARG" << std::endl;
	gROOT->SetStyle("Plain");
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetTitleTextColor(kBlack);
	obj->Draw("hist");
      }


    }
    sprintf ( namestring, "%d", nHist );
    std::cout << namestring << std::endl;

    strcat(namestring,".png");
    std::cout << namestring << std::endl;

    c1->SaveAs(namestring,"png");
    nHist++;
    ss.clear();
  delete obj;
 }
  
  PsSubType_Vertex->GetXaxis()->SetBinLabel(1,"pPs_unperturbed");
  PsSubType_Vertex->GetXaxis()->SetBinLabel(2,"oPs_unperturbed");
  PsSubType_Vertex->GetXaxis()->SetBinLabel(3,"psi00_2gamma");
  PsSubType_Vertex->GetXaxis()->SetBinLabel(4,"psi10_2gamma");
  PsSubType_Vertex->GetXaxis()->SetBinLabel(5,"psi00_3gamma");
  PsSubType_Vertex->GetXaxis()->SetBinLabel(6,"psi10_3gamma");
  PsSubType_Vertex->Draw("hist");
  c1->SetLogy();
  c1->SaveAs("psSubType.png");



  asymmetryParameter = -(hAP1_Vertex->Integral(1,100) - hAP1_Vertex->Integral(101,200))/hAP1_Vertex->Integral();
  asymmetryParameter_Error = 1./TMath::Sqrt(hAP1_Vertex->Integral());
  asymmetryParameter_wS = -(hAP1_wS->Integral(1,100) - hAP1_wS->Integral(101,200))/hAP1_wS->Integral();
  asymmetryParameter_wS_Error = 1./TMath::Sqrt(hAP1_wS->Integral());
  asymmetryParameter_APEX_YM = -(hAP1_APEX_YM->Integral(1,100) - hAP1_APEX_YM->Integral(101,200))/hAP1_APEX_YM->Integral();
  asymmetryParameter_APEX_YM_Error = 1./TMath::Sqrt(hAP1_APEX_YM->Integral());
  asymmetryParameter_APEX_YM_Cut1 = -(hAP1_APEX_YM_Cut1->Integral(1,100) - hAP1_APEX_YM_Cut1->Integral(101,200))/hAP1_APEX_YM_Cut1->Integral();
  asymmetryParameter_APEX_YM_Cut1_Error = 1./TMath::Sqrt(hAP1_APEX_YM_Cut1->Integral());
  asymmetryParameter_APEX_YM_Cut2 = -(hAP1_APEX_YM_Cut2->Integral(1,100) - hAP1_APEX_YM_Cut2->Integral(101,200))/hAP1_APEX_YM_Cut2->Integral();
  asymmetryParameter_APEX_YM_Cut2_Error = 1./TMath::Sqrt(hAP1_APEX_YM_Cut2->Integral());
  asymmetryParameter_APEX_YM_Cut3 = -(hAP1_APEX_YM_Cut3->Integral(1,100) - hAP1_APEX_YM_Cut3->Integral(101,200))/hAP1_APEX_YM_Cut3->Integral();
  asymmetryParameter_APEX_YM_Cut3_Error = 1./TMath::Sqrt(hAP1_APEX_YM_Cut3->Integral());

  rh_asymmetryParameter = -(rh_asymmetry_Vertex->Integral(1,100) - rh_asymmetry_Vertex->Integral(101,200))/rh_asymmetry_Vertex->Integral();
  rh_asymmetryParameter_Error = 1./TMath::Sqrt(rh_asymmetry_Vertex->Integral());
  rh_asymmetryParameter_wS = -(rh_asymmetry_wS->Integral(1,100) - rh_asymmetry_wS->Integral(101,200))/rh_asymmetry_wS->Integral();
  rh_asymmetryParameter_wS_Error = 1./TMath::Sqrt(rh_asymmetry_wS->Integral());
  rh_asymmetryParameter_APEX_RR = -(rh_asymmetry_APEX_RR->Integral(1,100) - rh_asymmetry_APEX_RR->Integral(101,200))/rh_asymmetry_APEX_RR->Integral();
  rh_asymmetryParameter_APEX_RR_Error = 1./TMath::Sqrt(rh_asymmetry_APEX_RR->Integral());
  rh_asymmetryParameter_APEX_RR_Cut1 = -(rh_asymmetry_APEX_RR_Cut1->Integral(1,100) - rh_asymmetry_APEX_RR_Cut1->Integral(101,200))/rh_asymmetry_APEX_RR_Cut1->Integral();
  rh_asymmetryParameter_APEX_RR_Cut1_Error = 1./TMath::Sqrt(rh_asymmetry_APEX_RR_Cut1->Integral());
  rh_asymmetryParameter_APEX_RR_Cut2 = -(rh_asymmetry_APEX_RR_Cut2->Integral(1,100) - rh_asymmetry_APEX_RR_Cut2->Integral(101,200))/rh_asymmetry_APEX_RR_Cut2->Integral();
  rh_asymmetryParameter_APEX_RR_Cut2_Error = 1./TMath::Sqrt(rh_asymmetry_APEX_RR_Cut2->Integral());
  rh_asymmetryParameter_APEX_RR_Cut3 = -(rh_asymmetry_APEX_RR_Cut3->Integral(1,100) - rh_asymmetry_APEX_RR_Cut3->Integral(101,200))/rh_asymmetry_APEX_RR_Cut3->Integral();
  rh_asymmetryParameter_APEX_RR_Cut3_Error = 1./TMath::Sqrt(rh_asymmetry_APEX_RR_Cut3->Integral());


  CPAsymmetryParameter = 0.405597;
  CPAsymmetryParameter_WS = 0.382321;
  CPAsymmetryParameter_APEX = -0.0366094;
  CPAsymmetryParameter_APEX_CUT1 = -0.144826;
  CPAsymmetryParameter_APEX_CUT2 = -0.119088;
  CPAsymmetryParameter_APEX_CUT3 = -0.194133;
  
  C_CP_UPPER = (rh_asymmetryParameter + 1.645*rh_asymmetryParameter_Error)/(CPAsymmetryParameter);
  C_CP_LOWER = (rh_asymmetryParameter - 1.645*rh_asymmetryParameter_Error)/(CPAsymmetryParameter);
  C_CP = (C_CP_UPPER-C_CP_LOWER)/2.;
  C_CP_UPPER_WS = (rh_asymmetryParameter_wS + 1.645*rh_asymmetryParameter_wS_Error)/(CPAsymmetryParameter_WS);
  C_CP_LOWER_WS = (rh_asymmetryParameter_wS - 1.645*rh_asymmetryParameter_wS_Error)/(CPAsymmetryParameter_WS);
  C_CP_WS = (C_CP_UPPER_WS-C_CP_LOWER_WS)/2.;
  C_CP_UPPER_APEX = (rh_asymmetryParameter_APEX_RR + 1.645*rh_asymmetryParameter_APEX_RR_Error)/(CPAsymmetryParameter_APEX);
  C_CP_LOWER_APEX = (rh_asymmetryParameter_APEX_RR - 1.645*rh_asymmetryParameter_APEX_RR_Error)/(CPAsymmetryParameter_APEX);
  C_CP_APEX = (C_CP_UPPER_WS-C_CP_LOWER_WS)/2.;
  C_CP_UPPER_APEX_Cut1 = (rh_asymmetryParameter_APEX_RR_Cut1 + 1.645*rh_asymmetryParameter_APEX_RR_Cut1_Error)/(CPAsymmetryParameter_APEX_CUT1);
  C_CP_LOWER_APEX_Cut1 = (rh_asymmetryParameter_APEX_RR_Cut1 - 1.645*rh_asymmetryParameter_APEX_RR_Cut1_Error)/(CPAsymmetryParameter_APEX_CUT1);
  C_CP_APEX_CUT1 = (C_CP_UPPER_APEX_Cut1-C_CP_LOWER_APEX_Cut1)/2.;
  C_CP_UPPER_APEX_Cut2 = (rh_asymmetryParameter_APEX_RR_Cut2 + 1.645*rh_asymmetryParameter_APEX_RR_Cut2_Error)/(CPAsymmetryParameter_APEX_CUT2);
  C_CP_LOWER_APEX_Cut2 = (rh_asymmetryParameter_APEX_RR_Cut2 - 1.645*rh_asymmetryParameter_APEX_RR_Cut2_Error)/(CPAsymmetryParameter_APEX_CUT2);
  C_CP_APEX_CUT2 = (C_CP_UPPER_APEX_Cut2-C_CP_LOWER_APEX_Cut2)/2.;
  C_CP_UPPER_APEX_Cut3 = (rh_asymmetryParameter_APEX_RR_Cut3 + 1.645*rh_asymmetryParameter_APEX_RR_Cut3_Error)/(CPAsymmetryParameter_APEX_CUT3);
  C_CP_LOWER_APEX_Cut3 = (rh_asymmetryParameter_APEX_RR_Cut3 - 1.645*rh_asymmetryParameter_APEX_RR_Cut3_Error)/(CPAsymmetryParameter_APEX_CUT3);
  C_CP_APEX_CUT3 = (C_CP_UPPER_APEX_Cut3-C_CP_LOWER_APEX_Cut3)/2.;

  std::cout << "C_CP_UPPER (Vertex): " << C_CP_UPPER << std::endl;
  std::cout << "C_CP_LOWER (Vertex): " << C_CP_LOWER << std::endl;
  std::cout << "C_CP (Vertex): " << C_CP << std::endl;
  std::cout << "C_CP_UPPER (Outside Source Holder): " << C_CP_UPPER_WS << std::endl;
  std::cout << "C_CP_LOWER (Outside Source Holder): " << C_CP_LOWER_WS << std::endl;
  std::cout << "C_CP (Outside Source Holder: " << C_CP_WS << std::endl;
  std::cout << "C_CP_UPPER (APEX): " << C_CP_UPPER_APEX << std::endl;
  std::cout << "C_CP_LOWER (APEX): " << C_CP_LOWER_APEX << std::endl;
  std::cout << "C_CP (APEX): " << C_CP_APEX << std::endl;
  std::cout << "C_CP_UPPER (APEX Cut1): " << C_CP_UPPER_APEX_Cut1 << std::endl;
  std::cout << "C_CP_LOWER (APEX Cut1): " << C_CP_LOWER_APEX_Cut1 << std::endl;
  std::cout << "C_CP (APEX Cut1): " << C_CP_APEX_CUT1 << std::endl;
  std::cout << "C_CP_UPPER (APEX Cut2): " << C_CP_UPPER_APEX_Cut2 << std::endl;
  std::cout << "C_CP_LOWER (APEX Cut2): " << C_CP_LOWER_APEX_Cut2 << std::endl;
  std::cout << "C_CP (APEX Cut2): " << C_CP_APEX_CUT2 << std::endl;
  std::cout << "C_CP_UPPER (APEX Cut3): " << C_CP_UPPER_APEX_Cut3 << std::endl;
  std::cout << "C_CP_LOWER (APEX Cut3): " << C_CP_LOWER_APEX_Cut3 << std::endl;
  std::cout << "C_CP (APEX Cut3): " << C_CP_APEX_CUT3 << std::endl;




  
  std::cout << "Asymmetry Parameter at Vertex: " << asymmetryParameter << " +- " << 1./TMath::Sqrt(hAP1_Vertex->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter at Outside Source Holder: " << asymmetryParameter_wS << " +- " << 1./TMath::Sqrt(hAP1_wS->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter with APEX: " << asymmetryParameter_APEX_YM << " +- " << 1./TMath::Sqrt(hAP1_APEX_YM->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter with APEX Cut1: " << asymmetryParameter_APEX_YM_Cut1 << " +- " << 1./TMath::Sqrt(hAP1_APEX_YM_Cut1->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter with APEX Cut2: " << asymmetryParameter_APEX_YM_Cut2 << " +- " << 1./TMath::Sqrt(hAP1_APEX_YM_Cut2->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter with APEX Cut3: " << asymmetryParameter_APEX_YM_Cut3 << " +- " << 1./TMath::Sqrt(hAP1_APEX_YM_Cut3->Integral()) << std::endl;

  
  std::cout << "RH Asymmetry Parameter at Vertex: " << rh_asymmetryParameter << " +- " << 1./TMath::Sqrt(rh_asymmetry_Vertex->Integral()) << std::endl;
  std::cout << "RH Asymmetry Parameter at Outside Source Holder: " << rh_asymmetryParameter_wS << " +- " << 1./TMath::Sqrt(rh_asymmetry_wS->Integral()) << std::endl;
  std::cout << "RH Asymmetry Parameter with APEX: " << rh_asymmetryParameter_APEX_RR << " +- " << 1./TMath::Sqrt(rh_asymmetry_APEX_RR->Integral()) << std::endl;
    std::cout << "RH Asymmetry Parameter with APEX Cut1: " << rh_asymmetryParameter_APEX_RR_Cut1 << " +- " << 1./TMath::Sqrt(rh_asymmetry_APEX_RR_Cut1->Integral()) << std::endl;
    std::cout << "RH Asymmetry Parameter with APEX Cut2: " << rh_asymmetryParameter_APEX_RR_Cut2 << " +- " << 1./TMath::Sqrt(rh_asymmetry_APEX_RR_Cut2->Integral()) << std::endl;
    std::cout << "RH Asymmetry Parameter with APEX Cut3: " << rh_asymmetryParameter_APEX_RR_Cut3 << " +- " << 1./TMath::Sqrt(rh_asymmetry_APEX_RR_Cut3->Integral()) << std::endl;


    
}
