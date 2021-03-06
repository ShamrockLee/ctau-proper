// example code to check the proper decay length of the chi2 in the long-lived
// extension of monoH

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>

using namespace std;
void xAna_monoH_ctau(std::string inputFile,std::string outputFile, bool debug=false){

  //get TTree from file ...
  TreeReader data(inputFile.data());
  
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};

  //TH1F* h_ctau=new TH1F("h_ctau","",100,0,1);
  //TH1F* h_ctau2=new TH1F("h_ctau2","",100000,0,10000);
  //TH1F* h_ctau_lab  = (TH1F*)h_ctau2->Clone("h_ctau_lab");
  //TH1F* h_chi2_gamma  = (TH1F*)h_ctau2->Clone("h_ctau_proper");
  //TH1F* h_chi2_gamma  = (TH1F*)h_ctau2->Clone("h_chi2_gamma");
  
  // Suitable for Mx2-1_Mx1-0p1
  /*
  TH1F* h_betatau_lab =new TH1F("h_betatau_lab","",1000,0,500);
  TH1F* h_ctau_lab =new TH1F("h_ctau_lab","", 1000, 0, 500);
  TH1F* h_ctau_proper =new TH1F("h_ctau_proper","", 1000, 0, 2);
  TH1F* h_chi2_beta =new TH1F("h_chi2_beta","", 1000, 0.9999, 1.0001);
  TH1F* h_chi2_gamma =new TH1F("h_chi2_gamma","", 1000, 0, 10000);
  */
  
  // Suitable for Mx2-150_Mx1-1
  TH1F* h_betatau_lab =new TH1F("h_betatau_lab","",500,0,5);
  TH1F* h_ctau_lab =new TH1F("h_ctau_lab","", 500, 0, 5);
  TH1F* h_ctau_proper =new TH1F("h_ctau_proper","", 500, 0, 2);
  TH1F* h_chi2_beta =new TH1F("h_chi2_beta","", 500, 0, 2);
  TH1F* h_chi2_gamma =new TH1F("h_chi2_gamma","", 500, 0, 20);

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal ++;

    // 0. check the generator-level information and make sure there is a Z->e+e-
    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genMomParId      = data.GetPtrInt("genMomParId");
    Int_t* genMo1      = data.GetPtrInt("genMo1");

    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");
    TClonesArray* genParVtx = (TClonesArray*) data.GetPtrTObject("genParVtx");


    // 1. find chi2 and chi1

    for(int ig=0; ig < nGenPar; ig++){
      nPass[0]++;
      if(genParId[ig]!=5000522)continue; // chi1 added by Pythia8
      nPass[1]++;
      int mo1 = genMo1[ig];
      if(mo1<0)continue;
      nPass[2]++;
      int momId=genMomParId[ig];
      if(abs(momId)!=18)continue;  // chi2
      nPass[3]++;
      //      TLorentzVector* chi1_p4 = (TLorentzVector*)genParP4->At(ig);      
      TLorentzVector* chi2_p4 = (TLorentzVector*)genParP4->At(mo1);

      // note, the position saved in CMSSW has a unit of cm
      // while in Pythia8 the setting has a unit of mm
      
      TVector3* chi1_vtx = (TVector3*)genParVtx->At(ig);
      TVector3* chi2_vtx = (TVector3*)genParVtx->At(mo1);

      TVector3 dist =  *chi1_vtx - *chi2_vtx;
      
      double betatau_lab = dist.Mag();
      double chi2_gamma = chi2_p4->Gamma();
      double chi2_beta = chi2_p4->Beta();
      double ctau_lab = betatau_lab / chi2_beta;
      double ctau_proper = ctau_lab / chi2_gamma;


      if(debug){
	cout << "ctau_lab = " << ctau_lab;
      }
      
      h_betatau_lab -> Fill(betatau_lab);
      h_ctau_lab -> Fill(ctau_lab);
      h_chi2_beta -> Fill(chi2_beta);
      h_chi2_gamma -> Fill(chi2_gamma);
      h_ctau_proper -> Fill(ctau_proper);
    }
  }
    

  Double_t ctauLabMaximum = h_ctau_lab->GetMaximum();
  Int_t ctauLabMaximumBin = h_ctau_lab->GetMaximumBin();
  Double_t ctauLabMinimum = h_ctau_lab->GetMinimum();
  Int_t ctauLabMinimumBin = h_ctau_lab->GetMinimumBin();
  Double_t ctauProperMaximum = h_ctau_proper->GetMaximum();
  Int_t ctauProperMaximumBin = h_ctau_proper->GetMaximumBin();
  Double_t ctauProperMinimum = h_ctau_proper->GetMinimum();
  Int_t ctauProperMinimumBin = h_ctau_proper->GetMinimumBin();
  Double_t chi2GammaMaximum = h_chi2_gamma->GetMaximum();
  Int_t chi2GammaMaximumBin = h_chi2_gamma->GetMaximumBin();
  Double_t chi2GammaMinimum = h_chi2_gamma->GetMinimum();
  Int_t chi2GammaMinimumBin = h_chi2_gamma->GetMinimumBin();
  
  
  TFile* outFile = new TFile(outputFile.data(),"recreate");
  h_betatau_lab->Write();
  h_ctau_lab->Write();
  h_ctau_proper->Write();
  h_chi2_beta->Write();
  h_chi2_gamma->Write();
  outFile->Close();
  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
    
  
  

}
