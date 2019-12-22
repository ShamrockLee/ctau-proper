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

  TH1F* h_ctau=new TH1F("h_ctau","",100,0,1);
  TH1F* h_ctau2=new TH1F("h_ctau2","",100000,0,10000);
  TH1F* h_ctau_lab  = (TH1F*)h_ctau2->Clone("h_ctau_lab");
  TH1F* h_ctau_proper  = (TH1F*)h_ctau2->Clone("h_ctau_proper");
  TH1F* h_chi2_gamma  = (TH1F*)h_ctau2->Clone("h_chi2_gamma");

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
      
      double ctau_lab = dist.Mag();
      double chi2_gamma = chi2_p4->Gamma();
      double ctau_proper = ctau_lab / (double)(chi2_p4->Gamma());


      if(debug){
	cout << "ctau_lab = " << ctau_lab;
      }
      
      h_ctau_lab -> Fill(ctau_lab);
      h_chi2_gamma -> Fill(chi2_gamma);
      h_ctau_proper -> Fill(ctau_proper);
    }
  }
    


  TFile* outFile = new TFile(outputFile.data(),"recreate");
  h_ctau_lab->Write();
  h_ctau_proper->Write();
  h_chi2_gamma->Write();
  outFile->Close();
  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
    
  
  

}
