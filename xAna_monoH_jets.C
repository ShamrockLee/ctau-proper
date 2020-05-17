#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include "untuplizer.h"

void xAna_monoH_jets(std::string inputFile,std::string outputFile, bool toRecreateOutFile=true, bool debug=false) {
  bool isHerwigpp=(inputFile.find("herwigpp")!= std::string::npos);
  
  TreeReader data(inputFile.data());
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  
  TH1F* hFATjetPairP4M = new TH1F("hFATjetPairP4M", "Static Mass of FATjet Pairs", 50, 0, 500);
  TH1F* hFATjetPairP4Pt = new TH1F("hFATjetPairP4Pt", "Pt of FATjet Pairs", 100, 0, 1000);
  TH1F* hFATjetPairP4Rho = new TH1F("hFATjetPairP4Rho", "Rho of FATjet Pairs", 500, 0, 5000);
  TH1F* hFATjetPairP4Phi = new TH1F("hFATjetPairP4Phi", "Phi of FATjet Pairs", 200, -M_PI, M_PI);
  TH1F* hFATjetPairP4Theta = new TH1F("hFATjetPairP4Theta", "Theta of FATjet Pairs", 100, 0, M_PI);
  TH1F* hFATjetPairP4Eta = new TH1F("hFATjetPairP4Eta", "Pseudorapidity (Eta) of FATjet Pairs", 100, -5, 5);
  TH1F* hFATjetPairP4Rapidity = (TH1F*)hFATjetPairP4Eta->Clone("hFATjetPairP4Rapidity");
  hFATjetPairP4Rapidity->SetTitle("Rapidity of FATjet Pairs");
  
  
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genParSt      = data.GetPtrInt("genParSt");
    Int_t* genMomParId   = data.GetPtrInt("genMomParId");


    int genHIndex[2]={-1,-1};


    for(int ig=0; ig < nGenPar; ig++){


      if(genParId[ig]!=25)continue;


      if(isHerwigpp && genMomParId[ig]!=25)continue;


      if(genHIndex[0]<0)
          genHIndex[0]=ig;


      else if(genHIndex[1]<0)
          genHIndex[1]=ig;


    }    


    if(genHIndex[0]<0 || genHIndex[1]<0)continue;
    nPass[0]++;


    if(genHIndex[0]==genHIndex[1])continue;
    nPass[1]++;


    TLorentzVector genH_l4[2];
    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");


    for(int ih=0; ih<2; ih++)
    genH_l4[ih] = *((TLorentzVector*)genParP4->At(genHIndex[ih]));




    if(debug){
    std::cout << genHIndex[0] << "\t" << genHIndex[1] << std::endl;
    genH_l4[0].Print();
    genH_l4[1].Print();


    }
    int nFATJet         = data.GetInt("FATnJet");
    const int nJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");


    // check matching first


    TLorentzVector *fatjetPairP4 = nullptr;
    
    bool findAMatch=false;
    const float dRMax=0.4;
    int matchedHJetIndex[2]={-1,-1};
            
    for(int ij=0; ij<nJets; ij++)
    {
    TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);


    for(int jj=0; jj<nJets; jj++)
    {


      if(ij==jj)continue;
      TLorentzVector* thatJet = (TLorentzVector*)fatjetP4->At(jj);
      
      if(thisJet->DeltaR(genH_l4[0])<dRMax && 
      thatJet->DeltaR(genH_l4[1])<dRMax)
      {
      matchedHJetIndex[0]=ij;
      matchedHJetIndex[1]=jj;
      findAMatch=true;
      
      *fatjetPairP4 = *thisJet + *thatJet;
      
      break;
      }


      if(findAMatch)break;


    }	


    if(findAMatch)break;


    }


    if(!findAMatch)continue;
    
    hFATjetPairP4M->Fill(fatjetPairP4->M());
    hFATjetPairP4Pt->Fill(fatjetPairP4->Pt());
    hFATjetPairP4Rho->Fill(fatjetPairP4->Rho());
    hFATjetPairP4Phi->Fill(fatjetPairP4->Phi());
    hFATjetPairP4Theta->Fill(fatjetPairP4->Theta());
    hFATjetPairP4Eta->Fill(fatjetPairP4->Eta());
    hFATjetPairP4Rapidity->Fill(fatjetPairP4->Rapidity());
    
    nTotal++;
    
    
  }
  
  TFile* outFile = new TFile(outputFile.data(), toRecreateOutFile ? "recreate" : "update");
  hFATjetPairP4M->Write();
  hFATjetPairP4Pt->Write();
  hFATjetPairP4Rho->Write();
  hFATjetPairP4Phi->Write();
  hFATjetPairP4Theta->Write();
  hFATjetPairP4Eta->Write();
  hFATjetPairP4Rapidity->Write();

}
