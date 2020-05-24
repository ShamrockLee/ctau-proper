#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include "untuplizer.h"

void xAna_monoH_jets(std::string inputFile,std::string outputFile, bool toRecreateOutFile=true, bool debug=false) {
  // bool isHerwigpp=(inputFile.find("herwigpp")!= std::string::npos);
  
  TreeReader data(inputFile.data());
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  
  TH1F* hTHINNJetMatched = new TH1F("hTHINNJetMatched", "number of matched thin jets (d or dbar from chi2)", 8, 0-0.5f, 8-0.5f); 
  // TH1F* hTHINjetPairP4M = new TH1F("hTHINjetPairP4M", "Static Mass of THINjet Pairs", 50, 0, 500);
  // TH1F* hTHINjetPairP4Pt = new TH1F("hTHINjetPairP4Pt", "Pt of THINjet Pairs", 100, 0, 1000);
  // TH1F* hTHINjetPairP4Rho = new TH1F("hTHINjetPairP4Rho", "Rho of THINjet Pairs", 500, 0, 5000);
  // TH1F* hTHINjetPairP4Phi = new TH1F("hTHINjetPairP4Phi", "Phi of THINjet Pairs", 200, -M_PI, M_PI);
  // TH1F* hTHINjetPairP4Theta = new TH1F("hTHINjetPairP4Theta", "Theta of THINjet Pairs", 100, 0, M_PI);
  // TH1F* hTHINjetPairP4Eta = new TH1F("hTHINjetPairP4Eta", "Pseudorapidity (Eta) of THINjet Pairs", 100, -5, 5);
  // TH1F* hTHINjetPairP4Rapidity = (TH1F*)hTHINjetPairP4Eta->Clone("hTHINjetPairP4Rapidity");
  // hTHINjetPairP4Rapidity->SetTitle("Rapidity of THINjet Pairs");
  
  
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++) {
    data.GetEntry(jEntry);
    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genParSt      = data.GetPtrInt("genParSt");
    Int_t* genMomParId   = data.GetPtrInt("genMomParId");


    // int genHIndex[2]={-1,-1};
    std::vector<Int_t> vIndexesJetMatched;
    vIndexesJetMatched.clear();


    for(int ig=0; ig < nGenPar; ig++){


      if(abs(genParId[ig])!=1)continue;


      if(abs(genMomParId[ig])!=18)continue;


      // if(genHIndex[0]<0)
      //     genHIndex[0]=ig;


      // else if(genHIndex[1]<0)
      //     genHIndex[1]=ig;
      
      vIndexesJetMatched.push_back(ig);


    }    


    // if(genHIndex[0]<0 || genHIndex[1]<0)continue;
    // nPass[0]++;


    // if(genHIndex[0]==genHIndex[1])continue;
    // nPass[1]++;


    // TLorentzVector genH_l4[2];
    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");


    // for(int ih=0; ih<2; ih++)
    // genH_l4[ih] = *((TLorentzVector*)genParP4->At(genHIndex[ih]));




    if(debug){
    // std::cout << genHIndex[0] << "\t" << genHIndex[1] << std::endl;
    // genH_l4[0].Print();
    // genH_l4[1].Print();
    std::cout << vIndexesJetMatched.size() << " jets: ";
    for (Int_t ig: vIndexesJetMatched) std::cout << ig << " ";
    std::cout << std::endl;
    
    for (Int_t ig: vIndexesJetMatched) genParP4->At(ig)->Print();
    
    std::cout << std::endl;

    }
    int nTHINJet         = data.GetInt("THINnJet");
    const int nJets=nTHINJet;
    TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");


    // check matching first


    TLorentzVector *thinjetPairP4 = nullptr;
    
    // bool findAMatch=false;
    // const float dRMax=0.4;
    // int matchedHJetIndex[2]={-1,-1};
            
    for(int ij=0; ij<nJets; ij++)
    {
    // TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);


    // for(int jj=0; jj<nJets; jj++)
    // {


      // if(ij==jj)continue;
      // TLorentzVector* thatJet = (TLorentzVector*)thinjetP4->At(jj);
      
      // if(thisJet->DeltaR(genH_l4[0])<dRMax && 
      // thatJet->DeltaR(genH_l4[1])<dRMax)
      // {
      // matchedHJetIndex[0]=ij;
      // matchedHJetIndex[1]=jj;
      // findAMatch=true;
      
      // *thinjetPairP4 = *thisJet + *thatJet;
      
      // break;
      // }


      // if(findAMatch)break;


    // }	


    // if(findAMatch)break;


    }


    // if(!findAMatch)continue;
    
    hTHINNJetMatched->Fill(vIndexesJetMatched.size());
    // hTHINjetPairP4M->Fill(thinjetPairP4->M());
    // hTHINjetPairP4Pt->Fill(thinjetPairP4->Pt());
    // hTHINjetPairP4Rho->Fill(thinjetPairP4->Rho());
    // hTHINjetPairP4Phi->Fill(thinjetPairP4->Phi());
    // hTHINjetPairP4Theta->Fill(thinjetPairP4->Theta());
    // hTHINjetPairP4Eta->Fill(thinjetPairP4->Eta());
    // hTHINjetPairP4Rapidity->Fill(thinjetPairP4->Rapidity());
    
    nTotal++;
    
    
  }
  
  TFile* outFile = new TFile(outputFile.data(), toRecreateOutFile ? "recreate" : "update");
  hTHINNJetMatched->Write();
  // hTHINjetPairP4M->Write();
  // hTHINjetPairP4Pt->Write();
  // hTHINjetPairP4Rho->Write();
  // hTHINjetPairP4Phi->Write();
  // hTHINjetPairP4Theta->Write();
  // hTHINjetPairP4Eta->Write();
  // hTHINjetPairP4Rapidity->Write();
  outFile->Close();

}
