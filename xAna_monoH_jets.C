#include <iostream>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include "untuplizer.h"
#include <algorithm> // for std::next_permutation(v)
#include <assert.h>

void xAna_monoH_jets(std::string inputFile,std::string outputFile, bool toRecreateOutFile=true, bool debug=false) {
  // bool isHerwigpp=(inputFile.find("herwigpp")!= std::string::npos);
  
  TreeReader data(inputFile.data());
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  
  TH1F* hNumDWanted = new TH1F("hNumDWanted", "number of d or dbar from chi2", 8, 0-0.5f, 8-0.5f); 
  TH1F* hIsDsFoundSignCorrect = new TH1F("hIsDsFoundSignCorrect", "Entrys with and without a d-dbar pair from chi2 and another pair from chi2-bar", 2, 0-0.5f, 2-0.5f);
  TH1F* hnTHINMatchedSeparately = new TH1F("hnTHINMatchedSeparately", "Number of seperately matched thin jets (d or dbar from chi2)", 5, 0-0.5f, 5-0.5f);
  TH1F* hnTHINMatched = new TH1F("hnTHINMatched", "Number of matched thin jets (d or dbar from chi2)", 5, 0-0.5f, 5-0.5f);
  TH1F* hDeltaRTHINjetPairsFromChi2ordi = new TH1F("hDeltaRTHINjetPairsFromChi2ordi", "deltaR of d-pairs from (ordinary) chi2", 100, 0, M_PI);
  TH1F* hDeltaRTHINjetPairsFromChi2bar = (TH1F*)hDeltaRTHINjetPairsFromChi2ordi->Clone("hDeltaRTHINjetPairsFromChi2bar");
  hDeltaRTHINjetPairsFromChi2bar->SetTitle("deltaR of d-pairs from chi2-bar");
  TH1F* hDeltaRBetweenTwoTHINjetPairs = (TH1F*)hDeltaRTHINjetPairsFromChi2ordi->Clone("hDeltaRBetweenTwoTHINjetPairs");
  hDeltaRBetweenTwoTHINjetPairs->SetTitle("deltaR between two d-pairs");
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
    std::vector<Int_t> vIndexesDWated; // It's GOOD to declear variaable INSIDE the loop.
    vIndexesDWated.clear();

    for(int ig=0; ig < nGenPar; ig++){

      if(abs(genParId[ig])!=1)continue;

      if(abs(genMomParId[ig])!=18)continue;

      // if(genHIndex[0]<0)
      //     genHIndex[0]=ig;


      // else if(genHIndex[1]<0)
      //     genHIndex[1]=ig;
      
      vIndexesDWated.push_back(ig);
      
    }
    
    std::vector<Int_t> vIndexesDWatedSignCorrectSeperated;
    vIndexesDWatedSignCorrectSeperated.clear();
    bool isDsFoundSignCorrect = true;
    {
      bool boolsSignCheckResult[4];
      for (int k=0; k<4; k++) boolsSignCheckResult[k] = false;
      for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++) {
          for (int ig=0; ig<nGenPar; ig++) {
            if ((genParId[ig] * (i==1 ? -1 : 1) > 0) &&
                (genMomParId[ig] * (j==1 ? -1 : 1) > 0)) {
              boolsSignCheckResult[i + (j<<1)] = true;
              vIndexesDWatedSignCorrectSeperated.push_back(ig);
            }
          }
        }
      }
      for (int k=0; k<4; k++) {
        if (!boolsSignCheckResult[k]) isDsFoundSignCorrect = false;
      }
    }
    // hNumDWanted->Fill(vIndexesDWated.size()); //This is temporarily -- to test if it is affected by the bow operations


    // if(genHIndex[0]<0 || genHIndex[1]<0)continue;
    // nPass[0]++;


    // if(genHIndex[0]==genHIndex[1])continue;
    // nPass[1]++;


    // Compare THINjetP4 (RECO) and genParP4 (GEN) to check if the direction matches.
    // TLorentzVector genH_l4[2];
    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");


    // for(int ih=0; ih<2; ih++)
    // genH_l4[ih] = *((TLorentzVector*)genParP4->At(genHIndex[ih]));




    if(debug){
    // std::cout << genHIndex[0] << "\t" << genHIndex[1] << std::endl;
    // genH_l4[0].Print();
    // genH_l4[1].Print();
    std::cout << vIndexesDWated.size() << " jets: ";
    for (Int_t ig: vIndexesDWated) std::cout << ig << " ";
    std::cout << std::endl;
    
    for (Int_t ig: vIndexesDWated) genParP4->At(ig)->Print();
    
    std::cout << std::endl;

    }
    int nTHINJet         = data.GetInt("THINnJet");
    const int nJets=nTHINJet;
    TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");


    // check matching first


    // TLorentzVector *thinjetPairP4 = nullptr;
    
    // bool findAMatch=false;
    const float dRMax=0.4; // for THINjet (0.4 rad) (whereas the angle of FATjet is 0.8 rad)
    // int matchedHJetIndex[2]={-1,-1};
            
    // for(int ij=0; ij<nJets; ij++)
    // {
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


    // }
    
    std::vector<Int_t> vIndexesJetMatchedSeparately;
    vIndexesJetMatchedSeparately.clear();
    Int_t nJetMatchedSeparately = 0;
    std::vector<Int_t> vIndexesDWatedSwapped;
    if (isDsFoundSignCorrect)
      vIndexesDWatedSwapped = vIndexesDWatedSignCorrectSeperated;
    else
      vIndexesDWatedSwapped = vIndexesDWated;
    std::sort(vIndexesDWatedSwapped.begin(), vIndexesDWatedSwapped.end());
    
    // Check the signs of genParId's and genMomParId's
    Int_t nParticlesToMatch = 4;
    bool findAMatchSeparately = false;
    if (isDsFoundSignCorrect || (Int_t)vIndexesDWated.size() == nParticlesToMatch) {
      
      // Match uning DeltaR
      do {
        int indexJetMatchedLast = -1;
        for (int i=0; i<nParticlesToMatch; i++) {
          for (int j=indexJetMatchedLast + 1; j<nJets; j++) {
            if (((TLorentzVector*)genParP4->At(vIndexesDWatedSwapped[i]))
                ->DeltaR(*(TLorentzVector*)thinjetP4->At(j)) < dRMax) {
               vIndexesJetMatchedSeparately.push_back(j);
              if (i == nParticlesToMatch-1) {
                findAMatchSeparately = true;
                nJetMatchedSeparately = nParticlesToMatch;
                break;
              }
              indexJetMatchedLast = j;
            }
          }
          if (findAMatchSeparately) break;
          if (indexJetMatchedLast==nJets-1) {
            if (i>nJetMatchedSeparately) {
              nJetMatchedSeparately = i;
              break;
            }
          }
        }
        if (findAMatchSeparately)break;
      } while (std::next_permutation(vIndexesDWatedSwapped.begin(), vIndexesDWatedSwapped.end()));
    }
    
    std::vector<Int_t> vIndexesJetMatched;
    vIndexesJetMatchedSeparately.clear();
    Int_t nJetMatched = 0;
    for (int i=0; i<nParticlesToMatch; i++) {
      for (int j=0; j<nJets; j++) {
        if (((TLorentzVector*)genParP4->At(vIndexesDWatedSwapped[i]))
              ->DeltaR(*(TLorentzVector*)thinjetP4->At(j)) < dRMax) {
              vIndexesJetMatched.push_back(j);
        }
      }
    }
    nJetMatched = vIndexesJetMatched.size();
    
    if (isDsFoundSignCorrect && findAMatchSeparately) {
      // Check id signs
      if (debug) {
        bool boolsSignCheckResult[4];
        for (int k=0; k<4; k++) boolsSignCheckResult[k] = false;
        for (int k=0; k<4; k++) {
          int i = (genParId[vIndexesDWatedSwapped[k]] < 0) ? 1 : 0;
          int j = (genMomParId[vIndexesDWatedSwapped[k]] < 0) ? 1 : 0;
          boolsSignCheckResult[i + (j<<1)] = true;
        }
        for (int k=0; k<4; k++) 
          if (!boolsSignCheckResult[k]) {
            printf("boolsSignCheckResult: {%d,\t%d,\t%d,\t%d}\n", 
                 boolsSignCheckResult[0], 
                 boolsSignCheckResult[1], 
                 boolsSignCheckResult[2], 
                 boolsSignCheckResult[3]);
            std::cout << "Makeup of vIndexesDWatedSwapped is not correct!" << std::endl;
            throw "Makeup of vIndexesDWatedSwapped is not correct!";
          }
      }
      if (debug) std::cout<< "Find a match!" << std::endl;
      std::vector<Int_t> vIndexesDPairsOrdered;
      vIndexesDPairsOrdered.resize(4);
      for (int i=0; i<4; i++) {
        vIndexesDPairsOrdered[
        (genParId[vIndexesDWatedSwapped[i]]<0 ? 0b01 : 0)
        + (genMomParId[vIndexesDWatedSwapped[i]]<0 ? 0b10 : 0)] =
        vIndexesJetMatchedSeparately[i];
      }
      if (debug) {
          printf("vIndexesDPairsOrdered: {%d,\t%d,\t%d,\t%d}\n", 
                 vIndexesDPairsOrdered[0], 
                 vIndexesDPairsOrdered[1], 
                 vIndexesDPairsOrdered[2], 
                 vIndexesDPairsOrdered[3]);
      }
      hDeltaRTHINjetPairsFromChi2ordi->Fill(
          ((TLorentzVector*)thinjetP4->At(vIndexesDPairsOrdered[0]))
          ->DeltaR(*(TLorentzVector*)thinjetP4->At(vIndexesDPairsOrdered[1]))
                                 );
      hDeltaRTHINjetPairsFromChi2bar->Fill(
          ((TLorentzVector*)thinjetP4->At(vIndexesDPairsOrdered[2]))
          ->DeltaR(*(TLorentzVector*)thinjetP4->At(vIndexesDPairsOrdered[3]))
                                 );
      hDeltaRBetweenTwoTHINjetPairs->Fill(
          (*(TLorentzVector*)thinjetP4->At(vIndexesDPairsOrdered[0]) 
          + *(TLorentzVector*)thinjetP4->At(vIndexesDPairsOrdered[1]))
          .DeltaR(*(TLorentzVector*)thinjetP4->At(vIndexesDPairsOrdered[2])
              + *(TLorentzVector*)thinjetP4->At(vIndexesDPairsOrdered[3]))
                                 );
    }

    // if(!findAMatch)continue;
    
    hNumDWanted->Fill(vIndexesDWated.size());
    hIsDsFoundSignCorrect->Fill((int)isDsFoundSignCorrect);
    hnTHINMatchedSeparately->Fill(nJetMatchedSeparately);
    hnTHINMatched->Fill(nJetMatched);
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
  hNumDWanted->Write();
  hIsDsFoundSignCorrect->Write();
  hnTHINMatchedSeparately->Write();
  hnTHINMatched->Write();
  hDeltaRTHINjetPairsFromChi2ordi->Write();
  hDeltaRTHINjetPairsFromChi2bar->Write();
  hDeltaRBetweenTwoTHINjetPairs->Write();
  // hTHINjetPairP4M->Write();
  // hTHINjetPairP4Pt->Write();
  // hTHINjetPairP4Rho->Write();
  // hTHINjetPairP4Phi->Write();
  // hTHINjetPairP4Theta->Write();
  // hTHINjetPairP4Eta->Write();
  // hTHINjetPairP4Rapidity->Write();
  outFile->Close();

}
