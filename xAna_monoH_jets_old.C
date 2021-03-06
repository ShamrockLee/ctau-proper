#include <TClonesArray.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <assert.h>

#include <algorithm>  // for std::next_permutation(v)
#include <fstream>
#include <iostream>

#include "untuplizer.h"

template <typename E>
size_t countUnique(std::vector<E>& vec) {
  if (vec.size() <= 1) {
    return vec.size();
  }
  std::vector<E> vecSorted(vec);
  size_t result = 1;
  for (typename std::vector<E>::iterator curr = vecSorted.begin(),
                                         next = vecSorted.begin();
       ++next != vecSorted.end(); curr++) {
    if (*curr != *next) {
      result++;
    }
  }
  return result;
}

void xAna_monoH_jets(std::string inputFile, std::string outputFile,
                     bool toRecreateOutFile = true, bool debug = false) {
  // bool isHerwigpp=(inputFile.find("herwigpp")!= std::string::npos);

  TreeReader data(inputFile.data());
  Long64_t nTotal = 0;
  Long64_t nPass[20] = {0};

  TH1F* hNumDWanted = new TH1F("hNumDWanted", "number of d or dbar from chi2",
                               8, 0 - 0.5f, 8 - 0.5f);
  TH1F* hIsDsFoundSignCorrect =
      new TH1F("hIsDsFoundSignCorrect",
               "Entrys with and without a d-dbar pair from chi2 and another "
               "pair from chi2-bar",
               2, 0 - 0.5f, 2 - 0.5f);
  TH1F* hnTHINMatchedSeparately =
      new TH1F("hnTHINMatchedSeparately",
               "Number of seperately matched thin jets (d or dbar from chi2)",
               5, 0 - 0.5f, 5 - 0.5f);
  TH1F* hnTHINMatched = new TH1F(
      "hnTHINMatched", "Number of matched thin jets (d or dbar from chi2)", 5,
      0 - 0.5f, 5 - 0.5f);
  TH1F* hnTHINMatchedUnique =
      new TH1F("hnTHINMatchedUnique",
               "Number of unique thin jets matched (d or dbar from chi2)", 5,
               0 - 0.5f, 5 - 0.5f);
  TH1F* hDeltaRTHINjetPairsFromChi2ordi =
      new TH1F("hDeltaRTHINjetPairsFromChi2ordi",
               "deltaR of d-pairs from (ordinary) chi2", 100, 0, 6);
  TH1F* hDeltaRTHINjetPairsFromChi2bar =
      (TH1F*)hDeltaRTHINjetPairsFromChi2ordi->Clone(
          "hDeltaRTHINjetPairsFromChi2bar");
  hDeltaRTHINjetPairsFromChi2bar->SetTitle("deltaR of d-pairs from chi2-bar");
  TH1F* hDeltaRBetweenTwoTHINjetPairs = new TH1F(
      "hDeltaRBetweenTwoTHINjetPairs", "deltaR between two d-pairs", 100, 0, 8);
  // TH1F* hTHINjetPairP4M = new TH1F("hTHINjetPairP4M", "Static Mass of THINjet
  // Pairs", 50, 0, 500);
  TH1F* hTHINjetP4Pt00 =
      new TH1F("hTHINjetP4Pt00", "Pt of the first (ordinary-ordinary) THINjet",
               100, 0, 1000);
  TH1F* hDP4Pt00Mismatched = new TH1F(
      "hDP4Pt00Mismatched",
      "Pt of the first (ordinary-ordinary) d-quark which isn't matched", 100, 0,
      500);
  TH1F* hDP4Pt01Mismatched =
      (TH1F*)hDP4Pt00Mismatched->Clone("hDP4Pt01Mismatched");
  hDP4Pt01Mismatched->SetTitle(
      "Pt of the second (ordinary-anti) d-quark which isn't matched");
  TH1F* hDP4Pt10Mismatched =
      (TH1F*)hDP4Pt00Mismatched->Clone("hDP4Pt10Mismatched");
  hDP4Pt10Mismatched->SetTitle(
      "Pt of the third (anti-ordinary) d-quark which isn't matched");
  TH1F* hDP4Pt11Mismatched =
      (TH1F*)hDP4Pt00Mismatched->Clone("hDP4Pt11Mismatched");
  hDP4Pt11Mismatched->SetTitle(
      "Pt of the second (anti-anti) d-quark which isn't matched");
  TH1F* arrHDP4PtMismatched[] = {hDP4Pt00Mismatched, hDP4Pt01Mismatched,
                                 hDP4Pt10Mismatched, hDP4Pt11Mismatched};
  TH1F* hTHINjetPairP4Pt0 =
      new TH1F("hTHINjetPairP4Pt0", "Pt of the first (ordinary) THINjet Pairs",
               100, 0, 1000);
  // TH1F* hTHINjetPairP4Rho = new TH1F("hTHINjetPairP4Rho", "Rho of THINjet
  // Pairs", 500, 0, 5000); TH1F* hTHINjetPairP4Phi = new
  // TH1F("hTHINjetPairP4Phi", "Phi of THINjet Pairs", 200, -M_PI, M_PI); TH1F*
  // hTHINjetPairP4Theta = new TH1F("hTHINjetPairP4Theta", "Theta of THINjet
  // Pairs", 100, 0, M_PI);
  TH1F* hTHINjetP4Eta00 =
      new TH1F("hTHINjetP4Eta00", "Pseudorapidity (Eta) of the first THINjet",
               100, -5, 5);
  TH1F* hTHINjetPairP4Eta0 =
      new TH1F("hTHINjetPairP4Eta0",
               "Pseudorapidity (Eta) of the first THINjet Pairs", 100, -5, 5);
  // TH1F* hTHINjetPairP4Rapidity =
  // (TH1F*)hTHINjetPairP4Eta->Clone("hTHINjetPairP4Rapidity");
  // hTHINjetPairP4Rapidity->SetTitle("Rapidity of THINjet Pairs");

  std::ofstream file_mismatched;
  file_mismatched.open("mismatch_high_Mx2-150_Mx1-1_ctau-1_20200723.txt");

  for (Long64_t jEntry = 0; jEntry < data.GetEntriesFast(); jEntry++) {
    data.GetEntry(jEntry);
    Int_t nGenPar = data.GetInt("nGenPar");
    Int_t* genParId = data.GetPtrInt("genParId");
    Int_t* genParSt = data.GetPtrInt("genParSt");
    Int_t* genMomParId = data.GetPtrInt("genMomParId");

    // int genHIndex[2]={-1,-1};
    std::vector<Int_t>
        vIndexesDWated;  // It's GOOD to declear variaable INSIDE the loop.
    vIndexesDWated.clear();

    for (int ig = 0; ig < nGenPar; ig++) {
      if (abs(genParId[ig]) != 1) continue;

      if (abs(genMomParId[ig]) != 18) continue;

      // if(genHIndex[0]<0)
      //     genHIndex[0]=ig;

      // else if(genHIndex[1]<0)
      //     genHIndex[1]=ig;

      vIndexesDWated.push_back(ig);
    }

    std::vector<Int_t> vIndexesDWatedSignCorrect;
    vIndexesDWatedSignCorrect.clear();
    bool isDsFoundSignCorrect = false;
    bool boolsSignCheckResult[4];
    for (int k = 0; k < 4; k++) boolsSignCheckResult[k] = false;
    // Int_t vIndexesDWatedOrdered[4];
    // if (vIndexesDWated.size() == 4) {
    //   for (int i = 0; i < 2; i++) {
    //     for (int j = 0; j < 2; j++) {
    //       for (int ig = 0; ig < nGenPar; ig++) {
    //         if ((genParId[ig] * (i == 1 ? -1 : 1) > 0) &&
    //             (genMomParId[ig] * (j == 1 ? -1 : 1) > 0)) {
    //           boolsSignCheckResult[i + (j << 1)] = true;
    //           vIndexesDWatedSignCorrect.push_back(ig);
    //           vIndexesDWatedOrdered[i + (j << 1)] = ig;
    //         }
    //       }
    //     }
    //   }
    //   isDsFoundSignCorrect = true;
    //   for (int k = 0; k < 4; k++) {
    //     if (!boolsSignCheckResult[k]) isDsFoundSignCorrect = false;
    //   }
    // }
    // hNumDWanted->Fill(vIndexesDWated.size()); //This is temporarily -- to
    // test if it is affected by the bow operations

    // if(genHIndex[0]<0 || genHIndex[1]<0)continue;
    // nPass[0]++;

    // if(genHIndex[0]==genHIndex[1])continue;
    // nPass[1]++;

    // Compare THINjetP4 (RECO) and genParP4 (GEN) to check if the direction
    // matches. TLorentzVector genH_l4[2];
    TClonesArray* genParP4 = (TClonesArray*)data.GetPtrTObject("genParP4");

    // for(int ih=0; ih<2; ih++)
    // genH_l4[ih] = *((TLorentzVector*)genParP4->At(genHIndex[ih]));

    if (debug) {
      // std::cout << genHIndex[0] << "\t" << genHIndex[1] << std::endl;
      // genH_l4[0].Print();
      // genH_l4[1].Print();
      std::cout << vIndexesDWated.size() << " jets: ";
      for (Int_t ig : vIndexesDWated) std::cout << ig << " ";
      std::cout << std::endl;

      for (Int_t ig : vIndexesDWated) genParP4->At(ig)->Print();

      std::cout << std::endl;
    }
    int nTHINJet = data.GetInt("THINnJet");
    const int nJets = nTHINJet;
    TClonesArray* thinjetP4 = (TClonesArray*)data.GetPtrTObject("THINjetP4");

    // check matching first

    // TLorentzVector *thinjetPairP4 = nullptr;

    // bool findAMatch=false;
    const float dRMax =
        0.4;  // for THINjet (0.4 rad) (whereas the angle of FATjet is 0.8 rad)
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
    if (isDsFoundSignCorrect) {
      vIndexesDWatedSwapped = vIndexesDWatedSignCorrect;
    } else {
      vIndexesDWatedSwapped = vIndexesDWated;
    }
    std::sort(vIndexesDWatedSwapped.begin(), vIndexesDWatedSwapped.end());

    // Check the signs of genParId's and genMomParId's
    Int_t nParticlesToMatch = 4;
    bool areAllMatchedSeparately = false;

    std::vector<Int_t> vIndexesJetMatched;
    vIndexesJetMatched.clear();
    bool findAMatch = false;
    Int_t nJetMatched = 0;
    Int_t nJetMatchedUnique = 0;

    if (isDsFoundSignCorrect ||
        ((Int_t)vIndexesDWated.size() == nParticlesToMatch)) {
      if (debug) std::cout << "Entered if statement." << std::endl;

      // Match using DeltaR

      // seperately
      // Swap jet's permutation.
      /*
      do {
        int indexJetMatchedLast = -1;
        // int indexJLast = -1;
        // Iterate particles
        if (debug) {
          std::cout << "First is looping!" << std::endl;
          std::cout << "{";
          if (nJets > 0) {
            std::cout << vIndexesDWatedSwapped[0];
          }
          for (int k=1; k<nJets; k++) {
            std::cout << ", " << vIndexesDWatedSwapped[k];
          }
          std::cout << "}" << std::endl;
        }
        for (int i=0; i<nParticlesToMatch; i++) {
          // Iterate jets
          bool isThisParticleMatched = false;
          for (int j=indexJetMatchedLast + 1; j<nJets; j++) {
            // See if matched
            if (((TLorentzVector*)genParP4->At(vIndexesDWatedSwapped[i]))
                ->DeltaR(*(TLorentzVector*)thinjetP4->At(j)) < dRMax) {
              vIndexesJetMatchedSeparately.push_back(j);
              isThisParticleMatched = true;
              // See if all matched (if it is the last one)
              indexJetMatchedLast = j;
              if (i == nParticlesToMatch-1) {
                areAllMatchedSeparately = true;
                nJetMatchedSeparately = nParticlesToMatch;
              }
              break;
            }
          }
          if (areAllMatchedSeparately) break;
          if (!isThisParticleMatched) break;
        }
        if (areAllMatchedSeparately)break;
        // if one of them does not mach, swap the permutation and find again.
      } while (std::next_permutation(vIndexesDWatedSwapped.begin(),
      vIndexesDWatedSwapped.end()));
      */
      // Iterate the particles
      // if (debug) {
      //   std::cout << "vIndexesDWatedOrdered:";
      //   for (int idebug = 0; idebug < 4; idebug++)
      //     std::cout << "\t" << vIndexesDWatedOrdered[idebug];
      //   std::cout << std::endl;
      // }
      for (int i = 0; i < nParticlesToMatch; i++) {
        // Iterate the jets
        Double_t DeltaRMin = 1000000000.;
        for (int j = 0; j < nJets; j++) {
          // if (debug) std::cout<< "Second is looping!" << std::endl;
          Double_t DeltaRCurrent =
              ((TLorentzVector*)genParP4->At(vIndexesDWatedOrdered[i]))
                  ->DeltaR(*(TLorentzVector*)thinjetP4->At(j));
          if (DeltaRCurrent < DeltaRMin) {
            DeltaRMin = DeltaRCurrent;
          };
          // if (((TLorentzVector*)genParP4->At(vIndexesDWated[i]))
          //         ->DeltaR(*(TLorentzVector*)thinjetP4->At(j)) < dRMax) {
          //   vIndexesJetMatched.push_back(j);
          //   break;
          // }
          if (DeltaRMin <= dRMax) {
            vIndexesJetMatched.push_back(j);
          } else {
            int ig = vIndexesDWatedOrdered[j];
            if (debug) std::cout << "ig: " << "vIndexesDWantedOrdered" << ig << std::endl;
            if (((TLorentzVector*)(genParP4->At(ig)))->Pt() >= 100) {
              file_mismatched
                  << jEntry << "\t" << ((j >= 2) ? "-" : "+")
                  << (((j & 1) == 1) ? "-" : "+") << " :\t"
                  << "Pt: " << ((TLorentzVector*)(genParP4->At(ig)))->Pt()
                  << "\t"
                  << "DeltaRMin: " << DeltaRMin << std::endl;
            }
          }
        }
      }
      nJetMatched = vIndexesJetMatched.size();
      nJetMatchedUnique = countUnique(vIndexesJetMatched);
      if (nJetMatched == nParticlesToMatch) {
        findAMatch = true;
      }

      // Check id signs
      std::vector<Int_t> vIndexesJetMatchedOrdered;
      vIndexesJetMatchedOrdered.resize(4);
      if (debug)
        for (int k = 0; k < 4; k++) vIndexesJetMatchedOrdered[k] = -1;
      bool boolsSignCheckResult[4];
      for (int k = 0; k < 4; k++) boolsSignCheckResult[k] = false;
      for (int k = 0; k < 4; k++) vIndexesJetMatchedOrdered[k] = -1;
      for (int i = 0; i < (Int_t)vIndexesJetMatched.size(); i++) {
        Int_t k = (genParId[vIndexesDWated[i]] < 0 ? 0b01 : 0) +
                  (genMomParId[vIndexesDWated[i]] < 0 ? 0b10 : 0);
        vIndexesJetMatchedOrdered[k] = vIndexesJetMatched[i];
        boolsSignCheckResult[k] = true;
      }
      /*
      bool boolsSignCheckResult[4];
      for (int k = 0; k < 4; k++) boolsSignCheckResult[k] = false;
      for (int k = 0; k < 4; k++) {
        int i = (genParId[vIndexesDWatedSwapped[k]] < 0) ? 1 : 0;
        int j = (genMomParId[vIndexesDWatedSwapped[k]] < 0) ? 1 : 0;
        boolsSignCheckResult[i + (j << 1)] = true;
      }
      */
      if (isDsFoundSignCorrect && findAMatch) {
        if (debug) {
          for (int k = 0; k < 4; k++) {
            if (!boolsSignCheckResult[k]) {
              printf("boolsSignCheckResult: {%d,\t%d,\t%d,\t%d}\n",
                     boolsSignCheckResult[0], boolsSignCheckResult[1],
                     boolsSignCheckResult[2], boolsSignCheckResult[3]);
              std::cout << "Makeup of vIndexesDWatedSwapped is not correct!"
                        << std::endl;
              // throw "Makeup of vIndexesDWatedSwapped is not correct!";
            }
          }
        }
        if (debug) std::cout << "Find a match!" << std::endl;

        if (debug) {
          printf("vIndexesJetMatchedOrdered: {%d,\t%d,\t%d,\t%d}\n",
                 vIndexesJetMatchedOrdered[0], vIndexesJetMatchedOrdered[1],
                 vIndexesJetMatchedOrdered[2], vIndexesJetMatchedOrdered[3]);
        }
        hDeltaRTHINjetPairsFromChi2ordi->Fill(
            ((TLorentzVector*)thinjetP4->At(vIndexesJetMatchedOrdered[0]))
                ->DeltaR(*(TLorentzVector*)thinjetP4->At(
                    vIndexesJetMatchedOrdered[1])));
        hDeltaRTHINjetPairsFromChi2bar->Fill(
            ((TLorentzVector*)thinjetP4->At(vIndexesJetMatchedOrdered[2]))
                ->DeltaR(*(TLorentzVector*)thinjetP4->At(
                    vIndexesJetMatchedOrdered[3])));
        hTHINjetP4Pt00->Fill(
            ((TLorentzVector*)thinjetP4->At(vIndexesJetMatchedOrdered[0]))
                ->Pt());
        hTHINjetP4Eta00->Fill(
            ((TLorentzVector*)thinjetP4->At(vIndexesJetMatchedOrdered[0]))
                ->Eta());
        TLorentzVector jetPairsP4[2];
        jetPairsP4[0] =
            (*(TLorentzVector*)thinjetP4->At(vIndexesJetMatchedOrdered[0]) +
             *(TLorentzVector*)thinjetP4->At(vIndexesJetMatchedOrdered[1]));
        jetPairsP4[1] =
            *(TLorentzVector*)thinjetP4->At(vIndexesJetMatchedOrdered[2]) +
            *(TLorentzVector*)thinjetP4->At(vIndexesJetMatchedOrdered[3]);
        hDeltaRBetweenTwoTHINjetPairs->Fill(
            jetPairsP4[0].DeltaR(jetPairsP4[1]));
        hTHINjetPairP4Pt0->Fill(jetPairsP4[0].Pt());
        hTHINjetPairP4Eta0->Fill(jetPairsP4[0].Eta());
      };

      // for (int k = 0; k < 4; k++) {
      //   if (!boolsSignCheckResult[k]) {
      //     for (int ig = 0; ig < nGenPar; ig++) {
      //       if (genParId[ig] == (((k & 1) == 1) ? -1 : 1) &&
      //           genMomParId[ig] == ((k >= 2) ? -18 : 18)) {
      //         TLorentzVector* parP4 = (TLorentzVector*)genParP4->At(ig);
      //         Double_t parPt = parP4->Pt();
      //         arrHDP4PtMismatched[k]->Fill(parPt);
      //         if (parPt >= 100) {
      //           file_mismatched << jEntry << "\t" << ((k >= 2) ? "-" : "+")
      //                           << (((k & 1) == 1) ? "-" : "+")
      //                           << " :\tPt: " << parPt;
      //           Double_t DeltaRMin = 1000000000.;
      //           for (int i = 0; i < nJets; i++) {
      //             Double_t DeltaRCurrent =
      //                 ((TLorentzVector*)genParP4->At(ig))
      //                     ->DeltaR(*(TLorentzVector*)thinjetP4->At(i));
      //             if (DeltaRCurrent < DeltaRMin) {
      //               DeltaRMin = DeltaRCurrent;
      //             };
      //           };
      //           file_mismatched << "\t"
      //                           << "DeltaRMin: " << DeltaRMin;
      //           file_mismatched << std::endl;
      //         }
      //         break;
      //       }
      //     }
      //   }
      // }
    }

    // if(!findAMatch)continue;

    hNumDWanted->Fill(vIndexesDWated.size());
    hIsDsFoundSignCorrect->Fill((int)isDsFoundSignCorrect);
    hnTHINMatchedSeparately->Fill(nJetMatchedSeparately);
    hnTHINMatched->Fill(nJetMatched);
    hnTHINMatchedUnique->Fill(nJetMatchedUnique);
    // hTHINjetPairP4M->Fill(thinjetPairP4->M());
    // hTHINjetPairP4Pt->Fill(thinjetPairP4->Pt());
    // hTHINjetPairP4Rho->Fill(thinjetPairP4->Rho());
    // hTHINjetPairP4Phi->Fill(thinjetPairP4->Phi());
    // hTHINjetPairP4Theta->Fill(thinjetPairP4->Theta());
    // hTHINjetPairP4Eta->Fill(thinjetPairP4->Eta());
    // hTHINjetPairP4Rapidity->Fill(thinjetPairP4->Rapidity());

    nTotal++;
  }

  file_mismatched.close();

  TFile* outFile =
      new TFile(outputFile.data(), toRecreateOutFile ? "recreate" : "update");
  hNumDWanted->Write();
  hIsDsFoundSignCorrect->Write();
  // hnTHINMatchedSeparately->Write(); // The first loop hasn't been fixed.
  hnTHINMatched->Write();
  hnTHINMatchedUnique->Write();
  hDeltaRTHINjetPairsFromChi2ordi->Write();
  hDeltaRTHINjetPairsFromChi2bar->Write();
  hDeltaRBetweenTwoTHINjetPairs->Write();
  // hTHINjetPairP4M->Write();
  hTHINjetP4Pt00->Write();
  hTHINjetPairP4Pt0->Write();
  // hTHINjetPairP4Rho->Write();
  // hTHINjetPairP4Phi->Write();
  // hTHINjetPairP4Theta->Write();
  hTHINjetP4Eta00->Write();
  // hDP4Pt00Mismatched->Write();
  for (int k = 0; k < 4; k++) arrHDP4PtMismatched[k]->Write();
  hTHINjetPairP4Eta0->Write();
  // hTHINjetPairP4Rapidity->Write();
  outFile->Close();
}
