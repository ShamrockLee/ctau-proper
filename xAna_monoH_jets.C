#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>
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

void xAna_monoH_jets(std::string inputFile,
                     std::string outputFileHead, // "jets"
                     std::string outputFileVar, // "Mx2-150_Mx1-1_ctau-1"
                     std::string outputFileTail, // "20200730"
                     bool toRecreateOutFile = true, bool debug = false) {
  TreeReader data(inputFile.data());
  Long64_t nTotal = 0;
  Long64_t nPass[20] = {0};

  TH1F* hNumDWanted = new TH1F("hNumDWanted", "number of d or dbar from chi2",
                               8, 0 - 0.5f, 8 - 0.5f);

  // TH1F* hIsDsFoundSignCorrect =
  //     new TH1F("hIsDsFoundSignCorrect",
  //              "Entrys with and without a d-dbar pair from chi2 and another "
  //              "pair from chi2-bar",
  //              2, 0 - 0.5f, 2 - 0.5f);
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
  TH1F *arrHDeltaRMinTHINjet[4], *arrHDeltaRDPairs[4];
  TH2F *arrHDeltaRMinTHINjetVsPt[4];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      arrHDeltaRMinTHINjet[(i << 1) + j] = new TH1F(
          (TString) "hDeltaRMinTHINjet" + i + j,
          (TString) "Minimal DeltaR between d" + i + j + " and THIN jets", 100,
          0, 2.f);
      arrHDeltaRMinTHINjetVsPt[(i << 1) + j] = new TH2F(
          (TString) "hDeltaRMinTHINjetVsPt" + i + j,
          (TString) "Minimal DeltaR between d" + i + j + " and THIN jets vs. the Pt of d" + i + j, 100,
          0, 2.f, 100, 0, 200);
    }
    arrHDeltaRDPairs[i] =
        new TH1F((TString) "hDeltaDPairs" + i,
                 (TString) "DeltaR of d-pairs from " +
                     (i == 1 ? "chi2-bar" : "(ordinary) chi2"),
                 100, 0, 6);
  }
  TH1F* hDeltaRBetweenTwoDPairs = new TH1F(
      "hDeltaRBetweenTwoDPairs", "DeltaR between two d-pairs", 100, 0, 8);
  TH1F* hDeltaRTHINjetPairs0 =
      new TH1F("hDeltaRTHINjetPairs0",
               "deltaR of d-jet-pairs from (ordinary) chi2", 100, 0, 6);
  TH1F* hDeltaRTHINjetPairs1 =
      (TH1F*)hDeltaRTHINjetPairs0->Clone("hDeltaRTHINjetPairs1");
  hDeltaRTHINjetPairs1->SetTitle("deltaR of d-jet-pairs from chi2-bar");
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
  file_mismatched.open((TString) "mismatch_" + outputFileVar + "_" +
                       outputFileTail + ".txt");

  for (Long64_t jEntry = 0; jEntry < data.GetEntriesFast(); jEntry++) {
    data.GetEntry(jEntry);

    Int_t nGenPar = data.GetInt("nGenPar");
    Int_t* genParId = data.GetPtrInt("genParId");
    Int_t* genParSt = data.GetPtrInt("genParSt");
    Int_t* genMomParId = data.GetPtrInt("genMomParId");

    // int genHIndex[2]={-1,-1};
    // It's GOOD to declear variaable INSIDE the loop.
    std::vector<Int_t> vIndexesDWatedOrdered(4);
    bool boolsDsFoundSignCorrect[4];
    for (int k = 0; k < 4; k++) boolsDsFoundSignCorrect[k] = false;
    for (int ig = 0; ig < nGenPar; ig++) {
      if (abs(genParId[ig]) != 1) continue;
      if (abs(genMomParId[ig]) != 18) continue;

      // if(genHIndex[0]<0)
      //     genHIndex[0]=ig;

      // else if(genHIndex[1]<0)
      //     genHIndex[1]=ig;

      int k = ((genMomParId[ig] < 0) << 1) + (genParId[ig] < 0);
      vIndexesDWatedOrdered[k] = ig;
      boolsDsFoundSignCorrect[k] = true;
    }
    bool isDsFoundSignCorrect = true;
    for (bool isDFound : boolsDsFoundSignCorrect)
      if (!isDFound) isDsFoundSignCorrect = false;

    TClonesArray* genParP4 = (TClonesArray*)data.GetPtrTObject("genParP4");
    TClonesArray* thinjetP4 = (TClonesArray*)data.GetPtrTObject("THINjetP4");
    Int_t nTHINJet = data.GetInt("THINnJet");

    // for THINjet (0.4 rad) (whereas the angle of FATjet is 0.8 rad)
    float dRMaxTHIN = 0.4;

    if (isDsFoundSignCorrect) {
      TLorentzVector* p4DWanted[4];
      for (int k = 0; k < 4; k++)
        p4DWanted[k] = (TLorentzVector*)genParP4->At(vIndexesDWatedOrdered[k]);

      for (int i=0; i<2; i++) {
        arrHDeltaRDPairs[i]->Fill(p4DWanted[(i<<1)]->DeltaR(*p4DWanted[(i<<1) + 1]));
      }
      hDeltaRBetweenTwoDPairs->Fill((*p4DWanted[0] + *p4DWanted[1]).DeltaR(*p4DWanted[2] + *p4DWanted[3]));
      std::vector<Int_t> vIndexesTHINjetMatchedUnfiltered;
      vIndexesTHINjetMatchedUnfiltered.clear();
      std::vector<Double_t> vDeltaRTHINjetMatchedUnfiltered;
      vDeltaRTHINjetMatchedUnfiltered.clear();
      for (int i = 0; i < 4; i++) {
        Double_t dRThisMin = 10000.;
        Int_t indexTHINjetMatched;
        for (int j = 0; j < nTHINJet; j++) {
          Double_t dRThisCurrent =
              p4DWanted[i]->DeltaR(*(TLorentzVector*)thinjetP4->At(j));
          if (dRThisCurrent < dRThisMin) {
            dRThisMin = dRThisCurrent;
            indexTHINjetMatched = j;
          };
        }
        vIndexesTHINjetMatchedUnfiltered.push_back(indexTHINjetMatched);
        vDeltaRTHINjetMatchedUnfiltered.push_back(dRThisMin);
        arrHDeltaRMinTHINjet[i]->Fill(dRThisMin);
        arrHDeltaRMinTHINjetVsPt[i]->Fill(dRThisMin, p4DWanted[i]->Pt());
        if (dRThisMin > dRMaxTHIN) {
          arrHDP4PtMismatched[i]->Fill(p4DWanted[i]->Pt());
          file_mismatched << jEntry << "\t" << ((i >= 2) ? "-" : "+")
                          << (((i & 1) == 1) ? "-" : "+") << " :\t"
                          << "parPt: " << p4DWanted[i]->Pt() << "\t"
                          << "dRThisMin: " << dRThisMin << std::endl;
        }
      }
    }
  }

  file_mismatched.close();
  TString outputFile = (TString) "output_" + outputFileHead + "_" +
                       outputFileVar + "_" + outputFileTail + ".root";
  TFile* outFile =
      new TFile(outputFile.Data(), toRecreateOutFile ? "recreate" : "update");

  for (int i=0; i<2; i++) {
    arrHDeltaRDPairs[i]->Write();
  }
  hDeltaRBetweenTwoDPairs->Write();
  for (int k=0; k<4; k++) {
    arrHDeltaRMinTHINjet[k]->Write();
    arrHDeltaRMinTHINjetVsPt[k]->Write();
    arrHDP4PtMismatched[k]->Write();
  }
}