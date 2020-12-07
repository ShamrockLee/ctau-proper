#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <assert.h>

#include <TCanvas.h>
#include <TStyle.h>

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

void xAna_monoZ_jets(std::string inputFile,
                     TString outputFile,
                    //  std::string outputFileHead,  // "jets-newdata"
                    //  std::string outputFileVar,   // "Mx2-150_Mx1-1_ctau-1"
                    //  std::string outputFileTail,  // "20200730"
                     Bool_t toRecreateOutFile = true, Bool_t debug = false) {
  TreeReader data(inputFile.data(), "Events");
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
  TH1F *arrHDeltaRMinTHINjet[4], *arrHDeltaRDPairs[4], *arrHDP4Pt[4];
  TH2F* arrHDeltaRMinTHINjetVsPt[4];
  std::vector<Int_t> vDeltaRMaxX100;
  vDeltaRMaxX100.clear();
  const Int_t nDeltaRMax = 6;
  TClonesArray tcflHTHINjetP4Pt("TH1F", (nDeltaRMax<<2)), tcflHTHINjetP4PtMismatched("TH1F", (nDeltaRMax << 2));
  for (Int_t val = 40; val <= 50; val += 2) vDeltaRMaxX100.push_back(val);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      int k = (i << 1) + j;
      arrHDeltaRMinTHINjet[k] = new TH1F(
          (TString) "hDeltaRMinTHINjet" + i + j,
          (TString) "Minimal DeltaR between d" + i + j + " and THIN jets", 100,
          0, 2.f);
      arrHDeltaRMinTHINjetVsPt[k] =
          new TH2F((TString) "hDeltaRMinTHINjetVsPt" + i + j,
                   (TString) "Minimal DeltaR between d" + i + j +
                       " and THIN jets vs. the Pt of d" + i + j,
                   100, 0, 2.f, 100, 0, 200);
      arrHDP4Pt[k] = new TH1F((TString) "hDP4Pt" + i + j,
                              (TString) "Pt of d" + i + j, 100, 0, 1000);
    }
    arrHDeltaRDPairs[i] =
        new TH1F((TString) "hDeltaRDPairs" + i,
                 (TString) "DeltaR of d-pairs from " +
                     (i == 1 ? "chi2-bar" : "(ordinary) chi2"),
                 100, 0, 6);
  }
  for (uint m=0; m<nDeltaRMax; m++) {
    Int_t deltaRMaxX100 = vDeltaRMaxX100[m];
    //TH1F *arrHTHINjetP4Pt[4], *arrHTHINjetP4PtMismatched[4];
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        int k = (i << 1) + j;
        TString nameHead = (TString) "hTHINjetP4Pt" + i + j;
        TString titleHead = (TString) "The Pt of d" + i + j;
        Int_t indexTcfl = (m<<2) + k;
        if (debug) std::cout << "testing: " << "k: " << k << std::endl;
        new(tcflHTHINjetP4Pt[indexTcfl]) TH1F(nameHead + "Matched" + deltaRMaxX100,
                     titleHead + "that matches (dRMax = " + deltaRMaxX100 + ")",
                     100, 0, 500);
        if (debug) std::cout << "testing " << "generating the next one";
        // tcflHTHINjetP4PtMismatched[indexTcfl] = (TH1F*)(tcflHTHINjetP4Pt[indexTcfl]->Clone(
        //     nameHead + "Mismatched" + deltaRMaxX100));
        // ((TH1F*)tcflHTHINjetP4PtMismatched[indexTcfl])->SetTitle(
        //     titleHead + "that doesn't match (dRMax = " + deltaRMaxX100 + ")");
        new(tcflHTHINjetP4PtMismatched[indexTcfl]) TH1F(nameHead + "Mismatched" + deltaRMaxX100,
                     titleHead + "that doesn't match (dRMax = " + deltaRMaxX100/100 + ")",
                     100, 0, 500);
      }
    }
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

  // TH1F* hTHINjetP4Pt00 =
  //     new TH1F("hTHINjetP4Pt00", "Pt of the first (ordinary-ordinary)
  //     THINjet",
  //              100, 0, 1000);
  // TH1F* hDP4Pt00Mismatched = new TH1F(
  //     "hDP4Pt00Mismatched",
  //     "Pt of the first (ordinary-ordinary) d-quark which isn't matched", 100,
  //     0, 500);
  // TH1F* hDP4Pt01Mismatched =
  //     (TH1F*)hDP4Pt00Mismatched->Clone("hDP4Pt01Mismatched");
  // hDP4Pt01Mismatched->SetTitle(
  //     "Pt of the second (ordinary-anti) d-quark which isn't matched");
  // TH1F* hDP4Pt10Mismatched =
  //     (TH1F*)hDP4Pt00Mismatched->Clone("hDP4Pt10Mismatched");
  // hDP4Pt10Mismatched->SetTitle(
  //     "Pt of the third (anti-ordinary) d-quark which isn't matched");
  // TH1F* hDP4Pt11Mismatched =
  //     (TH1F*)hDP4Pt00Mismatched->Clone("hDP4Pt11Mismatched");
  // hDP4Pt11Mismatched->SetTitle(
  //     "Pt of the second (anti-anti) d-quark which isn't matched");
  // TH1F* arrHDP4PtMismatched[] = {hDP4Pt00Mismatched, hDP4Pt01Mismatched,
  //                                hDP4Pt10Mismatched, hDP4Pt11Mismatched};
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

  // std::ofstream file_mismatched;
  // file_mismatched.open((TString) "mismatch_" + outputFileVar + "_" +
  //                      outputFileTail + ".txt");

  for (Long64_t jEntry = 0; jEntry < data.GetEntriesFast(); jEntry++) {
    data.GetEntry(jEntry);

    Int_t nGenPar = data.GetInt("nGenPart");
    Int_t* genParId = data.GetPtrInt("GenPart_pdgId");
    // Int_t* genParSt = data.GetPtrInt("genParSt");
    Int_t* genMomParIndex = data.GetPtrInt("GenPart_genPartIdxMother");

    // int genHIndex[2]={-1,-1};
    // It's GOOD to declear variable INSIDE the loop.
    std::vector<Int_t> vIndexesDWatedOrdered(4);
    bool boolsDsFoundSignCorrect[4];
    for (int k = 0; k < 4; k++) boolsDsFoundSignCorrect[k] = false;
    for (int ig = 0; ig < nGenPar; ig++) {
      if (abs(genParId[ig]) != 1) continue;
      if (abs(genParId[genMomParIndex[ig]]) != 18) continue;

      // if(genHIndex[0]<0)
      //     genHIndex[0]=ig;

      // else if(genHIndex[1]<0)
      //     genHIndex[1]=ig;

      int k = ((genParId[genMomParIndex[ig]] < 0) << 1) + (genParId[ig] < 0);
      vIndexesDWatedOrdered[k] = ig;
      boolsDsFoundSignCorrect[k] = true;
    }
    bool isDsFoundSignCorrect = true;
    for (bool isDFound : boolsDsFoundSignCorrect)
      if (!isDFound) isDsFoundSignCorrect = false;

    // TClonesArray* genParP4 = (TClonesArray*)data.GetPtrTObject("genParP4");
    // TClonesArray* thinjetP4 = (TClonesArray*)data.GetPtrTObject("THINjetP4");
    Float_t *genParPt = data.GetPtrFloat("GenPart_pt");
    Float_t *genPartEta = data.GetPtrFloat("GenPart_eta");
    Float_t *genPartPhi = data.GetPtrFloat("GenPart_phi");
    Float_t *thinjetPt = data.GetPtrFloat("Jet_pt");
    Float_t *thinjetEta = data.GetPtrFloat("Jet_eta");
    Float_t *thinjetPhi = data.GetPtrFloat("Jet_phi");
    Int_t nTHINJet = data.GetInt("nJet");

    // for THINjet (0.4 rad) (whereas the angle of FATjet is 0.8 rad)
    float dRMaxTHIN = 0.4;

    if (isDsFoundSignCorrect) {
      // TLorentzVector* p4DWanted[4];
      // for (int k = 0; k < 4; k++)
      //   p4DWanted[k] = (TLorentzVector*)genParP4->At(vIndexesDWatedOrdered[k]);

      for (int i = 0; i < 2; i++) {
        arrHDeltaRDPairs[i]->Fill(
            // p4DWanted[(i << 1)]->DeltaR(*p4DWanted[(i << 1) + 1])
            TMath::Sqrt(
              (genPartEta[vIndexesDWatedOrdered[i<<1]]-genPartEta[vIndexesDWatedOrdered[(i<<1)+1]]) * (genPartEta[vIndexesDWatedOrdered[i<<1]]-genPartEta[vIndexesDWatedOrdered[(i<<1)+1]])
              + (genPartPhi[vIndexesDWatedOrdered[i<<1]]-genPartPhi[vIndexesDWatedOrdered[(i<<1)+1]]) * (genPartPhi[vIndexesDWatedOrdered[i<<1]]-genPartPhi[vIndexesDWatedOrdered[(i<<1)+1]]))
        );
      }
      // hDeltaRBetweenTwoDPairs->Fill((*p4DWanted[0] + *p4DWanted[1])
      //                                   .DeltaR(*p4DWanted[2] + *p4DWanted[3]));
      std::vector<Int_t> vIndexesTHINjetMatchedUnfiltered;
      vIndexesTHINjetMatchedUnfiltered.clear();
      std::vector<Double_t> vDeltaRTHINjetMatchedUnfiltered;
      vDeltaRTHINjetMatchedUnfiltered.clear();
      for (int i = 0; i < 4; i++) {
        Double_t dRThisMin = 10000.;
        Int_t indexTHINjetMatched;
        for (int j = 0; j < nTHINJet; j++) {
          Double_t dRThisCurrent =
              // p4DWanted[i]->DeltaR(*(TLorentzVector*)thinjetP4->At(j));
              TMath::Sqrt(
                (genPartEta[vIndexesDWatedOrdered[i]]-thinjetEta[j]) * (genPartEta[vIndexesDWatedOrdered[i]]-thinjetEta[j])
                +(genPartPhi[vIndexesDWatedOrdered[i]]-thinjetPhi[j]) * (genPartPhi[vIndexesDWatedOrdered[i]]-thinjetPhi[j])
              );
          if (dRThisCurrent < dRThisMin) {
            dRThisMin = dRThisCurrent;
            indexTHINjetMatched = j;
          };
        }
        vIndexesTHINjetMatchedUnfiltered.push_back(indexTHINjetMatched);
        vDeltaRTHINjetMatchedUnfiltered.push_back(dRThisMin);
        arrHDeltaRMinTHINjet[i]->Fill(dRThisMin);
        // arrHDeltaRMinTHINjetVsPt[i]->Fill(dRThisMin, p4DWanted[i]->Pt());
        arrHDeltaRMinTHINjetVsPt[i]->Fill(dRThisMin, genParPt[vIndexesDWatedOrdered[i]]);
        for (uint m = 0; m < nDeltaRMax; m++) {
          Double_t ptJetCurrent =
              // ((TLorentzVector*)thinjetP4->At(indexTHINjetMatched))->Pt();
              thinjetPt[indexTHINjetMatched];
          if (debug) std::cout << "testing " << "m: " << m << "i: " << i << "dRThisMin: " << dRThisMin << "vDeltaRMaxX100[i]: " << vDeltaRMaxX100[i] <<  std::endl;
          Int_t indexTcfl = (m<<2) + i;
          if (dRThisMin <= (Double_t)vDeltaRMaxX100[m] / 100) {
            if (debug) std::cout << "testing " << "tcflHTHINjetP4Pt[indexTcfc]: " << tcflHTHINjetP4Pt[indexTcfl] << std::endl;
            ((TH1F*)tcflHTHINjetP4Pt[indexTcfl])->Fill(ptJetCurrent);
          } else {
            ((TH1F*)tcflHTHINjetP4PtMismatched[indexTcfl])->Fill(ptJetCurrent);
          }
        }
        // if (dRThisMin > dRMaxTHIN) {
        //   file_mismatched << jEntry << "\t" << ((i >= 2) ? "-" : "+")
        //                   << (((i & 1) == 1) ? "-" : "+") << " :\t"
        //                   // << "parPt: " << p4DWanted[i]->Pt() << "\t"
        //                   << "parPt: " << genParPt[vIndexesDWatedOrdered[i]] << "\t"
        //                   << "dRThisMin: " << dRThisMin << std::endl;
        // }
      }
    }
  }
  // file_mismatched.close();

  TH1F *arrHDeltaRMinTHINjetCumuBack[4];
  for (int k=0; k<4; k++) {
    TH1F *hDeltaRMinTHINjetCurrent = arrHDeltaRMinTHINjet[k];
    TH1F *hDeltaRMinTHINjetCumuBackCurrent =
      (TH1F*)hDeltaRMinTHINjetCurrent->GetCumulative(kFALSE, "CumuBack");
    hDeltaRMinTHINjetCumuBackCurrent->SetTitle(
      (TString)hDeltaRMinTHINjetCurrent->GetTitle() + " (Cumulative)");
    arrHDeltaRMinTHINjetCumuBack[k] = hDeltaRMinTHINjetCumuBackCurrent;
  }

  // TString outputFile = (TString) "output_" + outputFileHead + "_" +
  //                      outputFileVar + "_" + outputFileTail + ".root";
  TFile* outFile = TFile::Open(outputFile.Data(), toRecreateOutFile ? "recreate" : "update");

  for (int i = 0; i < 2; i++) {
    arrHDeltaRDPairs[i]->Write();
  }
  hDeltaRBetweenTwoDPairs->Write();
  for (int k = 0; k < 4; k++) {
    arrHDeltaRMinTHINjet[k]->Write();
    arrHDeltaRMinTHINjet[k]->GetCumulative(kFALSE, "CumuBack")->Write();
    arrHDeltaRMinTHINjetVsPt[k]->Write();
    for (int m=0; m<nDeltaRMax; m++) {
      Int_t indexTcfl = (m<<2) + k;
      ((TH1F*)tcflHTHINjetP4Pt[indexTcfl])->Write();
      ((TH1F*)tcflHTHINjetP4PtMismatched[indexTcfl])->Write();
    }
  }
  if (false) {
    TCanvas *c1 = new TCanvas;
    gStyle->SetOptStat(111111);
    TString outputFileHead = "jets";
    TString outputFileVar = "testvar";
    TString outputFileTail = "20201209";
    TString outDir = (TString)"../out_images/output_" + outputFileHead + "_" + outputFileTail;
    TString outImageNameHead = (TString)outputFileHead + "_" + outputFileVar;
    TString outImageCommonPath = outDir + "/" + outImageNameHead + "_";
    for (int k=0; k<4; k++) {
      TH1F *hDeltaRMinTHINjetCurrent = arrHDeltaRMinTHINjet[k];
      hDeltaRMinTHINjetCurrent->Draw();
      c1->Print(outImageCommonPath + hDeltaRMinTHINjetCurrent->GetName() + ".svg");
      hDeltaRMinTHINjetCurrent->GetCumulative(kFALSE, "CumuBack")->Draw();
      c1->Print(outImageCommonPath + hDeltaRMinTHINjetCurrent->GetName() + "CumuBack" + ".svg");
      TH2F *hDeltaRMinTHINjetVsPtCurrent = arrHDeltaRMinTHINjetVsPt[k];
      c1->Print(outImageCommonPath + hDeltaRMinTHINjetVsPtCurrent->GetName() + ".svg");
      c1->Clear();
      ((TH1F*)tcflHTHINjetP4Pt[k])->Draw();
      for (int m=1; m<nDeltaRMax; m++) {
        TH1F *hTHINjetP4PtCurrent = (TH1F*)tcflHTHINjetP4Pt[(m<<2)+k];
        hTHINjetP4PtCurrent->Draw("SAME");
      }
      c1->Print(outImageCommonPath + "hTHINjetP4Pt" + (k>=2) + ((k & 1) > 0) + "40to50" + ".svg");
      c1->Clear();
      ((TH1F*)tcflHTHINjetP4PtMismatched[k])->Draw();
      for (int m=1; m<nDeltaRMax; m++) {
        TH1F *hTHINjetP4PtMismatchedCurrent = (TH1F*)tcflHTHINjetP4PtMismatched[(m<<2)+k];
        hTHINjetP4PtMismatchedCurrent->Draw("SAME");
      }
      c1->Print(outImageCommonPath + "hTHINjetP4Pt" + (k>=2) + ((k & 1) > 0) + "Mismatched40to50" + ".svg");
    }
    c1->Close();
  }
}
