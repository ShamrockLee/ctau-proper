#include <Rtypes.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#include "HistMerger.C"

void xAna_monoZ_genhist(const TString nameCondorPack,
                        const TString nameDatagroup,
                        const TString nameClusterID,
                        TString dirCondorPacks = ".",
                        const Bool_t toRecreateOutFile = true,
                        const Bool_t debug = false,
                        const Bool_t allowMissing = false) {
  const TString seperatorPath = "/";
  const TString namesLepton[] = {"Electron", "Muon",
                                 "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon",
                                      "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};

  HistMerger* merger = new HistMerger;

  std::vector<TString> vNameTT;
  vNameTT.clear();
  UInt_t iTreeZMassCuttedElectron = 0;
  UInt_t iTreePreselectedMatchingElectron = 0;
  UInt_t iTreeAllMatchedElectron = iTreeZMassCuttedElectron = vNameTT.size();
  for (Byte_t i = 0; i < 2; i++) {
    TString nameLeptonCurrent = namesLepton[i];
    vNameTT.push_back((TString) "ZMassCutted" + nameLeptonCurrent);
  }
  if (nameDatagroup.SubString("signal").Start() == 0) {
    for (Byte_t i = 0; i < 2; i++) {
      vNameTT.push_back("Gen" + namesLepton[i]);
    }
    iTreePreselectedMatchingElectron = vNameTT.size();
    for (Byte_t i = 0; i < 2; i++) {
      vNameTT.push_back("PreselectedMatching" + namesLepton[i]);
    }
    for (Byte_t i = 0; i < 2; i++) {
      vNameTT.push_back("AllMatched" + namesLepton[i]);
    }
  }

  std::vector<UInt_t> vNumberFile;
  vNumberFile.clear();
  std::vector<Double_t> vCrossSection;
  vCrossSection.clear();

  TString dirCondorPackCurrent =
      dirCondorPacks + seperatorPath + nameCondorPack;
  TString pathTextFileNumbers = dirCondorPackCurrent + seperatorPath +
                                "fileNumbers_" + nameDatagroup + ".txt";
  TString pathTextCrossSections = dirCondorPackCurrent + seperatorPath +
                                  "crossSections_" + nameDatagroup + ".txt";

  {
    // std::ifstream infileText;
    // infileText.open(pathTextFileNumbers);
    // UInt_t valNumberFileCurrent;
    // if (!infileText) {
    //   std::cerr << "Unable to open file " << pathTextFileNumbers << std::endl;
    // }
    // while (!infileText) {
    //   infileText >> valNumberFileCurrent;
    //   // if (infileText >> valNumberFileCurrent) {
    //   vNumberFile.push_back(valNumberFileCurrent);
    //   // } else {
    //   //   std::cerr << "One line failed to parse." << std::endl;
    //   // }
    // }
    // infileText.close();
    // infileText.open(pathTextCrossSections);
    // Double_t valCrossSectionCurrent;
    // if (!infileText) {
    //   std::cerr << "Unable to open file " << pathTextCrossSections << std::endl;
    // }
    // while (!infileText) {
    //   infileText >> valCrossSectionCurrent;
    //   // if (infileText >> valCrossSectionCurrentt) {
    //   vCrossSection.push_back(valNumberFileCurrent);
    //   // } else {
    //   //   std::cerr << "One line failed to parse." << std::endl;
    //   // }
    // }
    // infileText.close();
    // TODO
    if (nameCondorPack == "preselect") {
      if (nameDatagroup == "TT") {
        vNumberFile = {
#include "preselect/fileNumbers_TT_withComma.txt"
        };
        vCrossSection = {730.6, 730.6};
      } else if (nameDatagroup == "DYJets") {
        vNumberFile = {
#include "preselect/fileNumbers_DYJets_withComma.txt"
        };
        vCrossSection = {
#include "preselect/crossSections_DYJets_withComma.txt"
        };
      } else if (nameDatagroup == "signal_Mx2-150_Mv-500_Mx1-1_ctau-1") {
        vNumberFile = {4};
        vCrossSection = {1};
      }
    }
  }

  std::function<TString(TString)> funNameTFTemp =
      [nameDatagroup, nameClusterID](TString keyword) -> TString {
    return "output_" + nameDatagroup + "_" + keyword + "_" + nameClusterID +
           "_hist.root";
  };

  TString dirTFInAndSlash = dirCondorPackCurrent + seperatorPath + "dataTest_" +
                            nameDatagroup + seperatorPath;

  std::function<TString(UInt_t)> funPathTFIn = [dirTFInAndSlash, nameDatagroup,
                                                nameClusterID](UInt_t iFile) {
    return dirTFInAndSlash + "output_" + nameDatagroup + "_" + nameClusterID +
           "_" + iFile + ".root";
  };
  std::function<TString(TString)> funNameTFOut =
      [nameDatagroup, nameClusterID](TString nameTT) -> TString {
    return "output_" + nameDatagroup + "_" + nameTT + "_" + nameClusterID +
           "_hist.root";
  };

  std::function<void(Int_t&, Double_t&, Double_t&, TString, TString, TString,
                     TString)>
      adjustHistSetting = [debug](Int_t& nBinCorrect, Double_t& lowerCorrect,
                                  Double_t& upperCorrect, TString nameTT,
                                  TString nameLeafModified,
                                  TString typeNameLeaf, TString titleLeaf) {
        if (typeNameLeaf.Contains("loat") || typeNameLeaf.Contains("ouble")) {
          if (nameLeafModified.EndsWith("_phi")) {
            if (debug) std::cout << "Found _phi" << std::endl;
            lowerCorrect = -TMath::Pi();
            upperCorrect = TMath::Pi();
            return;
          }
          if (nameLeafModified.BeginsWith("Jet_btag")) {
            if (nameLeafModified.BeginsWith("Jet_btagCMVA")) {
              if (debug) std::cout << "Found Jet_btagCMVA" << std::endl;
              lowerCorrect = -1;
              upperCorrect = 1;
            } else {
              if (debug) std::cout << "Found Jet_btag" << std::endl;
              lowerCorrect = 0;
              upperCorrect = 1;
            }
            return;
          }
          if (nameLeafModified.BeginsWith("Jet_ch")) {
            if (debug) std::cout << "Found Jet_ch" << std::endl;
            lowerCorrect = 0;
            upperCorrect = 1;
            return;
          }
          if (nameLeafModified.BeginsWith("Jet_qgl")) {
            if (debug) std::cout << "Found Jet_qgl" << std::endl;
            lowerCorrect = 0;
            upperCorrect = 1;
            return;
          }
          if (titleLeaf.EndsWith("Energy Fraction") ||
              titleLeaf.EndsWith("energy fraction")) {
            if (debug) std::cout << "Found Energy Fraction" << std::endl;
            lowerCorrect = 0;
            upperCorrect = 0;
            return;
          }
        }
        if (typeNameLeaf.Contains("Int") || typeNameLeaf.Contains("int")) {
          if (nameLeafModified.EndsWith("idxJet")) {
            if (debug) std::cout << "Found idxJet" << std::endl;
            lowerCorrect = 0;
            nBinCorrect = upperCorrect - lowerCorrect;
            return;
          }
          if (nameLeafModified.Contains("_rankJet")) {
            if (debug) std::cout << "Found rankJet" << std::endl;
            lowerCorrect = 0;
            nBinCorrect = upperCorrect - lowerCorrect;
            return;
          }
        }
        if (debug)
          std::cout << "No additional settings for " << nameLeafModified
                    << " to apply." << std::endl;
      };
  // mergeToHists(vNameTT, vNumberFile, vCrossSection,
  // funPathTFIn,
  // dirCondorPackCurrent, funNameTFTemp,
  // dirCondorPackCurrent, funNameTFOut,
  // seperatorPath,
  // adjustHistSetting, nullptr, nullptr,
  // toRecreateOutFile, debug, allowMissing);
  merger->vNameTT = vNameTT;
  merger->vNumberFile = vNumberFile;
  merger->vCrossSection = vCrossSection;
  merger->funPathTFIn = funPathTFIn;
  merger->dirTFTemp = dirCondorPackCurrent;
  merger->funNameTFTemp = funNameTFTemp;
  merger->dirTFOut = dirCondorPackCurrent;
  merger->funNameTFOut = funNameTFOut;
  merger->seperatorPath = seperatorPath;
  merger->adjustHistSettingPerLeafTreeExtra = adjustHistSetting;
  merger->toRecreateOutFile = toRecreateOutFile;
  merger->debug = debug;
  merger->allowMissingInputFiles = allowMissing;
  // merger->nLeavesToUseCorrectedTempFileMin = 0;
  merger->funIsToVetoLeaf = [](TString nameTT,
                               TString nameLeafModified) -> Bool_t {
    if (nameLeafModified.Contains("jEntr")) {
      return true;
    }
    return false;
    // return !(nameLeafModified.BeginsWith("GenDMatching_rankJetPassedPt") ||
    //          nameLeafModified.BeginsWith("GenDMatching_JetMatchedId"));
    // return !nameLeafModified.EqualTo("Jet_pt");
  };
  merger->funTitleLeaf = [](TLeaf* leaf) -> TString {
    return leaf->GetBranch()->GetTitle();
  };
  for (TString nameTT : vNameTT) {
    if (nameTT.Contains("Match")) {
      const TString tstrDJetDeltaR = "GenDMathcing_deltaRJet";
      const TString tstrMatchedId = "GenDMatching_JetMatchedId";
      const TString tstrDJetRankPt = "GenDMatching_rankJetPassedPt";
      const TString tstrJetRankPt = "JetMatched_rankJetPassedPt";
      for (UInt_t iD = 0; iD < 4; iD++) {
        auto ptrAnalyzerDJetDeltaR = new HistMerger::LeafAnalyzerDefault;
        ptrAnalyzerDJetDeltaR->SetNameTT(nameTT);
        ptrAnalyzerDJetDeltaR->SetExpressionCustom(
            TString::Format("%s%d", tstrDJetDeltaR.Data(), iD), "Float_t",
            TString::Format("DeltaR between each d quark no. %d and its closest jet", iD),
            TString::Format("%s[%d]", tstrDJetDeltaR.Data(), iD)
        );
        ptrAnalyzerDJetDeltaR->SetHasTarget({tstrDJetDeltaR});
        merger->vAnalyzerCustomByName.push_back((HistMerger::LeafAnalyzerAbstract*)ptrAnalyzerDJetDeltaR);
        auto ptrAnalyzerMatchedId = new HistMerger::LeafAnalyzerDefault;
        ptrAnalyzerMatchedId->SetNameTT(nameTT);
        ptrAnalyzerMatchedId->SetExpressionCustom(
            TString::Format("%s%d", tstrMatchedId.Data(), iD), "Bool_t",
            TString::Format("Whether d-quark no. %d matches a THIN jet", iD),
            TString::Format("%s[%d]", tstrMatchedId.Data(), iD));
        ptrAnalyzerMatchedId->SetHasTarget({tstrMatchedId});
        merger->vAnalyzerCustomByName.push_back(
            (HistMerger::LeafAnalyzerAbstract*)ptrAnalyzerMatchedId);
        auto ptrAnalyzerDJetRankPt = new HistMerger::LeafAnalyzerDefault;
        ptrAnalyzerDJetRankPt->SetNameTT(nameTT);
        ptrAnalyzerDJetRankPt->SetExpressionCustom(
            TString::Format("%s%d", tstrDJetRankPt.Data(), iD), "Int_t",
            TString::Format("Rank (0-indexed) of pt among  passed jet matching "
                            "quark no. %d",
                            iD),
            TString::Format("%s[%d]", tstrDJetRankPt.Data(), iD),
            TString::Format("%s[%d]>0", tstrMatchedId.Data(), iD).Data());
        ptrAnalyzerDJetRankPt->SetHasTarget({tstrDJetRankPt, tstrMatchedId});
        merger->vAnalyzerCustomByName.push_back(
            (HistMerger::LeafAnalyzerAbstract*)ptrAnalyzerDJetRankPt);
        auto ptrAnalyzerRankPtAllMatched =
            new HistMerger::LeafAnalyzerDefault(*ptrAnalyzerDJetRankPt);
        ptrAnalyzerRankPtAllMatched->SetSelection("");
        ptrAnalyzerRankPtAllMatched->SetHasTarget(
            {tstrDJetRankPt}, [tstrMatchedId](TTree* tree) -> Bool_t {
              return tree->GetLeaf(tstrMatchedId) == nullptr;
            });
        merger->vAnalyzerCustomByName.push_back(
            (HistMerger::LeafAnalyzerAbstract*)ptrAnalyzerRankPtAllMatched);
      }
    }
  }
  merger->Run();
}
