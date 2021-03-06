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
#include <numeric>

#include "HistMerger.C"

// std::vector<Int_t> vIota(Int_t iBegin, Int_t iEnd) {
//   std::vector<Int_t> result(iEnd - iBegin);
//   result.resize(result.capacity());
//   std::iota(result.begin(), result.end(), iBegin);
//   return result;
// }

void xAna_monoZ_genhist(const TString nameCondorPack,
                        const TString nameDatagroup,
                        const TString nameClusterID,
                        TString dirCondorPacks = ".",
                        const Bool_t toRecreateOutFile = true,
                        const Bool_t debug = false,
                        const Bool_t allowMissing = false,
                        const Bool_t updateCorrectedFile = false) {
  const TString seperatorPath = "/";
  const TString namesLepton[] = {"Electron", "Muon",
                                 "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon",
                                      "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};

  HistMerger* merger = new HistMerger;

  std::vector<TString> vNameTT;
  vNameTT.clear();
  Bool_t isSignal = nameDatagroup.BeginsWith("signal");
  {
    for (Byte_t i = 0; i < 2; i++) {
      TString nameLeptonCurrent = namesLepton[i];
      vNameTT.push_back((TString) "NumCorrect" + nameLeptonCurrent);
      vNameTT.push_back((TString) "ZMassCutted" + nameLeptonCurrent);
    }
    if (isSignal){
      for (Byte_t i = 0; i < 2; i++) {
        vNameTT.push_back("Gen" + namesLepton[i]);
      }
    } else {
      vNameTT.push_back("Original");
    }
    for (const char* nameJetCamel : {"Jet", "FatJet"}) {
      for (Byte_t i = 0; i < 2; i++) {
        vNameTT.push_back((TString) "PreselectedMatching" + nameJetCamel +
                          namesLepton[i]);
      }
      if (isSignal) for (Byte_t i = 0; i < 2; i++) {
        vNameTT.push_back((TString) "AllMatched" + nameJetCamel +
                          namesLepton[i]);
      }
    }
  }

  std::vector<UInt_t> vNumberFile;
  vNumberFile.clear();
  std::vector<Double_t> vCrossSection;
  vCrossSection.clear();

  std::function<Int_t (TDirectory *tdirIn, Long64_t &nEntriesOriginal)> funRefgetNEntriesOriginal = [isSignal](TDirectory *tdirIn, Long64_t &nEntriesOriginal) {
    TObject* tvdRaw = tdirIn->Get("tvdNEntryOriginal");
    if (!tvdRaw || tvdRaw->IsZombie()) {
      return 1;
    }
    TVectorD* tvd = static_cast<TVectorD *>(tvdRaw);
    if (tvd->IsZombie() || tvd->GetNoElements() < 1) {
      nEntriesOriginal = 0;
      return 1;
    }
    nEntriesOriginal = (*((TVectorD *)tdirIn->Get("tvdNEntryOriginal")))[0];
    if (isSignal) {
      nEntriesOriginal /= 3;
    }
    return 0;
  };

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
    //   std::cerr << "Unable to open file " << pathTextFileNumbers <<
    //   std::endl;
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
    //   std::cerr << "Unable to open file " << pathTextCrossSections <<
    //   std::endl;
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
      } else if (nameDatagroup == "signal_Mx2-1_Mv-500_Mx1-0p1_ctau-1") {
        vNumberFile = {3};
        vCrossSection = {1};
      } else if (nameDatagroup == "signal_Mx2-50_Mv-500_Mx1-1_ctau-10") {
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
      adjustHistSetting = [debug, nameCondorPack](Int_t& nBinsCorrect, Double_t& lowerCorrect,
                                  Double_t& upperCorrect, TString nameTT,
                                  TString nameLeafModified,
                                  TString typeNameLeaf, TString titleLeaf) {
        if (typeNameLeaf.Contains("loat") || typeNameLeaf.Contains("ouble")) {
          const Int_t binDensityUsual = 1;
          const Int_t binDensityHigh = 100;
          const Int_t binDensityMedium = 10;
          if (nameLeafModified.Contains("_eta")) {
            if (debug) std::cout << "Found _eta" << std::endl;
            nBinsCorrect = TMath::Nint((upperCorrect - lowerCorrect) * binDensityMedium);
            return;
          }
          if (nameLeafModified.Contains("_phi")) {
            if (debug) std::cout << "Found _phi" << std::endl;
            lowerCorrect = TMath::Floor(-TMath::Pi() * binDensityMedium) / binDensityMedium;
            upperCorrect = TMath::Ceil(TMath::Pi() * binDensityMedium) / binDensityMedium;
            nBinsCorrect = TMath::Nint((upperCorrect - lowerCorrect) * binDensityMedium);
            return;
          }
          if (nameLeafModified.Contains("_deltaR")) {
            if (debug) std::cout << "Found _deltaR" << std::endl;
            lowerCorrect = 0.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityMedium;
            return;
          }
          if (nameLeafModified.Contains("_area")) {
            if (debug) std::cout << "Found _area" << std::endl;
            lowerCorrect = 0.;
            upperCorrect = TMath::Min(1000., upperCorrect);
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            return;
          }
          if (nameLeafModified.Contains("_mass")) {
            if (debug) std::cout << "Found _mass" << std::endl;
            lowerCorrect = 0.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityMedium;
            return;
          }
          if (nameLeafModified.Contains("_significance")) {
            if (debug) std::cout << "Found _significance" << std::endl;
            lowerCorrect = 0.;
            upperCorrect = TMath::Min(1000., upperCorrect);
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
          }
          if (nameLeafModified.Contains("_msoftdrop")) {
            if (debug) std::cout << "Found _msoftdropt" << std::endl;
            lowerCorrect = 0.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityMedium;
          }
          if (nameLeafModified.Contains("_btag")) {
            if (nameLeafModified.Contains("_btagCMVA") ||
                nameLeafModified.Contains("_btagHbb")) {
              if (debug) std::cout << "Found _btagCMVA" << std::endl;
              lowerCorrect = -1;
              upperCorrect = 1;
              nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            } else {
              if (debug) std::cout << "Found _btag" << std::endl;
              lowerCorrect = 0.;
              upperCorrect = 1.;
              nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            }
            return;
          }
          if (nameLeafModified.Contains("_deepTag")) {
            if (debug) std::cout << "Found _deepTag" << std::endl;
            lowerCorrect = 0.;
            upperCorrect = 1.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            return;
          }
          if (nameLeafModified.Contains("_ch")) {
            if (debug) std::cout << "Found _ch" << std::endl;
            lowerCorrect = 0.;
            upperCorrect = 1.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            return;
          }
          if (nameLeafModified.Contains("_qgl")) {
            if (debug) std::cout << "Found _qgl" << std::endl;
            lowerCorrect = 0.;
            upperCorrect = 1.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            return;
          }
          if (nameLeafModified.Contains("_puIdDisc")) {
            // Pileup ID discriminant
            if (debug) std::cout << "Found _puIdDisc" << std::endl;
            lowerCorrect = -1.;
            upperCorrect = 1.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            return;
          }
          if (nameLeafModified.Contains("_lsf")) {
            //
            if (debug) std::cout << "Found _lsf" << std::endl;
            lowerCorrect = -1.;
            upperCorrect = 1.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            return;
          }
          if (titleLeaf.Contains("Energy Fraction") ||
              titleLeaf.Contains("Energy fraction") ||
              titleLeaf.Contains("energy fraction")) {
            if (debug) std::cout << "Found Energy Fraction" << std::endl;
            lowerCorrect = 0.;
            upperCorrect = 1.;
            nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
            return;
          }
          if (titleLeaf.Contains("regression")) {
            if (titleLeaf.Contains("res ") || titleLeaf.EndsWith("res")
                || titleLeaf.Contains("resolution")) {
              lowerCorrect = 0.;
              nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
              return;
            }
            if (titleLeaf.Contains("correction")) {
              nBinsCorrect = (upperCorrect - lowerCorrect) * binDensityHigh;
              return;
            }
          }
        }
        if (typeNameLeaf.Contains("Int") || typeNameLeaf.Contains("int")) {
          if (nameLeafModified.EndsWith("idxJet") || nameLeafModified.EndsWith("idxFatJet")) {
            if (debug) std::cout << "Found idxJet" << std::endl;
            lowerCorrect = -0.5;
            if (upperCorrect <= lowerCorrect) {
              upperCorrect = lowerCorrect + 1;
            }
            nBinsCorrect = upperCorrect - lowerCorrect;
            return;
          }
          if (nameLeafModified.Contains("_rank")) {
            if (debug) std::cout << "Found rank" << std::endl;
            lowerCorrect = -0.5;
            if (upperCorrect <= lowerCorrect) {
              upperCorrect = lowerCorrect + 1;
            }
            nBinsCorrect = upperCorrect - lowerCorrect;
            return;
          }
        }
        if (false && nameCondorPack.EqualTo("preselect")) {
          if (typeNameLeaf.Contains("loat") || typeNameLeaf.Contains("ouble")) {
            if (nameLeafModified.BeginsWith("Electron") && nameLeafModified.Contains("_pt")) {
              if (lowerCorrect < 200) {
                upperCorrect = TMath::Min(upperCorrect, 200.);
              }
            }
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
  merger->funRefgetNEntriesOriginal = funRefgetNEntriesOriginal;
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
  merger->updateCorrectedFile = updateCorrectedFile;
  // merger->nLeavesToUseCorrectedTempFileMin = 0;
  merger->nLeavesToUseCorrectedTempFileMin = 30;
  merger->funIsToVetoLeaf = [](TString nameTT,
                               TString nameLeafModified) -> Bool_t {
    if (nameLeafModified.Contains("jEntr")) {
      return true;
    }
    return false;
    // return !(nameLeafModified.BeginsWith("GenDMatching_rankJetPassedPt") ||
    //          nameLeafModified.BeginsWith("GenDMatching_JetMatchedId"));
    // return !nameLeafModified.EqualTo("GenDMathcing_deltaRJet") && !nameLeafModified.EqualTo("Jet_neEmEF");
    // return !nameLeafModified.EqualTo("Jet_neEmEF");
    // return !nameLeafModified.EqualTo("nGenDMatched");
    // return !nameLeafModified.EqualTo("FatJetMatched_rankFatJetPassedPt");
    // return !nameLeafModified.EqualTo("FatJet_nBHadrons");
  };
  merger->funTitleLeaf = [](TLeaf* leaf) -> TString {
    return leaf->GetBranch()->GetTitle();
  };
  if (nameCondorPack.EqualTo("preselect")) { 
    for (TString nameTT : vNameTT) {
      for (Byte_t i=0; i<2; i++) {
        for (Byte_t isPassedPtEta=0; isPassedPtEta < 2; isPassedPtEta++) for (Byte_t isExact = 0; isExact < 2; isExact++) {
          auto analyzer = new HistMerger::LeafAnalyzerDefault;
          analyzer->SetNameTT(nameTT);
          const TString nameVarIn = "n" + namesLepton[i] + (isPassedPtEta ? "PassedPtEta" : "");
          const TString nameOperatorOut = isExact ? "==" : ">=";
          analyzer->SetExpressionCustom((isExact ? "have2" : "have2More") + namesLepton[i] + (isPassedPtEta ? "PassedPtEta" : ""), "Bool_t",
          nameVarIn + " " + nameOperatorOut + " 2",
          nameVarIn + nameOperatorOut + "2");
          analyzer->SetHasTarget({nameVarIn});
          merger->vAnalyzerCustomByName.push_back(analyzer);
        }
      }
      {
        auto analyzer = new HistMerger::LeafAnalyzerDefault;
        analyzer->SetNameTT(nameTT);
        analyzer->SetExpressionCustom("have2MoreJetPassed", "Bool_t", "nJetPassed >= 2",
            "nJetPassed>=2");
        analyzer->SetHasTarget({"nJetPassed"});
        merger->vAnalyzerCustomByName.push_back(analyzer);
      }
      {
        auto analyzer = new HistMerger::LeafAnalyzerDefault;
        analyzer->SetNameTT(nameTT);
        analyzer->SetExpressionCustom("have2MoreFatJetPassed", "Bool_t", "nFatJetPassed >= 2",
            "nFatJetPassed>=2");
        analyzer->SetHasTarget({"nFatJetPassed"});
        merger->vAnalyzerCustomByName.push_back(analyzer);
      }
      if (nameTT.Contains("Match")) {
        // const TString tstrDJetDeltaR = "GenDMatching_deltaRJet";
        // const TString tstrMatchedId = "GenDMatching_JetMatchedId";
        // const TString tstrDJetRankPt = "GenDMatching_rankJetPassedPt";
        // const TString tstrJetRankPt = "JetMatched_rankJetPassedPt";
        // for (UInt_t iD = 0; iD < 4; iD++) {
        //   auto ptrAnalyzerDJetDeltaR = new HistMerger::LeafAnalyzerDefault;
        //   ptrAnalyzerDJetDeltaR->SetNameTT(nameTT);
        //   ptrAnalyzerDJetDeltaR->SetExpressionCustom(
        //       TString::Format("%s%d", tstrDJetDeltaR.Data(), iD), "Float_t",
        //       TString::Format(
        //           "DeltaR between each d quark no. %d and its closest jet", iD),
        //       TString::Format("%s[%d]", tstrDJetDeltaR.Data(), iD));
        //   ptrAnalyzerDJetDeltaR->SetHasTarget({tstrDJetDeltaR});
        //   merger->vAnalyzerCustomByName.push_back(
        //       (HistMerger::LeafAnalyzerAbstract*)ptrAnalyzerDJetDeltaR);
        //   auto ptrAnalyzerMatchedId = new HistMerger::LeafAnalyzerDefault;
        //   ptrAnalyzerMatchedId->SetNameTT(nameTT);
        //   ptrAnalyzerMatchedId->SetExpressionCustom(
        //       TString::Format("%s%d", tstrMatchedId.Data(), iD), "Bool_t",
        //       TString::Format("Whether d-quark no. %d matches a THIN jet", iD),
        //       TString::Format("%s[%d]", tstrMatchedId.Data(), iD));
        //   ptrAnalyzerMatchedId->SetHasTarget({tstrMatchedId});
        //   merger->vAnalyzerCustomByName.push_back(
        //       (HistMerger::LeafAnalyzerAbstract*)ptrAnalyzerMatchedId);
        //   auto ptrAnalyzerDJetRankPt = new HistMerger::LeafAnalyzerDefault;
        //   ptrAnalyzerDJetRankPt->SetNameTT(nameTT);
        //   ptrAnalyzerDJetRankPt->SetExpressionCustom(
        //       TString::Format("%s%d", tstrDJetRankPt.Data(), iD), "Int_t",
        //       TString::Format("Rank (0-indexed) of pt among  passed jet matching "
        //                       "quark no. %d",
        //                       iD),
        //       TString::Format("%s[%d]", tstrDJetRankPt.Data(), iD),
        //       TString::Format("%s[%d]>0", tstrMatchedId.Data(), iD).Data());
        //   ptrAnalyzerDJetRankPt->SetHasTarget({tstrDJetRankPt, tstrMatchedId});
        //   merger->vAnalyzerCustomByName.push_back(
        //       (HistMerger::LeafAnalyzerAbstract*)ptrAnalyzerDJetRankPt);
        //   auto ptrAnalyzerRankPtAllMatched =
        //       new HistMerger::LeafAnalyzerDefault(*ptrAnalyzerDJetRankPt);
        //   ptrAnalyzerRankPtAllMatched->SetSelection("");
        //   ptrAnalyzerRankPtAllMatched->SetHasTarget(
        //       {tstrDJetRankPt}, [tstrMatchedId](TTree* tree) -> Bool_t {
        //         return tree->GetLeaf(tstrMatchedId) == nullptr;
        //       });
        //   merger->vAnalyzerCustomByName.push_back(
        //       (HistMerger::LeafAnalyzerAbstract*)ptrAnalyzerRankPtAllMatched);
        // }
        auto ptrAnalyzerJetRankLast = new HistMerger::LeafAnalyzerDefault;
        ptrAnalyzerJetRankLast->SetNameTT(nameTT);
        ptrAnalyzerJetRankLast->SetExpressionCustom(
            "maxJetMatched_rankJetPassedPt",
            "Int_t", "Last rank of pt among passed jets",
            "(nJetMatched>0)?JetMatched_rankJetPassedPt[nJetMatched-1]:-10");
        ptrAnalyzerJetRankLast->SetHasTarget({"JetMatched_rankJetPassedPt", "nJetMatched"});
        merger->vAnalyzerCustomByName.push_back(ptrAnalyzerJetRankLast);
        auto ptrAnalyzerFatJetRankLast = new HistMerger::LeafAnalyzerDefault;
        ptrAnalyzerFatJetRankLast->SetNameTT(nameTT);
        ptrAnalyzerFatJetRankLast->SetExpressionCustom(
            "maxFatJetMatched_rankFatJetPassedPt",
            "Int_t", "Last rank of pt among passed jets",
            "(nFatJetMatched>0)?FatJetMatched_rankFatJetPassedPt[nFatJetMatched-1]:-10");
        ptrAnalyzerFatJetRankLast->SetHasTarget({"FatJetMatched_rankFatJetPassedPt", "nFatJetMatched"});
        merger->vAnalyzerCustomByName.push_back(ptrAnalyzerFatJetRankLast);
      }
    }
    merger->pushCustomAnalyzersWhenRun = [debug, namesLepton, namesLeptonLower](
        TTree *tree, const std::vector<HistMerger::LeafAnalyzerAbstract*> vAnalyzerLeafCustom,
        const std::vector<HistMerger::LeafAnalyzerAbstract*> vAnalyzerLeaf,
        const UInt_t nAnalyzerLeafOld,
        const std::function<void(HistMerger::LeafAnalyzerAbstract* analyzerNew)> pushbackNewAnalyzer) {
      for (auto citerVAnalyzerLeaf = vAnalyzerLeaf.cbegin() + nAnalyzerLeafOld;
          citerVAnalyzerLeaf != vAnalyzerLeaf.cend(); citerVAnalyzerLeaf++) {
        if (debug) std::cout << "On normal analyzer " << " index " << std::distance(vAnalyzerLeaf.cbegin(), citerVAnalyzerLeaf) << " (" << *citerVAnalyzerLeaf << ")" << std::endl;
        const TString nameLeafModified = (*citerVAnalyzerLeaf)->GetNameLeafModified();
        // if (nameLeafModified.BeginsWith("n")
        //     || nameLeafModified.BeginsWith("is") || nameLeafModified.BeginsWith("are")
        //     || nameLeafModified.BeginsWith("has") || nameLeafModified.BeginsWith("have")) {
        //   if (!TString(tree->GetName()).Contains("Fat") && !nameLeafModified.EqualTo("nJetPassed")) {
        //     for (Int_t nJetPassedExpected=1; nJetPassedExpected <= 8; nJetPassedExpected++) {
        //       auto analyzer = new HistMerger::LeafAnalyzerDefault;
        //       analyzer->SetExpressionCustom(
        //           nameLeafModified + "WithnJetPassedEq" + nJetPassedExpected, (*citerVAnalyzerLeaf)->GetTypeNameLeaf(),
        //           (*citerVAnalyzerLeaf)->GetTitleLeaf() + " with nJetPassed==" + nJetPassedExpected,
        //           nameLeafModified, Form("nJetPassed==%d", nJetPassedExpected));
        //       analyzer->SetHasTarget({nameLeafModified, "nJetPassed"});
        //       pushbackNewAnalyzer(analyzer);
        //     }
        //   }
        //   if (TString(tree->GetName()).Contains("Fat") && !nameLeafModified.EqualTo("nFatJetPassed")) {
        //     for (Int_t nFatJetPassedExpected=1; nFatJetPassedExpected <= 6; nFatJetPassedExpected++) {
        //       auto analyzer = new HistMerger::LeafAnalyzerDefault;
        //       analyzer->SetExpressionCustom(
        //           nameLeafModified + "WithnFatJetPassedEq" + nFatJetPassedExpected, (*citerVAnalyzerLeaf)->GetTypeNameLeaf(),
        //           (*citerVAnalyzerLeaf)->GetTitleLeaf() + " with nFatJetPassed==" + nFatJetPassedExpected,
        //           nameLeafModified, Form("nFatJetPassed==%d",nFatJetPassedExpected));
        //       analyzer->SetHasTarget({nameLeafModified, "nFatJetPassed"});
        //       pushbackNewAnalyzer(analyzer);
        //     }
        //   }
        // }
        if (nameLeafModified.BeginsWith("Electron_") && nameLeafModified.Contains("_WP")) {
          {
            auto analyzer = new HistMerger::LeafAnalyzerDefault;
            analyzer->SetExpressionCustom(
                TString("n") + nameLeafModified,
                "Int_t",
                "Number of electron passing " + (*citerVAnalyzerLeaf)->GetTitleLeaf(),
                "Sum$(" + nameLeafModified + ")"
            );
            analyzer->SetHasTarget({nameLeafModified});
            pushbackNewAnalyzer(analyzer);
          }
          {
            auto analyzer = new HistMerger::LeafAnalyzerDefault;
            analyzer->SetExpressionCustom(
                TString("have2More") + nameLeafModified,
                "Bool_t",
                "Whether there are more than 2 electron passing " + (*citerVAnalyzerLeaf)->GetTitleLeaf(),
                "Sum$(" + nameLeafModified + ")>=2"
            );
            analyzer->SetHasTarget({nameLeafModified});
            pushbackNewAnalyzer(analyzer);
          }
          {
            auto analyzer = new HistMerger::LeafAnalyzerDefault;
            analyzer->SetExpressionCustom(
                TString("have2") + nameLeafModified,
                "Bool_t",
                "Whether there are 2 electron passing " + (*citerVAnalyzerLeaf)->GetTitleLeaf(),
                "Sum$(" + nameLeafModified + ")==2"
            );
            analyzer->SetHasTarget({nameLeafModified});
            pushbackNewAnalyzer(analyzer);
          }
        }
        if (nameLeafModified.BeginsWith("Electron_") || nameLeafModified.BeginsWith("Muon_")) {
          Bool_t isMuon = nameLeafModified.BeginsWith("Muon_");
          {
            auto analyzer = new HistMerger::LeafAnalyzerDefault;
            analyzer->SetExpressionCustom(
                namesLepton[isMuon] +  "NumCorrect" + nameLeafModified(isMuon ? 4 : 8, nameLeafModified.Length()),
                (*citerVAnalyzerLeaf)->GetTypeNameLeaf(),
                (*citerVAnalyzerLeaf)->GetTitleLeaf(),
                nameLeafModified,
                Form("(%s_idxNumCorrect>=0)&&(%s_idxPassedPtEta==%s_idxNumCorrect)", namesLepton[isMuon].Data(), namesLepton[isMuon].Data(), namesLepton[isMuon].Data()));
            analyzer->SetHasTarget({nameLeafModified, namesLepton[isMuon] + "_idxNumCorrect", namesLepton[isMuon] + "_idxPassedPtEta"});
            pushbackNewAnalyzer(analyzer);
          }
        }
        Bool_t isFatJet = false;
        if (nameLeafModified.BeginsWith("FatJet_")) {
          isFatJet = true;
        } else if (nameLeafModified.BeginsWith("Jet_")) {
        } else {
          continue;
        }
        if (((TString)tree->GetName()).Contains("Match")) {
          auto analyzer = new HistMerger::LeafAnalyzerDefault;
          analyzer->SetExpressionCustom((TString)(isFatJet ? "FatJet" : "Jet") + "Matched" + nameLeafModified(isFatJet ? 6 : 3, nameLeafModified.Length()),
          (*citerVAnalyzerLeaf)->GetTypeNameLeaf(),
          (*citerVAnalyzerLeaf)->GetTitleLeaf() + " that is matched",
          TString::Format("%s[%sMatched_rank%sPassedPt]", nameLeafModified.Data(), isFatJet ? "FatJet" : "Jet", isFatJet ? "FatJet" : "Jet"));
          analyzer->SetHasTarget({nameLeafModified, TString::Format("%sMatched_rank%sPassedPt", isFatJet ? "FatJet" : "Jet", isFatJet ? "FatJet" : "Jet")});
          pushbackNewAnalyzer(analyzer);
        }
        // const UInt_t rankUpperMax = isFatJet ? 6 : 8;
        // for (UInt_t rankUpper = 1; rankUpper <= rankUpperMax; rankUpper++) {
        {
          UInt_t rankUpper = isFatJet ? 2 : 6;
          auto analyzer = new HistMerger::LeafAnalyzerDefault;
          analyzer->SetExpressionCustom((TString)(isFatJet ? "FatJet" : "Jet") + "First" + rankUpper + nameLeafModified(isFatJet ? 6 : 3, nameLeafModified.Length()),
              (*citerVAnalyzerLeaf)->GetTypeNameLeaf(),
              Form("First %d -- %s", rankUpper, (*citerVAnalyzerLeaf)->GetTitleLeaf().Data()),
              nameLeafModified, Form("Iteration$<%d", rankUpper));
          analyzer->SetHasTarget({nameLeafModified});
          pushbackNewAnalyzer(analyzer);
        }
        // }
      }
    };
  }
  merger->Run();
}
