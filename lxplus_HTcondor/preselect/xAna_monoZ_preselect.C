#include <Rtypes.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStyle.h>
#include <TVectorFfwd.h>
#include <TVectorT.h>

//#include <cstring>  // std::memcpy
#include <cstdio>  //printf
#include <iostream>
#include <numeric>
#include <vector>

#include "BranchMounter.h"
#include "BranchMounterForUntuplizer.h"
#include "untuplizer.h"

// template <typename TypeData>
// TBranch* branchAdder(TTree* tree, const char* name, TypeData* address, const
// char* typestring=nullptr, Int_t buffsize=32000) {
//   if (typestring == nullptr) {
//     return tree->Branch(name, address);
//   } else {
//     return tree->Branch(name, address, ((TString)name + "/" + typestring),
//     buffsize);
//   }
// }

// Bool_t getIsMatchAt(TString target, TString pattern, Ssiz_t startExpected) {
//   TSubString substring = target.SubString(pattern);
//   return substring.Start() == startExpected;
// }

/// Preselection
void xAna_monoZ_preselect(
    std::string inputFile,  //< Path of input file
    TString outputFileTree,
    // std::string outputFileHead,  //< Output filename head ex: "preselected"
    // std::string outputFileVar,   //< Variables of the data ex:
    //                              //"Mx2-150_Mv-500_Mx1-1_ctau-1"
    // std::string outputFileTail,  //< Output filename tail (date) ex:
    // "20200730"
    bool toRecreateOutFile = true, bool debug = false) {
  const Double_t massZ = 91.1876;  //< static mass of Z (constant)
  const Double_t massElectron =
      0.0005109989461;          //< static mass of electron (constant)
  const Double_t massDown = 0.0048;
  const Int_t pdgZ = 23;        //< pdgid of Z
  const Int_t pdgZp = 55;       //< pdgid of Z'
  const Int_t pdgX2 = 18;       //< pdgid of x2
  const Int_t pdgX1 = 5000522;  //< pdgid of x1
  const Int_t pdgDown = 1;      //< pdgid of down quark
  const Int_t pdgElectron = 11;
  const Int_t pdgMuon = 13;
  const Int_t pdgTau = 15;

  /// path to the output file to store the trees
  // TString outputFileTree = (TString) "output_" + outputFileHead + "_" +
  //                          outputFileVar + "_" + outputFileTail + "_tree" +
  //                          ".root";

  /// output fle to store the tree
  TFile* outFileTree = TFile::Open(outputFileTree.Data(),
                                   toRecreateOutFile ? "recreate" : "update");

  const char* nameTreeIn = "Events";              //< the input tree
  TreeReader data(inputFile.data(), nameTreeIn);  //< the TreeReader
  const Long64_t nEntry = data.GetEntriesFast();  //< number of entries
  // if (debug) std::cout<< "nEntry: " << nEntry << std::endl;

  const TString namesLepton[] = {"Electron", "Muon",
                                 "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon",
                                      "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};
  // TTree* ttGenElectron = new TTree("GenElectron", "GEN-level electron
  // events");

  TVectorD tvdNEntry(1);
  tvdNEntry[0] = nEntry;
  tvdNEntry.Write("tvdNEntryOriginal");

  /// Tree for GEN-correct  electron/muon/tau
  TTree* arrTTGen[3];

  /// Tree for events with correct numbers of electron/muon
  TTree* arrTTNumCorrect[2];

  /// Tree for events whose electron/muon pairs pass Z mass cut
  TTree* arrTTZMassCutted[2];

  /// Tree for the jet match results of events passing the preselections
  TTree* arrTTPreselectedMatching[2];

  /// Tree for the events whose four d-quarks are all matched to a jet.
  TTree* arrTTAllMatched[2];
  for (Byte_t i = 0; i < 3; i++) {
    arrTTGen[i] =
        new TTree((TString) "Gen" + namesLepton[i],
                  (TString) "GEN-level " + namesLeptonLower[i] + " events");
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTNumCorrect[i] = new TTree(
        (TString) "NumCorrect" + namesLepton[i],
        (TString) "Events with correct number of " + namesLeptonLower[i]);
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTZMassCutted[i] =
        new TTree((TString) "ZMassCutted" + namesLepton[i],
                  (TString) "Z mass cutted " + namesLeptonLower[i] + " events");
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTPreselectedMatching[i] =
        new TTree((TString) "PreselectedMatching" + namesLepton[i],
                  (TString) "Jet matching results of preselected " +
                      namesLeptonLower[i] + " events");
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTAllMatched[i] =
        new TTree((TString) "AllMatched" + namesLepton[i],
                  (TString) "Jet matching results of the all-matched " +
                      namesLeptonLower[i] + " events");
  }

  // TTree* ttNumCorrect = new TTree(
  //     "NumCorrect",
  //     "Variables of entries with 2 electrons/2 muons");  //< tree filled when
  //                                                        //the number is
  //                                                        right
  // TTree* ttZMassCutted = new TTree(
  //     "ZMassCutted",
  //     "Preselected variables");  //< tree filled when passing Z mass cut

  /// Index of entry, bounded to the trees.
  Long64_t jEntry;
  // ttNumCorrect->Branch("jEntryOriginal", &jEntry);
  // ttZMassCutted->Branch("jEntryOriginal", &jEntry);
  for (Byte_t i = 0; i < 3; i++) {
    arrTTGen[i]->Branch("jEntry", &jEntry);
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTNumCorrect[i]->Branch("jEntry", &jEntry);
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTZMassCutted[i]->Branch("jEntry", &jEntry);
  }
  Bool_t haveAllGenDMaching;
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "haveAllGenD";
    TString title = "Whether there are 4 d quarks generated by chi1";
    arrTTGen[i]
        ->Branch(name, &haveAllGenDMaching, name + "/O", 1)
        ->SetTitle(title);
  }
  std::vector<Int_t> vGenDMatchingIdx(4, 0);
  for (Byte_t i = 0; i < 2; i++) {
    const TString name = "GenDMatching_idx";
    const TString title = "Indexes of the d quarks";
    arrTTGen[i]->Branch(name, &vGenDMatchingIdx)->SetTitle(title);
  }
  std::vector<Float_t> vGenDMatchingPt(4, 0);
  std::vector<Float_t> vGenDMatchingEta(4, 0);
  std::vector<Float_t> vGenDMatchingPhi(4, 0);
  std::vector<Float_t> vGenDPairMatchingPt(2, 0);
  std::vector<Float_t> vGenDPairMatchingEta(2, 0);
  std::vector<Float_t> vGenDPairMatchingPhi(2, 0);
  for (Byte_t i = 0; i < 2; i++) {
    const char* nameTemplate = "GenD%sMatching_%s";
    const char* titleTemplate = "%s of each %s";
    arrTTGen[i]->Branch(TString::Format(nameTemplate, "", "pt"), &vGenDMatchingPt)->SetTitle(TString::Format(titleTemplate, "Pt", "d quark"));
    arrTTPreselectedMatching[i]->Branch(TString::Format(nameTemplate, "", "eta"), &vGenDMatchingPt)->SetTitle(TString::Format(titleTemplate, "Eta", "d quark"));
    arrTTAllMatched[i]->Branch(TString::Format(nameTemplate, "", "phi"), &vGenDMatchingPt)->SetTitle(TString::Format(titleTemplate, "Phi", "d quark"));
    arrTTGen[i]->Branch(TString::Format(nameTemplate, "Pair", "pt"), &vGenDMatchingPt)->SetTitle(TString::Format(titleTemplate, "Pt", "d-pair"));
    arrTTPreselectedMatching[i]->Branch(TString::Format(nameTemplate, "Pair", "eta"), &vGenDMatchingPt)->SetTitle(TString::Format(titleTemplate, "Eta", "d-pair"));
    arrTTAllMatched[i]->Branch(TString::Format(nameTemplate, "Pair", "phi"), &vGenDMatchingPt)->SetTitle(TString::Format(titleTemplate, "Phi", "d-pair"));
  }
  std::vector<Bool_t> vGenDMatchingIdPassed(4, false);
  UInt_t nGenDMatchingPassed = 0;
  Bool_t areAllGenDMatchingPassed = false;
  std::vector<Bool_t> vGenDPairMatchingIdPassed(2, false);
  UInt_t nGenDPairMatchingPassed = 0;
  Bool_t areAllGenDPairMatchingPassed;
  for (Byte_t i = 0; i < 2; i++) {
    const char* nameIdTemplate = "GenD%sMatching_idPassed";
    const char* titleIdTemplate =
        "Whether each of the %s pass the preselections";
    const char* nameNTemplate = "nGenD%sMatchingPassed";
    const char* titleNTemplate = "Number of %ss passing the preselections";
    const char* nameAreAllTemplate = "areAllGenD%sMatchingPassed";
    const char* titleAreAllTemplate =
        "Whether all the %ss have passed the preselections";
    arrTTGen[i]->Branch(TString::Format(nameIdTemplate, ""), &vGenDMatchingIdPassed)->SetTitle(TString::Format(titleIdTemplate, "d quark"));
    arrTTGen[i]
        ->Branch(TString::Format(nameNTemplate, ""), &nGenDMatchingPassed, TString::Format(nameNTemplate, "") + "/I", __SIZEOF_INT__)
        ->SetTitle(TString::Format(titleNTemplate, "d quark"));
    arrTTGen[i]
        ->Branch(TString::Format(nameAreAllTemplate, ""), &areAllGenDMatchingPassed, TString::Format(nameAreAllTemplate, "") + "/O", 1)
        ->SetTitle(TString::Format(titleAreAllTemplate, "d quark"));
    arrTTGen[i]->Branch(TString::Format(nameIdTemplate, "Pair"), &vGenDPairMatchingIdPassed)->SetTitle(TString::Format(titleIdTemplate, "d-pair"));
    arrTTGen[i]
        ->Branch(TString::Format(nameNTemplate, "Pair"), &nGenDPairMatchingPassed, TString::Format(nameNTemplate, "") + "/I", __SIZEOF_INT__)
        ->SetTitle(TString::Format(titleNTemplate, "d-pair"));
    arrTTGen[i]
        ->Branch(TString::Format(nameAreAllTemplate, "Pair"), &areAllGenDPairMatchingPassed, TString::Format(nameAreAllTemplate, "") + "/O", 1)
        ->SetTitle(TString::Format(titleAreAllTemplate, "d-pair"));
  }
  Bool_t arrHasGenLepton[3];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "hasGen" + namesLepton[i] + "Pair";
    TString title = "Whether this is a Z-" + namesLeptonNota[i] +
                    namesLeptonNota[i] + " event";
    arrTTGen[i]
        ->Branch(name, &(arrHasGenLepton[i]), name + "/O", 1)
        ->SetTitle(title);
  }
  Bool_t arrRecoIsLeptonNumCorrect[2];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "is" + namesLepton[i] + "NumCorrect";
    TString title = "Whether there are the right number of " +
                    namesLeptonLower[i] + " in RECO-level";
    arrTTGen[i]
        ->Branch(name, &(arrRecoIsLeptonNumCorrect[i]), name + "/O", 1)
        ->SetTitle(title);
  }
  Float_t arrLeptonPair_mass[2];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = namesLepton[i] + "Pair_mass";
    TString title = "The mass of the " + namesLeptonLower[i] + " pair";
    arrTTGen[i]
        ->Branch(name, &(arrLeptonPair_mass[i]), name + "/F", __SIZEOF_FLOAT__)
        ->SetTitle(title);
    arrTTNumCorrect[i]
        ->Branch(name, &(arrLeptonPair_mass[i]), name + "/F", __SIZEOF_FLOAT__)
        ->SetTitle(title);
    arrTTZMassCutted[i]
        ->Branch(name, &(arrLeptonPair_mass[i]), name + "/F", __SIZEOF_FLOAT__)
        ->SetTitle(title);
  }
  Bool_t arrIsPassingZMassCut[2];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = namesLepton[i] + "Pair_isPassingZMassCut";
    TString title = "Whether is the mass of the " + namesLeptonLower[i] +
                    " pair close to Z\'s";
    arrTTGen[i]
        ->Branch(name, &(arrIsPassingZMassCut[i]), name + "/O", 1)
        ->SetTitle(title);
    arrTTNumCorrect[i]
        ->Branch(name, &(arrIsPassingZMassCut[i]), name + "/O", 1)
        ->SetTitle(title);
  }

  const UInt_t nDQuarksExpected = 4;

  std::vector<Bool_t> vDJetMatchedId(nDQuarksExpected, false);
  for (Byte_t i = 0; i < 2; i++) {
    arrTTPreselectedMatching[i]
        ->Branch("GenDMatching_JetMatchedId", &vDJetMatchedId)
        ->SetTitle("Whether the d-quarks matches a thin jet");
  }
  UInt_t nDJetMatched;
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "nDJetMatched";
    arrTTPreselectedMatching[i]
        ->Branch(name, &nDJetMatched, name + "/I", __SIZEOF_INT__)
        ->SetTitle("Number of d-quarks matching a thin jet");
  }
  Bool_t areAllDJetMatched;
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "areAllDJetMatched";
    arrTTPreselectedMatching[i]
        ->Branch(name, &areAllDJetMatched, name + "/O", __SIZEOF_INT__)
        ->SetTitle("Whether all the d-quarks matches a thin jet");
  }
  UInt_t nJetMatched;
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "nJetMatched";
    TString title = "Number of unique thin jets matching a d-quark";
    arrTTPreselectedMatching[i]
        ->Branch(name, &nJetMatched, name + "/I", __SIZEOF_INT__)
        ->SetTitle(title);
    arrTTAllMatched[i]
        ->Branch(name, &nJetMatched, name + "/I", __SIZEOF_INT__)
        ->SetTitle(title);
  }
  auto lambdaFillTheTrees = [&]() -> void {
    const Bool_t toCollectGen = true;
    const Bool_t toCollectNumCorrect = false;
    const Bool_t toCollectZMassCutted = true;
    const Bool_t toCollectJetMatching = true;
    const Bool_t toCollectAllMatched = true;
    if (toCollectGen)
      for (Byte_t i = 0; i < 3; i++) {
        if (arrHasGenLepton[i]) {
          arrTTGen[i]->Fill();
        }
      }
    for (Byte_t i = 0; i < 2; i++) {
      if (toCollectNumCorrect && arrRecoIsLeptonNumCorrect[i]) {
        arrTTNumCorrect[i]->Fill();
      }
      if (toCollectZMassCutted && arrIsPassingZMassCut[i]) {
        arrTTZMassCutted[i]->Fill();
      }
      if (areAllGenDMatchingPassed && toCollectJetMatching &&
          arrIsPassingZMassCut[i]) {
        arrTTPreselectedMatching[i]->Fill();
      }
      if (areAllGenDMatchingPassed && toCollectAllMatched &&
          areAllDJetMatched && arrIsPassingZMassCut[i]) {
        arrTTAllMatched[i]->Fill();
      }
    }
  };

  UInt_t nJet, nJetPassed;
  UInt_t nFatJet, nFatJetPassed;
  for (Byte_t i = 0; i < 2; i++) { // Originally [0, 3)
    TString nameJetPrefix = "nJet";
    TString titleJetPrefix = "Number of jets";
    TString nameFatJetPrefix = "nFatJet";
    TString titleFatJetPrefix = "Number of fat jets";
    TString namePassSuffix = "Passed";
    TString titlePassSuffix = " passing preselections";
    arrTTGen[i]
        ->Branch(nameJetPrefix, &nJet, nameJetPrefix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix);
    arrTTGen[i]
        ->Branch(nameJetPrefix + namePassSuffix, &nJetPassed, nameJetPrefix + namePassSuffix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix + titlePassSuffix);
    arrTTPreselectedMatching[i]
        ->Branch(nameJetPrefix + namePassSuffix, &nJetPassed, nameJetPrefix + namePassSuffix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix + titlePassSuffix);
    arrTTAllMatched[i]
        ->Branch(nameJetPrefix + namePassSuffix, &nJetPassed, nameJetPrefix + namePassSuffix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix + titlePassSuffix);
    arrTTGen[i]
        ->Branch(nameFatJetPrefix, &nFatJet, nameFatJetPrefix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix);
    arrTTGen[i]
        ->Branch(nameFatJetPrefix + namePassSuffix, &nFatJetPassed, nameFatJetPrefix + namePassSuffix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix + titlePassSuffix);
    arrTTPreselectedMatching[i]
        ->Branch(nameFatJetPrefix + namePassSuffix, &nFatJetPassed, nameFatJetPrefix + namePassSuffix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix + titlePassSuffix);
    arrTTAllMatched[i]
        ->Branch(nameFatJetPrefix + namePassSuffix, &nFatJetPassed, nameFatJetPrefix + namePassSuffix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix + titlePassSuffix);
  }

  std::vector<UInt_t> vIdxJetPassed;
  // std::vector<UInt_t> vRankJetPassedPt;
  {
    TString nameIdx = "Jet_idxPassed";
    TString titleIdx = "Index of jets passing preselections";
    TString nameRank = "Jet_rankJetPassedPt";
    TString titleRank = "Rank (0-indexed) of pt among passed jets";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTGen[i]->Branch(nameIdx, &vIdxJetPassed)->SetTitle(titleIdx);
      // arrTTGen[i]->Branch(nameRank, &vRankJetPassedPt)->SetTitle(titleRank);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTNumCorrect[i]->Branch(nameIdx, &vIdxJetPassed)->SetTitle(titleIdx);
      // arrTTNumCorrect[i]->Branch(nameRank,
      // &vRankJetPassedPt)->SetTitle(titleRank);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTZMassCutted[i]->Branch(nameIdx, &vIdxJetPassed)->SetTitle(titleIdx);
      // arrTTZMassCutted[i]->Branch(nameRank,
      // &vRankJetPassedPt)->SetTitle(titleRank);
    }
  }
  std::vector<Int_t> vDIdxJetClosest(4, -1);
  std::vector<Int_t> vDJetClosestRankJetPassedPt(4, -10);
  std::vector<UInt_t> vJetMatchedRankJetPassedPt;
  {
    TString nameIdx = "GenDMatching_idxJet";
    TString titleIdx = "Index of closest jets of each D quarks";
    TString nameRank = "GenDMatching_rankJetPassedPt";
    TString titleRank =
        "Rank (0-indexed) of pt of the closest jet of each D quark among "
        "passed jets";
    TString nameRankUnique = "JetMatched_rankJetPassedPt";
    TString titleRankUnique =
        "Rank (0-indexed) of pt of each matched jet among passed jets";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatching[i]
          ->Branch(nameIdx, &vDIdxJetClosest)
          ->SetTitle(titleIdx);
      arrTTPreselectedMatching[i]
          ->Branch(nameRank, &vDJetClosestRankJetPassedPt)
          ->SetTitle(titleRank);
      arrTTPreselectedMatching[i]
          ->Branch(nameRankUnique, &vJetMatchedRankJetPassedPt)
          ->SetTitle(titleRankUnique);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTAllMatched[i]->Branch(nameIdx, &vDIdxJetClosest)->SetTitle(titleIdx);
      arrTTAllMatched[i]
          ->Branch(nameRank, &vDJetClosestRankJetPassedPt)
          ->SetTitle(titleRank);
      arrTTAllMatched[i]
          ->Branch(nameRankUnique, &vJetMatchedRankJetPassedPt)
          ->SetTitle(titleRankUnique);
    }
  }
  Bool_t areAllJetsMatchedLeading;
  {
    TString name = "areAllJetsMatchedLeading";
    TString title =
        "Whether all the jets matched are leading jets (rankMax == "
        "nJetsMatched-1)";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatching[i]
          ->Branch(name, &areAllJetsMatchedLeading, name + "/O", 1)
          ->SetTitle(title);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTAllMatched[i]
          ->Branch(name, &areAllJetsMatchedLeading, name + "/O", 1)
          ->SetTitle(title);
    }
  }
  std::vector<Float_t> vDeltaRDJet(4, -1);
  std::vector<Float_t> vDeltaRJetPair(2, -1);
  std::vector<Float_t> vDeltaRBetweenJetPairs(1, -1);
  {
    TString name = "GenDMathcing_deltaRJet";
    TString title = "DeltaR between each d quark and its closest jet";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatching[i]->Branch(name, &vDeltaRDJet)->SetTitle(title);
    }
  }
  {
    TString name = "GenDMathcingPair_deltaRJet";
    TString title = "DeltaR of each jet pair most likely to match the d quarks";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatching[i]
          ->Branch(name, &vDeltaRJetPair)
          ->SetTitle(title);
    }
  }
  {
    TString name = "GenDMathcingPairPair_deltaRJet";
    TString title =
        "DeltaR between jet pairs most likely to match the d quarks";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatching[i]
          ->Branch(name, &vDeltaRBetweenJetPairs)
          ->SetTitle(title);
    }
  }

  DescriptionCollectorForUntuplizer collectorGenDMatching;
  DescriptionCollectorForUntuplizer collectorHT;
  DescriptionCollectorForUntuplizer collectorJet;
  DescriptionCollectorForUntuplizer collectorFatJet;
  UInt_t nLeafOriginal;
  {
    TFile* tfIn = TFile::Open(inputFile.data(), "READ");
    TTree* ttIn = (TTree*)tfIn->Get(nameTreeIn);
    TObjArray* tarrLeafOriginal = ttIn->GetListOfLeaves();

    nLeafOriginal = tarrLeafOriginal->GetSize();

    // std::vector<Bool_t> arrVIdxLeafLeptonID[2];
    std::vector<Int_t> vIdxLeafJet;
    vIdxLeafJet.clear();
    for (UInt_t i = 0; i < nLeafOriginal; i++) {
      TLeaf* tlfCurrent = (TLeaf*)tarrLeafOriginal->At(i);
      // ObjectDescription descriptionCurrent;
      BranchDescription descriptionCurrent;
      descriptionCurrent.name = (TString)tlfCurrent->GetName();
      descriptionCurrent.title = (TString)tlfCurrent->GetBranch()->GetTitle();
      // TString typeNameCurrent = tlfCurrent->GetTypeName();
      descriptionCurrent.typeName = tlfCurrent->GetTypeName();
      if (descriptionCurrent.name.BeginsWith("Jet_")) {
        collectorJet.BranchFor(descriptionCurrent);
      }
      if (descriptionCurrent.name.BeginsWith("FatJet_")) {
        collectorFatJet.BranchFor(descriptionCurrent);
      }
      if (descriptionCurrent.name.BeginsWith("SoftActivityJetHT")) {
        // if (typeNameCurrent == "Float_t") {
        if (descriptionCurrent.typeName.EqualTo("Float_t")) {
          if (debug)
            std::cout << "Found " << descriptionCurrent.name
                      << " with type "
                      // << typeNameCurrent << std::endl;
                      << descriptionCurrent.typeName << std::endl;
          // vLeafDescriptionHTFloat.push_back(descriptionCurrent);
          collectorHT.BranchFor(descriptionCurrent);
        }
      }
      if (descriptionCurrent.name.BeginsWith("GenPart_") &&
          !descriptionCurrent.name.EndsWith("_idx") &&
          !descriptionCurrent.name.EndsWith("_pt") &&
          !descriptionCurrent.name.EndsWith("_eta") &&
          !descriptionCurrent.name.EndsWith("_phi") &&
          !descriptionCurrent.name.EndsWith("_mass")) {
        collectorGenDMatching.BranchFor(descriptionCurrent);
      }
    }
  }

  collectorGenDMatching.Prepare();
  collectorHT.Prepare();
  collectorJet.Prepare();
  collectorFatJet.Prepare();
  BranchMounterIterableForUntuplizer
  mounterGenDMatching(collectorGenDMatching, data, [](BranchDescription
  descriptionCurrent){
    Ssiz_t lenName = descriptionCurrent.name.Length();
    const Ssiz_t lenPrefOriginal = 8;
    return BranchDescription(
        "GenDMatching_" +
            descriptionCurrent.name(lenPrefOriginal, lenName - 1),
        descriptionCurrent.title, descriptionCurrent.typeName);
  });
  // BranchMounterIterableForUntuplizer mounterGenDMatching(collectorGenDMatching, data);
  BranchMounterScalarForUntuplizer mounterHT(collectorHT, data);
  BranchMounterIterableForUntuplizer mounterJet(collectorJet, data);
  BranchMounterIterableForUntuplizer mounterFatJet(collectorFatJet, data);
  {
    auto funDo = [&mounterGenDMatching, &mounterHT, &mounterJet, &mounterFatJet](TTree* tree) {
      mounterGenDMatching.BranchOn(tree);
      mounterJet.BranchOn(tree);
      mounterHT.BranchOn(tree);
      mounterFatJet.BranchOn(tree);
    };
    std::for_each(std::begin(arrTTGen), std::end(arrTTGen), funDo);
    // std::for_each(std::begin(arrTTZMassCutted), std::end(arrTTZMassCutted),
    // mounter.BranchOn);
    std::for_each(std::begin(arrTTPreselectedMatching),
                  std::end(arrTTPreselectedMatching), funDo);
    std::for_each(std::begin(arrTTAllMatched), std::end(arrTTAllMatched),
                  funDo);
  }

  for (jEntry = 0; jEntry < nEntry; jEntry++) {
    data.GetEntry(jEntry);

    nJet = data.GetInt("nJet");
    nFatJet =  data.GetInt("nFatJet");
    // if (debug) std::cout << "nJet: " << nJet << std::endl;

    Bool_t hasGenLepton = false;
    UInt_t nGenPart = data.GetInt("nGenPart");
    Int_t* ptrGenPart_pdgId = data.GetPtrInt("GenPart_pdgId");
    Int_t* ptrGenPart_genPartIdxMother =
        data.GetPtrInt("GenPart_genPartIdxMother");
    for (Byte_t i = 0; i < 3; i++) arrHasGenLepton[i] = false;

    haveAllGenDMaching = true;
    Bool_t arrHasGenD[4];
    std::fill(std::begin(arrHasGenD), std::end(arrHasGenD), false);
    std::fill(std::begin(vGenDMatchingIdx), std::end(vGenDMatchingIdx), 0);
    for (UInt_t ig = 0; ig < nGenPart; ig++) {
      if (ptrGenPart_pdgId[ptrGenPart_genPartIdxMother[ig]] == pdgZ) {
        Int_t genparIdAbs = TMath::Abs(ptrGenPart_pdgId[ig]);
        if (genparIdAbs == pdgElectron) {
          arrHasGenLepton[0] = true;
          hasGenLepton = true;
          // break;
        } else if (genparIdAbs == pdgMuon) {
          arrHasGenLepton[1] = true;
          hasGenLepton = true;
          // break;
        } else if (genparIdAbs == pdgTau) {
          arrHasGenLepton[2] = true;
          hasGenLepton = true;
          // break;
        }
      }
      if (TMath::Abs(ptrGenPart_pdgId[ptrGenPart_genPartIdxMother[ig]]) ==
              pdgX2 &&
          TMath::Abs(ptrGenPart_pdgId[ig]) == pdgDown) {
        Byte_t iD =
            (signbit(ptrGenPart_pdgId[ptrGenPart_genPartIdxMother[ig]]) << 1) +
            signbit(ptrGenPart_pdgId[ig]);
        arrHasGenD[iD] = true;
        vGenDMatchingIdx[iD] = ig;
      }
      for (Byte_t i = 0; i < 4; i++) {
        if (!arrHasGenD[i]) {
          haveAllGenDMaching = false;
          break;
        }
      }
    }
    Bool_t isGenRightEvent = hasGenLepton && haveAllGenDMaching;

    if (BranchMounterHelper::RunForEachNIdxSorted<
            BranchMounterIterableForUntuplizer, std::vector<Int_t>::iterator,
            Int_t>(mounterGenDMatching, vGenDMatchingIdx.begin(), 4) < 4 &&
        debug) {
      Error("BranchMounterHelper::RunForEachIdxSorted",
            "Not all d indices found!");
    };

    Int_t nElectron = data.GetInt("nElectron");
    Bool_t* ptrElectron_mvaFall17V2Iso_WPL =
        data.GetPtrBool("Electron_mvaFall17V2Iso_WPL");
    Bool_t* ptrElectron_mvaFall17V2Iso_WP90 =
        data.GetPtrBool("Electron_mvaFall17V2Iso_WP90");
    Bool_t* ptrElectron_mvaFall17V2Iso_WP80 =
        data.GetPtrBool("Electron_mvaFall17V2Iso_WP80");
    Bool_t* ptrElectron_isSoft = ptrElectron_mvaFall17V2Iso_WP80;
    Bool_t* ptrElectron_isTight = ptrElectron_mvaFall17V2Iso_WP90;
    Int_t nMuon = data.GetInt("nMuon");
    Bool_t* ptrMuon_softId = data.GetPtrBool("Muon_softId");
    Bool_t* ptrMuon_tightId = data.GetPtrBool("Muon_tightId");
    Bool_t* ptrMuon_isSoft = ptrMuon_softId;
    Bool_t* ptrMuon_isTight = ptrMuon_tightId;

    Int_t nElectronSoft, nElectronTight;
    {
      Int_t i;
      for (i = 0, nElectronSoft = 0, nElectronTight = 0; i < nElectron;
           nElectronSoft += ptrElectron_isSoft[i],
          nElectronTight += ptrElectron_isTight[i], i++)
        ;
    }
    Bool_t isNumElectronCorrect = (nElectronSoft == 2 && nElectronTight == 2);
    // nElectronPair = isNumElectronCorrect;

    Int_t nMuonSoft, nMuonTight;
    {
      Int_t i;
      for (i = 0, nMuonSoft = 0, nMuonTight = 0; i < nElectron;
           nMuonSoft += ptrMuon_isSoft[i], nMuonTight += ptrMuon_isTight[i],
          i++)
        ;
    }
    Bool_t isNumMuonCorrect = (nMuonSoft == 2 && nMuonTight == 2);
    // nMuonPair = isNumMuonCorrect;

    Bool_t isNumCorrect = (isNumElectronCorrect != isNumMuonCorrect);
    arrRecoIsLeptonNumCorrect[0] = isNumCorrect && isNumElectronCorrect;
    arrRecoIsLeptonNumCorrect[1] = isNumCorrect && isNumMuonCorrect;

    // if (!isNumCorrect) {
    //   lambdaFillTheTrees();
    //   continue;
    // }
    if (debug)
      std::cout << "elepair, mupair: " << (Int_t)isNumElectronCorrect
                << (Int_t)isNumMuonCorrect << std::endl;

    for (Byte_t i = 0; i < 2; i++) arrIsPassingZMassCut[i] = 0;

    Float_t* ptrElectron_pt = data.GetPtrFloat("Electron_pt");
    Float_t* ptrElectron_phi = data.GetPtrFloat("Electron_phi");
    Float_t* ptrElectron_eta = data.GetPtrFloat("Electron_eta");
    Float_t* ptrElectron_eCorr = data.GetPtrFloat("Electron_eCorr");

    TLorentzVector* ppElectronP4NumberCorrect[nElectronTight];
    TLorentzVector* ptrElectronP4NumberCorrectSum;

    if (arrRecoIsLeptonNumCorrect[0]) {
      ptrElectronP4NumberCorrectSum = new TLorentzVector();
      for (Int_t i = 0; i < nElectronTight; i++) {
        ppElectronP4NumberCorrect[i] = new TLorentzVector();
        // ppElectronP4NumberCorrect[i]->SetPtEtaPhiE(
        //     ptrElectron_pt[i], ptrElectron_eta[i], ptrElectron_phi[i],
        //     ptrElectron_eCorr[i]);
        ppElectronP4NumberCorrect[i]->SetPtEtaPhiM(
            ptrElectron_pt[i], ptrElectron_eta[i], ptrElectron_phi[i],
            massElectron);
        *ptrElectronP4NumberCorrectSum += *ppElectronP4NumberCorrect[i];
      }
      // Double_t electronPair_mass = ptrElectronP4NumberCorrectSum->M();
      // vElectronPair_massNumCorrect.push_back(electronPair_mass);
      arrLeptonPair_mass[0] = ptrElectronP4NumberCorrectSum->M();

      // leptonPair_massNumCorrect = electronPair_mass;
      // Mass z cuts
      Double_t massZCutUpper = massZ + 20;
      Double_t massZCutLower = massZ - 20;
      // Z mass cut for electron pairs
      arrIsPassingZMassCut[0] = (massZCutLower < arrLeptonPair_mass[0] &&
                                 arrLeptonPair_mass[0] < massZCutUpper);
      if (debug)
        std::cout << "massLPair: " << arrLeptonPair_mass[0] << std::endl;
    }

    Float_t* ptrMuon_pt = data.GetPtrFloat("Muon_pt");
    Float_t* ptrMuon_phi = data.GetPtrFloat("Muon_phi");
    Float_t* ptrMuon_eta = data.GetPtrFloat("Muon_eta");
    Float_t* ptrMuon_mass = data.GetPtrFloat("Muon_mass");

    TLorentzVector* ppMuonP4NumberCorrect[nMuonTight];
    TLorentzVector* ptrMuonP4NumberCorrectSum;
    if (arrRecoIsLeptonNumCorrect[1]) {
      ptrMuonP4NumberCorrectSum = new TLorentzVector();
      for (Int_t i = 0; i < nMuonTight; i++) {
        ppMuonP4NumberCorrect[i] = new TLorentzVector();
        ppMuonP4NumberCorrect[i]->SetPtEtaPhiM(ptrMuon_pt[i], ptrMuon_eta[i],
                                               ptrMuon_phi[i], ptrMuon_mass[i]);
        *ptrMuonP4NumberCorrectSum += *ppMuonP4NumberCorrect[i];
      }
      // Double_t muonPair_mass = ptrMuonP4NumberCorrectSum->M();
      // vMuonPair_massNumCorrect.push_back(muonPair_mass);
      arrLeptonPair_mass[1] = ptrMuonP4NumberCorrectSum->M();

      // leptonPair_massNumCorrect = muonPair_mass;
      Double_t massZCutUpper = massZ + 20;
      Double_t massZCutLower = massZ - 20;
      // Z mass cut for muon pairs
      arrIsPassingZMassCut[1] = (massZCutLower < arrLeptonPair_mass[1] &&
                                 arrLeptonPair_mass[1] < massZCutUpper);
      if (debug)
        std::cout << "massLPair: " << arrLeptonPair_mass[1] << std::endl;
    }

    if (!arrIsPassingZMassCut[0] && !arrIsPassingZMassCut[1] &&
        !isGenRightEvent) {
      lambdaFillTheTrees();
      continue;  // TODO
    }

    // std::vector<Int_t> vIdxZ_all, vIdxZp_all;
    // vIdxZ_all.clear();
    // vIdxZp_all.clear();
    // for (Int_t i=0; i<nGenPart; i++) {
    //   switch (ptrGenPart_pdgId[i]) {
    //   case pdgZ:
    //     vIdxZ_all.push_back(i);
    //   case pdgZp:
    //     vIdxZp_all.push_back(i);
    //   }
    // }

    // ttNumCorrect->Fill();

    // Jet_*
    vIdxJetPassed.clear();
    Float_t* ptrJet_pt_original = data.GetPtrFloat("Jet_pt");
    Float_t* ptrJet_eta_original = data.GetPtrFloat("Jet_eta");
    Float_t* ptrJet_phi_original = data.GetPtrFloat("Jet_phi");
    Float_t* ptrJet_mass_original = data.GetPtrFloat("Jet_mass");
    // for (UInt_t idxJet = 0; idxJet < nJet; idxJet++) {
    //   Bool_t isPassed = false;
    //   if (ptrJet_pt_original[idxJet] < 30) continue;
    //   if (TMath::Abs(ptrJet_eta_original[idxJet]) > 3) continue;
    //   vIdxJetPassed.push_back(idxJet);
    // }
    mounterJet.PrepareForPushing(nJet);
    {
      auto ptrJet_pt_original_cloned = ptrJet_pt_original;
      auto ptrJet_eta_original_cloned = ptrJet_eta_original;
      // auto ptrJet_phi_original_cloned = ptrJet_phi_original;
      // auto ptrJet_mass_original_cloned = ptrJet_mass_original;
      for (UInt_t idxJet = 0; idxJet < nJet; idxJet++) {
        Bool_t result = *ptrJet_pt_original_cloned > 30 &&
                        TMath::Abs(*ptrJet_eta_original_cloned) < 3;
        if (result) {
          vIdxJetPassed.push_back(idxJet);
        }
        mounterJet.PushOrSkip(result);
        ptrJet_pt_original_cloned++;
        ptrJet_eta_original_cloned++;
        // ptrJet_phi_original_cloned++;
        // ptrJet_mass_original_cloned++;
      }
    }
    nJetPassed = vIdxJetPassed.size();
    mounterFatJet.PrepareForPushing(nFatJet);
    {
      for (UInt_t idxFatJet = 0; idxFatJet < nFatJet; idxFatJet++) {
        mounterFatJet.PushOrSkip(true);
      }
    }
    mounterHT.Push();

    Byte_t iLeptonFound =
        arrIsPassingZMassCut[0] ? 0 : 1;  // Electron event or muon event

    Float_t* ptrGenPart_pt = data.GetPtrFloat("GenPart_pt");
    Float_t* ptrGenPart_eta = data.GetPtrFloat("GenPart_eta");
    Float_t* ptrGenPart_phi = data.GetPtrFloat("GenPart_phi");
    std::fill(std::begin(vDIdxJetClosest), std::end(vDIdxJetClosest), -10);
    std::fill(std::begin(vDJetClosestRankJetPassedPt),
              std::end(vDJetClosestRankJetPassedPt), -10);
    std::vector<UInt_t> vDJetMatchedRankJetPassedPtActual(4);
    vDJetMatchedRankJetPassedPtActual.clear();
    TLorentzVector* p4DMatching[4];
    TLorentzVector* p4DJetClosest[4];
    TLorentzVector* p4DPairMatching[2];
    nGenDMatchingPassed = 0;
    for (Byte_t iDMatching = 0; iDMatching < 4; iDMatching++) {
      if (debug) std::printf("iDMatching: %d\n", iDMatching);
      p4DJetClosest[iDMatching] = nullptr;
      vDeltaRDJet[iDMatching] = __FLT_MAX__;
      p4DMatching[iDMatching] = new TLorentzVector;
      vGenDMatchingPt[iDMatching] =
          ptrGenPart_pt[vGenDMatchingIdx[iDMatching]];
      vGenDMatchingEta[iDMatching] =
          ptrGenPart_eta[vGenDMatchingIdx[iDMatching]];
      vGenDMatchingPhi[iDMatching] =
          ptrGenPart_phi[vGenDMatchingIdx[iDMatching]];
      vGenDMatchingIdPassed[iDMatching] =
          vGenDMatchingPt[iDMatching] > 30 && TMath::Abs(vGenDMatchingEta[iDMatching]) < 3;
      if (vGenDMatchingIdPassed[iDMatching]) {
        nGenDMatchingPassed++;
      }
      p4DMatching[iDMatching]->SetPtEtaPhiM(
          vGenDMatchingPt[iDMatching], vGenDMatchingEta[iDMatching], vGenDMatchingPhi[iDMatching],
          massDown);
      for (Byte_t rankJetPassed = 0; rankJetPassed < nJetPassed;
           rankJetPassed++) {
        if (debug) std::printf("rankJetPassed: %d\n", rankJetPassed);
        TLorentzVector* p4Jet = new TLorentzVector;
        p4Jet->SetPtEtaPhiM(ptrJet_pt_original[vIdxJetPassed[rankJetPassed]],
                            ptrJet_eta_original[vIdxJetPassed[rankJetPassed]],
                            ptrJet_phi_original[vIdxJetPassed[rankJetPassed]],
                            ptrJet_mass_original[vIdxJetPassed[rankJetPassed]]);
        Double_t deltaRDJet = p4DMatching[iDMatching]->DeltaR(*p4Jet);
        if (deltaRDJet < vDeltaRDJet[iDMatching]) {
          // p4DJetClosest[iDMatching]->Delete();
          // delete p4DJetClosest[iDMatching];
          p4DJetClosest[iDMatching] = p4Jet;
          vDeltaRDJet[iDMatching] = deltaRDJet;
          vDJetClosestRankJetPassedPt[iDMatching] = rankJetPassed;
          vDIdxJetClosest[iDMatching] = vIdxJetPassed[rankJetPassed];
        } else {
          // p4Jet->Delete();
          // delete p4Jet;
        }
      }
      if (vDeltaRDJet[iDMatching] < 0.4) {
        vDJetMatchedId[iDMatching] = true;
        vDJetMatchedRankJetPassedPtActual.push_back(
            vDJetClosestRankJetPassedPt[iDMatching]);
      } else {
        vDJetMatchedId[iDMatching] = false;
      }
    }
    // std::cout << "testtag" << std::endl;
    // for (UInt_t iDPairMatching=0; iDPairMatching<2; iDPairMatching++) {
    //   *(p4DPairMatching[iDPairMatching]) = *(p4DMatching[iDPairMatching*2]) + *(p4DMatching[(iDPairMatching*2)+1]);
    //   vGenDPairMatchingPt[iDPairMatching] = p4DPairMatching[iDPairMatching]->Pt();
    //   vGenDPairMatchingEta[iDPairMatching] = p4DPairMatching[iDPairMatching]->Eta();
    //   vGenDPairMatchingPhi[iDPairMatching] = p4DPairMatching[iDPairMatching]->Phi();
    // }
    areAllGenDMatchingPassed = nGenDMatchingPassed == 4;
    if (!areAllGenDMatchingPassed) {
      lambdaFillTheTrees();
      continue;
    }
    for (Byte_t i = 0; i < 2; i++) {
      vDeltaRJetPair[i] =
          p4DJetClosest[i << 1]->DeltaR(*p4DJetClosest[(i << 1) + 1]);
    }
    vDeltaRBetweenJetPairs[0] =
        (*p4DJetClosest[0] + *p4DJetClosest[1])
            .DeltaR(*p4DJetClosest[2] + *p4DJetClosest[3]);

    nDJetMatched = vDJetMatchedRankJetPassedPtActual.size();
    if (nDJetMatched) {
      std::sort(vDJetMatchedRankJetPassedPtActual.begin(),
                vDJetMatchedRankJetPassedPtActual.end());
      // UInt_t rankMax = vDJetMatchedRankJetPassedPtActual.back();
      typename std::vector<UInt_t>::iterator iterJetMatchedRankJetPassedPtEnd =
          std::unique(vDJetMatchedRankJetPassedPtActual.begin(),
                      vDJetMatchedRankJetPassedPtActual.end());
      nJetMatched = std::distance(vDJetMatchedRankJetPassedPtActual.begin(),
                                  iterJetMatchedRankJetPassedPtEnd);
      vJetMatchedRankJetPassedPt.clear();
      vJetMatchedRankJetPassedPt.assign(
          vDJetMatchedRankJetPassedPtActual.begin(),
          iterJetMatchedRankJetPassedPtEnd);
      // areAllJetsMatchedLeading = rankMax == nJetMatched - 1;
      areAllJetsMatchedLeading =
          vJetMatchedRankJetPassedPt.back() == nJetMatched - 1;
    } else {
      nJetMatched = 0;
      areAllJetsMatchedLeading = true;
    }

    areAllDJetMatched = (nDJetMatched == 4);
    if (!areAllDJetMatched) {
      lambdaFillTheTrees();
      continue;
    }

    lambdaFillTheTrees();
    // ttZMassCutted->Fill();
  }

  // ttNumCorrect->GetCurrentFile()->Write();
  // ttZMassCutted->GetCurrentFile()->Write();
  outFileTree->Write();

  // // TString outputFileHist = (TString) "output_" + outputFileHead + "_" +
  // //                          outputFileVar + "_" + outputFileTail + "_hist"
  // + ".root";

  // // TFile *outFileHist = new TFile(outputFileHist.Data(), toRecreateOutFile
  // ? "recreate" : "update");

  // TString outImageDir =
  //     (TString) "../out_images/output_" + outputFileHead + "_" +
  //     outputFileTail;

  // // TCanvas* c1 = new TCanvas;
  // gStyle->SetOptStat(111111);
  // auto lambdaPrintHistograms = [/*&c1,*/ /*&outFileHist,*/ toRecreateOutFile,
  // outImageDir, outputFileHead, outputFileVar, outputFileTail](TTree *tt)->
  // void {
  //   std::cout << "Printing " << tt->GetName() << std::endl;
  //   TString nameTree = tt->GetName();
  //   TString outImageNameHead = (TString)outputFileHead + "_" + outputFileVar
  //   + "_" + nameTree; TString outImageCommonPath = outImageDir + "/" +
  //   outImageNameHead + "_"; TString outputFileHist = (TString) "output_" +
  //   outputFileHead + "_" +
  //                          outputFileVar + "_" + outputFileTail + "_hist_" +
  //                          nameTree + ".root";
  //   TFile *outFileHist = new TFile(outputFileHist.Data(), toRecreateOutFile ?
  //   "recreate" : "update"); for (TObject* leafObject :
  //   *(tt->GetListOfLeaves())) {
  //     TLeaf* leaf = (TLeaf*)leafObject;
  //     TString nameLeaf = leaf->GetName();
  //     if (nameLeaf.First("[") >= 0) {
  //       nameLeaf.Resize(nameLeaf.First("["));
  //     }
  //     TString nameHist = "h" + nameLeaf;
  //     tt->Draw(nameLeaf + ">>" + nameHist);
  //     TH1 *hist = (TH1 *) gDirectory->Get(nameHist);
  //     hist->SetTitle((TString)leaf->GetBranch()->GetTitle() + " (" + nameTree
  //     + ")"); TString outImagePath = outImageCommonPath + nameLeaf + ".svg";
  //     // c1->Clear();
  //     // hist->Draw();
  //     // c1->Print(outImagePath);
  //     outFileHist->WriteObject(hist, nameHist);
  //   }
  //   outFileHist->Close();
  // };
  // if (true) {
  // for (Byte_t i=0; i<3; i++) {
  //   lambdaPrintHistograms(arrTTGen[i]);
  // }
  // }
  // // for (Byte_t i=0; i<2; i++) {
  // //   lambdaPrintHistograms(arrTTNumCorrect[i]);
  // // }
  // for (Byte_t i=0; i<2; i++) {
  //   lambdaPrintHistograms(arrTTZMassCutted[i]);
  // }
  // // c1->Close();

  outFileTree->Close();
  // delete outFileTree;
  // for (UInt_t i = 0; i < nLeafJetFloat; i++) {
  //   // arrHLeafJetFloat[i]->Write();
  //   // outFileHist->WriteObject(arrHLeafJetFloat[i],
  //   // arrHLeafJetFloat[i]->GetName());
  // }
  // for (UInt_t i = 0; i < nLeafJetUInt; i++) {
  //   // arrHLeafJetUInt[i]->Write();
  //   // outFileHist->WriteObject(arrHLeafJetUInt[i],
  //   // arrHLeafJetUInt[i]->GetName());
  // }
  // for (UInt_t i = 0; i < nLeafJetInt; i++) {
  //   // arrHLeafJetInt[i]->Write();
  //   // outFileHist->WriteObject(arrHLeafJetInt[i],
  //   // arrHLeafJetInt[i]->GetName());
  // }
  // for (UInt_t i = 0; i < nLeafJetBool; i++) {
  //   // arrHLeafJetBool[i]->Write();
  //   // outFileHist->WriteObject(arrHLeafJetBool[i],
  //   // arrHLeafJetBool[i]->GetName());
  // }

  // // outFileHist->Close();
}
