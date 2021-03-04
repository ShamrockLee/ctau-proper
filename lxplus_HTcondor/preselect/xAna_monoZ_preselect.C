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
#include <iterator>

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

#ifndef ARRAY_SIZE_FUNCTION
#define ARRAY_SIZE_FUNCTION

template <typename T, size_t n>
size_t array_size(T const (&)[n]) {
  return n;
}

#endif

#ifndef IOTA_ITERATOR_CLASS
#define IOTA_ITERATOR_CLASS
template<typename TypeValue = Int_t>
class IotaIterator {
  typedef std::ptrdiff_t TypeDifference;
 protected:
  TypeValue fVal;
 public:
  // typedef std::random_access_iterator_tag iterator_category;
  typedef TypeValue value_type;
  typedef TypeDifference difference_type;
  typedef TypeValue* pointer_type;
  typedef TypeValue& reference_type;
  IotaIterator<TypeValue>(const TypeValue val=0)
  : fVal(val) {
  }
  const TypeValue& operator*() const {
    // std::cout << "Iota iterator indirects " << fVal << std::endl;
    return fVal;
  }
  IotaIterator<TypeValue>& operator+=(const TypeDifference n) {
    fVal += n;
    return *this;
  }
  IotaIterator<TypeValue>& operator-=(const TypeDifference n) {
    fVal -= n;
    return *this;
  }
  IotaIterator<TypeValue>& operator++(){
    ++fVal;
    return *this;
  }
  IotaIterator<TypeValue>& operator--(){
    --fVal;
    return *this;
  }
  IotaIterator<TypeValue> operator+(const TypeDifference n) const{
    IotaIterator<TypeValue> result = *this;
    result += n;
    return result;
  }
  IotaIterator<TypeValue> operator-(const TypeDifference n) const {
    IotaIterator<TypeValue> result = *this;
    result -= n;
    return result;
  }
  IotaIterator<TypeValue> operator++(int) {
    IotaIterator<TypeValue> result = *this;
    ++*this;
    return result;
  }
  IotaIterator<TypeValue> operator--(int) {
    IotaIterator<TypeValue> result = *this;
    --*this;
    return result;
  }
};
#endif

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
      0.0005109989461;  //< static mass of electron (constant)
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

  Bool_t isSignal = outputFileTree.Contains("ignal");

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
  TTree* arrTTPreselectedMatchingJet[2];

  /// Tree for the events whose four d-quarks are all matched to a jet.
  TTree* arrTTAllMatchedJet[2];

  /// Tree for the fat jet matchisg results of events passing the preselections
  TTree* arrTTPreselectedMatchingFatJet[2];

  /// Tree for the events whose two d-pairs are all matched to a fat jet.
  TTree* arrTTAllMatchedFatJet[2];
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
    arrTTPreselectedMatchingJet[i] =
        new TTree((TString) "PreselectedMatchingJet" + namesLepton[i],
                  (TString) "Jet matching results of preselected " +
                      namesLeptonLower[i] + " events");
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTAllMatchedJet[i] =
        new TTree((TString) "AllMatchedJet" + namesLepton[i],
                  (TString) "Jet matching results of the all-matched " +
                      namesLeptonLower[i] + " events");
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTPreselectedMatchingFatJet[i] =
        new TTree((TString) "PreselectedMatchingFatJet" + namesLepton[i],
                  (TString) "Fat jet matching results of preselected " +
                      namesLeptonLower[i] + " events");
  }
  for (Byte_t i = 0; i < 2; i++) {
    arrTTAllMatchedFatJet[i] =
        new TTree((TString) "AllMatchedFatJet" + namesLepton[i],
                  (TString) "Fat jet matching results of the all-matched " +
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


  std::vector<Int_t> arrVGenLeptonIdx[3];
  for (Byte_t i = 0; i < 3; i++) {
    const TString name = "Gen" + namesLepton[i] + "_genPartIdx";
    const TString title = "Indexes of the d quarks";
    arrVGenLeptonIdx[i].clear();
    arrVGenLeptonIdx[i].resize(2);
    arrTTGen[i]->Branch(name, &(arrVGenLeptonIdx[i]))->SetTitle(title);
  }

  // std::vector<Int_t> arrVGenLeptonPairPtEtaPhiMass[4][3];
  // {
  //   auto arrCstrPtEtaPhiMassLower = {"Pt", "Eta", "Phi", "Mass"};
  //   for (Byte_t j = 0; j < 4; j++) {
  //     for (Byte_t i = 0; i < 3; i++) {
  //       {

  //       }
  //       {
  //         const TString name = "Gen" + namesLepton[i] + "_genPartIdx";
  //         const TString title = "Indexes of the GEN-level " + namesLeptonLower[i];
  //         arrVGenLeptonPairPtEtaPhiMass[i][j].clear();
  //         arrVGenLeptonPairPtEtaPhiMass[i][j].resize(1);
  //         arrTTGen[i]->Branch(name, &(arrVGenLeptonPairPtEtaPhiMass[i][j]))->SetTitle(title);
  //       }
  //     }
  //   }
  // }


  // Float_t* ptrElectron_pt = data.GetPtrFloat("Electron_pt");
  // Float_t* ptrElectron_phi = data.GetPtrFloat("Electron_phi");
  // Float_t* ptrElectron_eta = data.GetPtrFloat("Electron_eta");
  // Float_t* ptrElectron_eCorr = data.GetPtrFloat("Electron_eCorr");

  // auto collectorElectronDynamics = DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableFloat);
  // collectorElectronDynamics.BranchFor(BranchDescription("Electron_pt", "Pt of each electron", "Float_t"));
  // collectorElectronDynamics.BranchFor(BranchDescription("Electron_eta", "Eta of each electron", "Float_t"));
  // collectorElectronDynamics.BranchFor(BranchDescription("Electron_phi", "Phi of each electron", "Float_t"));
  // collectorElectronDynamics.BranchFor(BranchDescription("Electron_eCorr", "Corrected energy of each electron", "Float_t"));
  // collectorElectronDynamics.Prepare();

  // std::vector<DescriptionCollectorFromFun> vCollectorLeptonDynamics(2, DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableFloat));

  std::vector<BranchDescription> arrVDescriptionLeptonDynamics[2];
  for (Byte_t i=0; i<2; i++) {
    arrVDescriptionLeptonDynamics[i].clear();
    arrVDescriptionLeptonDynamics[i].reserve(4);
    for (TString&& tstrParameter: {"pt", "eta", "phi", "mass"}) {
      // vCollectorLeptonDynamics[i].BranchFor(BranchDescription(
      arrVDescriptionLeptonDynamics[i].push_back(BranchDescription(
          namesLepton[i] + "_" + tstrParameter,
          tstrParameter + " of each " + namesLeptonLower[i],
          static_cast<TString>("Float_t")));
    }
    // arrCollectorLeptonDynamics[i].Prepare();
  }
  arrVDescriptionLeptonDynamics[0].resize(arrVDescriptionLeptonDynamics[0].size()-1);

  // auto mounterElectronDynamics = BranchMounterVectorSingle<Float_t, Float_t *>(collectorElectronDynamics);
  // auto &vElectronPt = mounterElectronDynamics.vvE[0];
  // auto &vElectronEta = mounterElectronDynamics.vvE[1];
  // auto &vElectronPhi = mounterElectronDynamics.vvE[2];
  // auto &vElectronECorr = mounterElectronDynamics.vvE[3];

  std::vector<BranchMounterVectorSingle<Float_t, Float_t*>> vMounterLeptonDynamics = {};
  vMounterLeptonDynamics.reserve(2);
  std::vector<BranchMounterVectorSingle<Int_t, IotaIterator<Int_t> >> vMounterLeptonIdxPassedPtEta = {};
  vMounterLeptonIdxPassedPtEta.reserve(2);
  std::vector<BranchMounterIterableChained<>> vMounterLeptonChained = {};
  vMounterLeptonChained.reserve(2);
  {
    auto funIterEInPtrFloat = [&data, debug](BranchDescription description)->Float_t*{
      if (debug) std::cout << "Getting iterator for " << description.name << std::endl;
      return data.GetPtrFloat(description.name);
    };
    auto funIterEInIota = [](BranchDescription description)->IotaIterator<Int_t>{
      return IotaIterator<Int_t>(0);
    };
    for (Byte_t i=0; i<2; i++) {
      vMounterLeptonDynamics.push_back(BranchMounterVectorSingle<Float_t, Float_t*>(arrVDescriptionLeptonDynamics[i]));
      vMounterLeptonDynamics.back().SetFunIterEIn(funIterEInPtrFloat);
      vMounterLeptonIdxPassedPtEta.push_back(BranchMounterVectorSingle<Int_t, IotaIterator<Int_t>>({BranchDescription(namesLepton[i] + "_idxPassedPtEta", "Index of " + namesLepton[i] + " passing Pt and Eta cut", TString("Int_t"))}));
      vMounterLeptonIdxPassedPtEta.back().SetFunIterEIn(funIterEInIota);
      vMounterLeptonChained.push_back(BranchMounterIterableChained<>({&vMounterLeptonDynamics.back(), &vMounterLeptonIdxPassedPtEta.back()}));
    }
  }

  // mounterElectronDynamics.BranchOn(arrTTGen[0]);
  // mounterElectronDynamics.BranchOn(arrTTNumCorrect[0]);
  // mounterElectronDynamics.BranchOn(arrTTZMassCutted[0]);
  // mounterElectronDynamics.SetFunIterEIn([&data](BranchDescription description)->Float_t*{
  //   return data.GetPtrFloat(description.name);
  // });
  {
    for (Byte_t i=0; i<2; i++) {
      for (auto&& ppTT: {arrTTGen, arrTTNumCorrect, arrTTZMassCutted}) {
        vMounterLeptonChained[i].BranchOn(ppTT[i]);
      }
    }
  }

  std::vector<Int_t> arrVLeptonIdxNumCorrect[2];
  for (Byte_t i = 0; i < 2; i++) {
    arrVLeptonIdxNumCorrect[i].resize(2);
    for (TTree *const *const ppTT: {arrTTGen, arrTTNumCorrect, arrTTZMassCutted, arrTTPreselectedMatchingJet, arrTTPreselectedMatchingFatJet, arrTTAllMatchedJet, arrTTAllMatchedFatJet}) {
      ppTT[i]->Branch(namesLepton[i] + "_idxNumCorrect", &(arrVLeptonIdxNumCorrect[i]));
    }
  }

  Bool_t haveAllGenDMatching;
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "haveAllGenDMatching";
    TString title = "Whether there are 4 d quarks generated by chi1";
    arrTTGen[i]
        ->Branch(name, &haveAllGenDMatching, name + "/O", 1)
        ->SetTitle(title);
  }
  std::vector<Int_t> vGenDMatchingIdx(4, 0);
  for (Byte_t i = 0; i < 2; i++) {
    const TString name = "GenDMatching_genPartIdx";
    const TString title = "Indexes of the d quarks";
    arrTTGen[i]->Branch(name, &vGenDMatchingIdx)->SetTitle(title);
  }
  std::vector<Float_t> vGenDMatchingPt(4, 0);
  std::vector<Float_t> vGenDMatchingEta(4, 0);
  std::vector<Float_t> vGenDMatchingPhi(4, 0);
  std::vector<Float_t>* arrPVGenDMatchingP4Component[3] = {
      &vGenDMatchingPt, &vGenDMatchingEta, &vGenDMatchingPhi};
  std::vector<Float_t> vGenDPairMatchingPt(2, 0);
  std::vector<Float_t> vGenDPairMatchingEta(2, 0);
  std::vector<Float_t> vGenDPairMatchingPhi(2, 0);
  std::vector<Float_t> vGenDPairMatchingMass(2, 0);
  std::vector<Float_t>* arrPVGenDPairMatchingP4Component[4] = {
      &vGenDPairMatchingPt, &vGenDPairMatchingEta, &vGenDPairMatchingPhi,
      &vGenDPairMatchingMass};
  {
    const char* nameTemplate = "GenD%sMatching_%s";
    const char* titleTemplate = "%s of each %s";
    constexpr const char* const arrCstrPtEtaPhi[] = {"Pt", "Eta", "Phi"};
    for (Byte_t i = 0; i < 2; i++) {
      for (Byte_t j = 0; j < 3; j++) {
        TString tstrVarLower(arrCstrPtEtaPhi[j]);
        tstrVarLower.ToLower();
        arrTTGen[i]
            ->Branch(TString::Format(nameTemplate, "", tstrVarLower.Data()),
                     arrPVGenDMatchingP4Component[j])
            ->SetTitle(
                TString::Format(titleTemplate, arrCstrPtEtaPhi[j], "d quark"));
        arrTTPreselectedMatchingJet[i]
            ->Branch(TString::Format(nameTemplate, "", tstrVarLower.Data()),
                     arrPVGenDMatchingP4Component[j])
            ->SetTitle(
                TString::Format(titleTemplate, arrCstrPtEtaPhi[i], "d quark"));
        arrTTAllMatchedJet[i]
            ->Branch(TString::Format(nameTemplate, "", tstrVarLower.Data()),
                     arrPVGenDMatchingP4Component[j])
            ->SetTitle(
                TString::Format(titleTemplate, arrCstrPtEtaPhi[i], "d quark"));
        arrTTGen[i]
            ->Branch(TString::Format(nameTemplate, "Pair", tstrVarLower.Data()),
                     arrPVGenDPairMatchingP4Component[j])
            ->SetTitle(
                TString::Format(titleTemplate, arrCstrPtEtaPhi[i], "d-pair"));
        arrTTPreselectedMatchingFatJet[i]
            ->Branch(TString::Format(nameTemplate, "Pair", tstrVarLower.Data()),
                     arrPVGenDPairMatchingP4Component[j])
            ->SetTitle(
                TString::Format(titleTemplate, arrCstrPtEtaPhi[i], "d-pair"));
        arrTTAllMatchedFatJet[i]
            ->Branch(TString::Format(nameTemplate, "Pair", tstrVarLower.Data()),
                     arrPVGenDPairMatchingP4Component[j])
            ->SetTitle(
                TString::Format(titleTemplate, arrCstrPtEtaPhi[i], "d-pair"));
      }
      arrTTGen[i]
          ->Branch(TString::Format(nameTemplate, "Pair", "mass"),
                   &vGenDPairMatchingMass)
          ->SetTitle(TString::Format(titleTemplate, "Mass", "d-pair"));
      arrTTAllMatchedFatJet[i]
          ->Branch(TString::Format(nameTemplate, "Pair", "mass"),
                   &vGenDPairMatchingMass)
          ->SetTitle(TString::Format(titleTemplate, "Mass", "d-pair"));
    }
  }
  std::vector<Bool_t> vGenDMatchingIdPassed(4, false);
  UInt_t nGenDMatchingPassed = 0;
  Bool_t areAllGenDMatchingPassed = false;
  std::vector<Bool_t> vGenDPairMatchingIdPassed(2, false);
  UInt_t nGenDPairMatchingPassed = 0;
  Bool_t areAllGenDPairMatchingPassed;
  for (TTree *const *const ppTT: {arrTTGen, arrTTNumCorrect, arrTTZMassCutted}) {
  for (Byte_t i = 0; i < 2; i++) {
    const char* nameIdTemplate = "GenD%sMatching_idPassed";
    const char* titleIdTemplate =
        "Whether each of the %s pass the preselections";
    const char* nameNTemplate = "nGenD%sMatchingPassed";
    const char* titleNTemplate = "Number of %ss passing the preselections";
    const char* nameAreAllTemplate = "areAllGenD%sMatchingPassed";
    const char* titleAreAllTemplate =
        "Whether all the %ss have passed the preselections";
    ppTT[i]
        ->Branch(TString::Format(nameIdTemplate, ""), &vGenDMatchingIdPassed)
        ->SetTitle(TString::Format(titleIdTemplate, "d quark"));
    ppTT[i]
        ->Branch(TString::Format(nameNTemplate, ""), &nGenDMatchingPassed,
                 TString::Format(nameNTemplate, "") + "/I", __SIZEOF_INT__)
        ->SetTitle(TString::Format(titleNTemplate, "d quark"));
    ppTT[i]
        ->Branch(TString::Format(nameAreAllTemplate, ""),
                 &areAllGenDMatchingPassed,
                 TString::Format(nameAreAllTemplate, "") + "/O", 1)
        ->SetTitle(TString::Format(titleAreAllTemplate, "d quark"));
    ppTT[i]
        ->Branch(TString::Format(nameIdTemplate, "Pair"),
                 &vGenDPairMatchingIdPassed)
        ->SetTitle(TString::Format(titleIdTemplate, "d-pair"));
    ppTT[i]
        ->Branch(TString::Format(nameNTemplate, "Pair"),
                 &nGenDPairMatchingPassed,
                 TString::Format(nameNTemplate, "Pair") + "/I", __SIZEOF_INT__)
        ->SetTitle(TString::Format(titleNTemplate, "d-pair"));
    ppTT[i]
        ->Branch(TString::Format(nameAreAllTemplate, "Pair"),
                 &areAllGenDPairMatchingPassed,
                 TString::Format(nameAreAllTemplate, "Pair") + "/O", 1)
        ->SetTitle(TString::Format(titleAreAllTemplate, "d-pair"));
  }
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
  Float_t arrLeptonPair_pt[2];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = namesLepton[i] + "Pair_pt";
    TString title = "Pt of the " + namesLeptonLower[i] + " pair";
    for (TTree **ppTT: {arrTTNumCorrect, arrTTZMassCutted, arrTTPreselectedMatchingJet, arrTTPreselectedMatchingFatJet, arrTTAllMatchedJet, arrTTAllMatchedFatJet}) {
      ppTT[i]
          ->Branch(name, &(arrLeptonPair_pt[i]), name + "/F", __SIZEOF_FLOAT__)
          ->SetTitle(title);
    }
  }
  Float_t arrLeptonPair_eta[2];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = namesLepton[i] + "Pair_eta";
    TString title = "Eta of the " + namesLeptonLower[i] + " pair";
    for (TTree **ppTT: {arrTTNumCorrect, arrTTZMassCutted, arrTTPreselectedMatchingJet, arrTTPreselectedMatchingFatJet, arrTTAllMatchedJet, arrTTAllMatchedFatJet}) {
      ppTT[i]
          ->Branch(name, &(arrLeptonPair_eta[i]), name + "/F", __SIZEOF_FLOAT__)
          ->SetTitle(title);
    }
  }
  Float_t arrLeptonPair_phi[2];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = namesLepton[i] + "Pair_phi";
    TString title = "Phi of the " + namesLeptonLower[i] + " pair";
    for (TTree **ppTT: {arrTTNumCorrect, arrTTZMassCutted, arrTTPreselectedMatchingJet, arrTTPreselectedMatchingFatJet, arrTTAllMatchedJet, arrTTAllMatchedFatJet}) {
      ppTT[i]
          ->Branch(name, &(arrLeptonPair_phi[i]), name + "/F", __SIZEOF_FLOAT__)
          ->SetTitle(title);
    }
  }
  Float_t arrLeptonPair_mass[2];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = namesLepton[i] + "Pair_mass";
    TString title = "Mass of the " + namesLeptonLower[i] + " pair";
    for (TTree **ppTT: {arrTTNumCorrect, arrTTZMassCutted, arrTTPreselectedMatchingJet, arrTTPreselectedMatchingFatJet, arrTTAllMatchedJet, arrTTAllMatchedFatJet}) {
      ppTT[i]
          ->Branch(name, &(arrLeptonPair_mass[i]), name + "/F", __SIZEOF_FLOAT__)
          ->SetTitle(title);
    }
  }
  Bool_t arrIsPassingZMassCut[2];
  for (Byte_t i = 0; i < 2; i++) {
    TString name =
        TString::Format("is%sPairPassingZMassCut", namesLepton[i].Data());
    TString title = "Whether is the mass of the " + namesLeptonLower[i] +
                    " pair close to Z\'s";
    arrTTGen[i]
        ->Branch(name, &(arrIsPassingZMassCut[i]), name + "/O", 1)
        ->SetTitle(title);
    arrTTNumCorrect[i]
        ->Branch(name, &(arrIsPassingZMassCut[i]), name + "/O", 1)
        ->SetTitle(title);
  }

  constexpr const UInt_t nDQuarksExpected = 4;

  std::vector<Bool_t> vGenDMatchedId(nDQuarksExpected, false);
  std::vector<Bool_t> vGenDPairMatchedId(nDQuarksExpected >> 1, false);
  for (Byte_t i = 0; i < 2; i++) {
    arrTTPreselectedMatchingJet[i]
        ->Branch("GenDMatching_idMatched", &vGenDMatchedId)
        ->SetTitle("Whether each d-quarks matches a jet");
    arrTTPreselectedMatchingFatJet[i]
        ->Branch("GenDPairMatching_idMatched", &vGenDPairMatchedId)
        ->SetTitle("Whether each d-pairs matches a fat jet");
  }
  UInt_t nGenDMatched;
  UInt_t nGenDPairMatched;
  for (Byte_t i = 0; i < 2; i++) {
    TString nameTemplate = "nGenD%sMatched";
    TString titleTemplate = "Number of %ss matching a %s";
    arrTTPreselectedMatchingJet[i]
        ->Branch(TString::Format(nameTemplate, ""), &nGenDMatched,
                 TString::Format(nameTemplate, "") + "/I", __SIZEOF_INT__)
        ->SetTitle(TString::Format(titleTemplate, "d-quark", "jet"));
    arrTTPreselectedMatchingFatJet[i]
        ->Branch(TString::Format(nameTemplate, "Pair"), &nGenDPairMatched,
                 TString::Format(nameTemplate, "Pair") + "/I", __SIZEOF_INT__)
        ->SetTitle(TString::Format(titleTemplate, "d-pair", "fat jet"));
  }
  Bool_t areAllGenDMatched;
  Bool_t areAllGenDPairMatched;
  for (Byte_t i = 0; i < 2; i++) {
    TString nameD = "areAllGenDMatched";
    arrTTPreselectedMatchingJet[i]
        ->Branch(nameD, &areAllGenDMatched, nameD + "/O", __SIZEOF_INT__)
        ->SetTitle("Whether all the d-quarks matches a jet");
    TString nameDPair = "areAllGenDPairMatched";
    arrTTPreselectedMatchingFatJet[i]
        ->Branch(nameDPair, &areAllGenDPairMatched, nameDPair + "/O",
                 __SIZEOF_INT__)
        ->SetTitle("Whether all the d-pairs matches a fat jet");
  }
  UInt_t nJetMatched;
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "nJetMatched";
    TString title = "Number of unique jets matching a d-quark";
    arrTTPreselectedMatchingJet[i]
        ->Branch(name, &nJetMatched, name + "/I", __SIZEOF_INT__)
        ->SetTitle(title);
    arrTTAllMatchedJet[i]
        ->Branch(name, &nJetMatched, name + "/I", __SIZEOF_INT__)
        ->SetTitle(title);
  }
  UInt_t nFatJetMatched;
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "nFatJetMatched";
    TString title = "Number of unique fat jets matching a d-pair";
    arrTTPreselectedMatchingFatJet[i]
        ->Branch(name, &nFatJetMatched, name + "/I", __SIZEOF_INT__)
        ->SetTitle(title);
    arrTTAllMatchedFatJet[i]
        ->Branch(name, &nFatJetMatched, name + "/I", __SIZEOF_INT__)
        ->SetTitle(title);
  }
  auto lambdaFillTheTrees = [&]() -> void {
    const Bool_t toCollectGen = isSignal;
    const Bool_t toCollectNumCorrect = true;
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
      if (toCollectJetMatching && ((!isSignal) || areAllGenDMatchingPassed) &&
          arrIsPassingZMassCut[i]) {
        arrTTPreselectedMatchingJet[i]->Fill();
      }
      if (toCollectAllMatched && ((!isSignal) || areAllGenDMatchingPassed) &&
          areAllGenDMatched && arrIsPassingZMassCut[i]) {
        arrTTAllMatchedJet[i]->Fill();
      }
      if (toCollectJetMatching && ((!isSignal) || areAllGenDPairMatchingPassed) &&
          arrIsPassingZMassCut[i]) {
        arrTTPreselectedMatchingFatJet[i]->Fill();
      }
      if (toCollectAllMatched && ((!isSignal) || areAllGenDPairMatchingPassed) &&
          areAllGenDPairMatched &&
          arrIsPassingZMassCut[i]) {
        arrTTAllMatchedFatJet[i]->Fill();
      }
    }
  };

  UInt_t nJet, nJetPassed;
  UInt_t nFatJet, nFatJetPassed;
  for (Byte_t i = 0; i < 2; i++) {  // Originally [0, 3)
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
        ->Branch(nameJetPrefix + namePassSuffix, &nJetPassed,
                 nameJetPrefix + namePassSuffix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix + titlePassSuffix);
    arrTTNumCorrect[i]
        ->Branch(nameJetPrefix, &nJet, nameJetPrefix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix);
    arrTTNumCorrect[i]
        ->Branch(nameJetPrefix + namePassSuffix, &nJetPassed,
                 nameJetPrefix + namePassSuffix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix + titlePassSuffix);
    arrTTPreselectedMatchingJet[i]
        ->Branch(nameJetPrefix, &nJet, nameJetPrefix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix);
    arrTTPreselectedMatchingJet[i]
        ->Branch(nameJetPrefix + namePassSuffix, &nJetPassed,
                 nameJetPrefix + namePassSuffix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix + titlePassSuffix);
    arrTTAllMatchedJet[i]
        ->Branch(nameJetPrefix, &nJet, nameJetPrefix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix);
    arrTTAllMatchedJet[i]
        ->Branch(nameJetPrefix + namePassSuffix, &nJetPassed,
                 nameJetPrefix + namePassSuffix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleJetPrefix + titlePassSuffix);
    arrTTGen[i]
        ->Branch(nameFatJetPrefix, &nFatJet, nameFatJetPrefix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix);
    arrTTGen[i]
        ->Branch(nameFatJetPrefix + namePassSuffix, &nFatJetPassed,
                 nameFatJetPrefix + namePassSuffix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix + titlePassSuffix);
    arrTTNumCorrect[i]
        ->Branch(nameFatJetPrefix, &nFatJet, nameFatJetPrefix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix);
    arrTTNumCorrect[i]
        ->Branch(nameFatJetPrefix + namePassSuffix, &nFatJetPassed,
                 nameFatJetPrefix + namePassSuffix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix + titlePassSuffix);
    arrTTPreselectedMatchingFatJet[i]
        ->Branch(nameFatJetPrefix, &nFatJet, nameFatJetPrefix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix);
    arrTTPreselectedMatchingFatJet[i]
        ->Branch(nameFatJetPrefix + namePassSuffix, &nFatJetPassed,
                 nameFatJetPrefix + namePassSuffix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix + titlePassSuffix);
    arrTTAllMatchedFatJet[i]
        ->Branch(nameFatJetPrefix, &nFatJet, nameFatJetPrefix + "/I",
                 __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix);
    arrTTAllMatchedFatJet[i]
        ->Branch(nameFatJetPrefix + namePassSuffix, &nFatJetPassed,
                 nameFatJetPrefix + namePassSuffix + "/I", __SIZEOF_INT__)
        ->SetTitle(titleFatJetPrefix + titlePassSuffix);
  }

  std::vector<UInt_t> vIdxJetPassed;
  std::vector<UInt_t> vIdxFatJetPassed;
  // std::vector<UInt_t> vRankJetPassedPt;
  {
    TString nameIdx = "Jet_jetIdxPassed";
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
  {
    TString nameIdx = "FatJet_fatJetIdxPassed";
    TString titleIdx = "Index of fat jets passing preselections";
    TString nameRank = "FatJet_rankJetPassedPt";
    TString titleRank = "Rank (0-indexed) of pt among passed fat jets";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTGen[i]->Branch(nameIdx, &vIdxFatJetPassed)->SetTitle(titleIdx);
      // arrTTGen[i]->Branch(nameRank, &vRankJetPassedPt)->SetTitle(titleRank);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTNumCorrect[i]
          ->Branch(nameIdx, &vIdxFatJetPassed)
          ->SetTitle(titleIdx);
      // arrTTNumCorrect[i]->Branch(nameRank,
      // &vRankJetPassedPt)->SetTitle(titleRank);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTZMassCutted[i]
          ->Branch(nameIdx, &vIdxFatJetPassed)
          ->SetTitle(titleIdx);
      // arrTTZMassCutted[i]->Branch(nameRank,
      // &vRankJetPassedPt)->SetTitle(titleRank);
    }
  }
  std::vector<Int_t> vDGenIdxJetClosest(4, -1);
  std::vector<Int_t> vDRankJetPassedPt(4, -10);
  std::vector<UInt_t> vJetMatchedRankJetPassedPt;
  {
    TString nameIdx = "GenDMatching_jetIdx";
    TString titleIdx = "Index of closest jets of each D quarks";
    TString nameRank = "GenDMatching_rankJetPassedPt";
    TString titleRank =
        "Rank (0-indexed) of pt of the closest jet of each D quark among "
        "passed jets";
    TString nameRankUnique = "JetMatched_rankJetPassedPt";
    TString titleRankUnique =
        "Rank (0-indexed) of pt of each matched jet among passed jets";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatchingJet[i]
          ->Branch(nameIdx, &vDGenIdxJetClosest)
          ->SetTitle(titleIdx);
      arrTTPreselectedMatchingJet[i]
          ->Branch(nameRank, &vDRankJetPassedPt)
          ->SetTitle(titleRank);
      arrTTPreselectedMatchingJet[i]
          ->Branch(nameRankUnique, &vJetMatchedRankJetPassedPt)
          ->SetTitle(titleRankUnique);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTAllMatchedJet[i]
          ->Branch(nameIdx, &vDGenIdxJetClosest)
          ->SetTitle(titleIdx);
      arrTTAllMatchedJet[i]
          ->Branch(nameRank, &vDRankJetPassedPt)
          ->SetTitle(titleRank);
      arrTTAllMatchedJet[i]
          ->Branch(nameRankUnique, &vJetMatchedRankJetPassedPt)
          ->SetTitle(titleRankUnique);
    }
  }
  std::vector<Int_t> vGenDPairIdxFatJetClosest(2, -1);
  std::vector<Int_t> vDPairRankFatJetPassedPt(2, -10);
  std::vector<UInt_t> vFatJetMatchedRankFatJetPassedPt;
  {
    TString nameIdx = "GenDPairMatching_fatJetIdx";
    TString titleIdx = "Index of closest fat jets of each d-pairs";
    TString nameRank = "GenDPairMatching_rankFatJetPassedPt";
    TString titleRank =
        "Rank (0-indexed) of pt of the closest jet of each d-pair among "
        "passed fat jets";
    TString nameRankUnique = "FatJetMatched_rankFatJetPassedPt";
    TString titleRankUnique =
        "Rank (0-indexed) of pt of each matched fat jet among passed fat jets";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatchingFatJet[i]
          ->Branch(nameIdx, &vGenDPairIdxFatJetClosest)
          ->SetTitle(titleIdx);
      arrTTPreselectedMatchingFatJet[i]
          ->Branch(nameRank, &vDPairRankFatJetPassedPt)
          ->SetTitle(titleRank);
      arrTTPreselectedMatchingFatJet[i]
          ->Branch(nameRankUnique, &vFatJetMatchedRankFatJetPassedPt)
          ->SetTitle(titleRankUnique);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTAllMatchedFatJet[i]
          ->Branch(nameIdx, &vGenDPairIdxFatJetClosest)
          ->SetTitle(titleIdx);
      arrTTAllMatchedFatJet[i]
          ->Branch(nameRank, &vDPairRankFatJetPassedPt)
          ->SetTitle(titleRank);
      arrTTAllMatchedFatJet[i]
          ->Branch(nameRankUnique, &vFatJetMatchedRankFatJetPassedPt)
          ->SetTitle(titleRankUnique);
    }
  }
  Bool_t areAllJetsMatchedLeading;
  {
    TString name = "areAllJetsMatchedLeading";
    TString title =
        "Whether all the jets matched are leading jets ( (rankMax == "
        "nJetsMatched-1)";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatchingJet[i]
          ->Branch(name, &areAllJetsMatchedLeading, name + "/O", 1)
          ->SetTitle(title);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTAllMatchedJet[i]
          ->Branch(name, &areAllJetsMatchedLeading, name + "/O", 1)
          ->SetTitle(title);
    }
  }
  Bool_t areAllFatJetsMatchedLeading;
  {
    TString name = "areAllFatJetsMatchedLeading";
    TString title =
        "Whether all the fat jets matched are leading jets (rankMax == "
        "nJetsMatched-1)";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatchingFatJet[i]
          ->Branch(name, &areAllFatJetsMatchedLeading, name + "/O", 1)
          ->SetTitle(title);
    }
    for (Byte_t i = 0; i < 2; i++) {
      arrTTAllMatchedJet[i]
          ->Branch(name, &areAllFatJetsMatchedLeading, name + "/O", 1)
          ->SetTitle(title);
    }
  }
  std::vector<Float_t> vDeltaRDJet(4, -1);
  std::vector<Float_t> vDeltaRJetPair(2, -1);
  std::vector<Float_t> vDeltaRBetweenJetPairs(1, -1);
  {
    TString name = "GenDMatching_deltaRJet";
    TString title = "DeltaR between each d quark and its closest jet";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatchingJet[i]
          ->Branch(name, &vDeltaRDJet)
          ->SetTitle(title);
    }
  }
  {
    TString name = "GenDMatchingPair_deltaRJet";
    TString title = "DeltaR of each jet pair most likely to match the d quarks";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatchingJet[i]
          ->Branch(name, &vDeltaRJetPair)
          ->SetTitle(title);
    }
  }
  {
    TString name = "GenDMatchingPairPair_deltaRJet";
    TString title =
        "DeltaR between jet pairs most likely to match the d quarks";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatchingJet[i]
          ->Branch(name, &vDeltaRBetweenJetPairs)
          ->SetTitle(title);
    }
  }
  std::vector<Float_t> vDeltaRDPairFatJet(2, -1);
  {
    const TString name = "GenDPairMatching_deltaRFatJet";
    const TString title = "DeltaR between each d-pair and its closest fat jet";
    for (Byte_t i = 0; i < 2; i++) {
      arrTTPreselectedMatchingFatJet[i]
          ->Branch(name, &vDeltaRDPairFatJet)
          ->SetTitle(title);
    }
  }
  DescriptionCollectorForUntuplizer collectorGenDMatching;
  DescriptionCollectorForUntuplizer collectorHT;
  DescriptionCollectorForUntuplizer collectorMET;
  BranchDescription descriptionNMETJet;
  DescriptionCollectorForUntuplizer collectorMETJet;
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
      if (descriptionCurrent.name.Contains("MET_")) {
        collectorMET.BranchFor(descriptionCurrent);
      }
      if (descriptionCurrent.name.BeginsWith("CorrT1METJet_")) {
        collectorMETJet.BranchFor(descriptionCurrent);
      }
      if (descriptionCurrent.name.EqualTo("nCorrT1METJet")) {
        descriptionNMETJet = descriptionCurrent;
      }
      if (descriptionCurrent.name.BeginsWith("GenPart_") &&
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
  collectorMET.Prepare();
  collectorMETJet.Prepare();
  collectorJet.Prepare();
  collectorFatJet.Prepare();
  BranchMounterIterableForUntuplizer mounterGenDMatching(
      collectorGenDMatching, data, [](BranchDescription descriptionCurrent) {
        Ssiz_t lenName = descriptionCurrent.name.Length();
        const Ssiz_t lenPrefOriginal = 8;
        return BranchDescription(
            "GenDMatching_" +
                descriptionCurrent.name(lenPrefOriginal, lenName - 1),
            descriptionCurrent.title, descriptionCurrent.typeName);
      });
  // BranchMounterIterableForUntuplizer
  // mounterGenDMatching(collectorGenDMatching, data);
  BranchMounterScalarForUntuplizer mounterHT(collectorHT, data);
  BranchMounterScalarForUntuplizer mounterMET(collectorMET, data);
  BranchMounterScalarSingle<Int_t> mounterNMETJet(std::vector<BranchDescription>({descriptionNMETJet}), [&data](BranchDescription description){return data.GetInt(description.name);});
  BranchMounterIterableForUntuplizer mounterMETJet(collectorMETJet, data);
  BranchMounterIterableForUntuplizer mounterJet(collectorJet, data);
  BranchMounterIterableForUntuplizer mounterFatJet(collectorFatJet, data);
  {
    auto funDoIntersection = [&mounterHT, &mounterMET, &mounterNMETJet, &mounterMETJet](TTree* tree) {
      mounterHT.BranchOn(tree);
      mounterMET.BranchOn(tree);
      mounterNMETJet.BranchOn(tree);
      mounterMETJet.BranchOn(tree);
    };
    auto funDoJet = [&mounterGenDMatching, &mounterJet](TTree* tree) {
      if (!TString(tree->GetName()).BeginsWith("NumCorrect")) {
        mounterGenDMatching.BranchOn(tree);
      }
      mounterJet.BranchOn(tree);
    };
    auto funDoFatJet = [&mounterFatJet](TTree* tree) {
      mounterFatJet.BranchOn(tree);
    };
    for (TTree* const tree : arrTTGen) {
      funDoIntersection(tree);
      funDoJet(tree);
      funDoFatJet(tree);
    }
    for (TTree* const tree : arrTTPreselectedMatchingJet) {
      funDoIntersection(tree);
      funDoJet(tree);
    }
    for (TTree* const tree : arrTTAllMatchedJet) {
      funDoIntersection(tree);
      funDoJet(tree);
    }
    for (TTree* const tree : arrTTPreselectedMatchingFatJet) {
      funDoIntersection(tree);
      funDoFatJet(tree);
    }
    for (TTree* const tree : arrTTAllMatchedFatJet) {
      funDoIntersection(tree);
      funDoFatJet(tree);
    }
  }

  for (jEntry = 0; jEntry < nEntry; jEntry++) {
    data.GetEntry(jEntry);

    nJet = data.GetInt("nJet");
    nFatJet = data.GetInt("nFatJet");
    // if (debug) std::cout << "nJet: " << nJet << std::endl;

    areAllGenDMatchingPassed = false;
    areAllGenDPairMatchingPassed = false;
    areAllGenDPairMatchingPassed = false;
    areAllGenDMatched = false;
    areAllGenDPairMatched = false;
    areAllJetsMatchedLeading = false;
    areAllFatJetsMatchedLeading = false;
    nGenDMatchingPassed = 0;
    nGenDPairMatchingPassed = 0;
    nGenDMatched = 0;
    nGenDPairMatched = 0;
    nJetPassed = 0;
    nFatJetPassed = 0;
    nJetMatched = 0;
    nFatJetMatched = 0;
    std::fill(std::begin(arrRecoIsLeptonNumCorrect), std::end(arrRecoIsLeptonNumCorrect), false);
    std::fill(std::begin(arrLeptonPair_pt), std::end(arrLeptonPair_pt), -1);
    std::fill(std::begin(arrLeptonPair_eta), std::end(arrLeptonPair_eta), -1000);
    std::fill(std::begin(arrLeptonPair_phi), std::end(arrLeptonPair_phi), -10);
    std::fill(std::begin(arrLeptonPair_mass), std::end(arrLeptonPair_mass), -1);
    // Allocated after declearation (2 or 4 elements)
    for (auto& vLeptonIdxNumCorrect: arrVLeptonIdxNumCorrect) {
      std::fill(vLeptonIdxNumCorrect.begin(), vLeptonIdxNumCorrect.end(), -1);
    }
    std::fill(vGenDMatchingPt.begin(), vGenDMatchingPt.end(), -1);
    std::fill(vGenDMatchingEta.begin(), vGenDMatchingEta.end(), -1000);
    std::fill(vGenDMatchingPhi.begin(), vGenDMatchingPhi.end(), -10);
    std::fill(vGenDPairMatchingPt.begin(), vGenDPairMatchingPt.end(), -1);
    std::fill(vGenDPairMatchingEta.begin(), vGenDPairMatchingEta.end(), -1000);
    std::fill(vGenDPairMatchingPhi.begin(), vGenDPairMatchingPhi.end(), -10);
    std::fill(vGenDPairMatchingMass.begin(), vGenDPairMatchingMass.end(), -1);
    std::fill(vDeltaRDJet.begin(), vDeltaRDJet.end(), -1);
    std::fill(vDeltaRJetPair.begin(), vDeltaRJetPair.end(), -1);
    std::fill(vDeltaRBetweenJetPairs.begin(), vDeltaRBetweenJetPairs.end(), -1);
    std::fill(vDeltaRDPairFatJet.begin(), vDeltaRDPairFatJet.end(), -1);
    std::fill(vDGenIdxJetClosest.begin(), vDGenIdxJetClosest.end(), -1);
    std::fill(vGenDPairIdxFatJetClosest.begin(),
              vGenDPairIdxFatJetClosest.end(), -1);
    std::fill(vGenDMatchedId.begin(), vGenDMatchedId.end(), false);
    std::fill(vGenDPairMatchedId.begin(), vGenDPairMatchedId.end(), false);
    vIdxJetPassed.clear();
    vIdxFatJetPassed.clear();
    vJetMatchedRankJetPassedPt.clear();
    vFatJetMatchedRankFatJetPassedPt.clear();

    areAllGenDMatchingPassed = false;
    areAllGenDPairMatchingPassed = false;

    Bool_t hasGenLepton = false;
    UInt_t nGenPart = data.GetInt("nGenPart");
    Int_t* ptrGenPart_pdgId = data.GetPtrInt("GenPart_pdgId");
    Int_t* ptrGenPart_genPartIdxMother =
        data.GetPtrInt("GenPart_genPartIdxMother");
    for (Byte_t i = 0; i < 3; i++) arrHasGenLepton[i] = false;

    haveAllGenDMatching = true;
    Bool_t arrHasGenD[4];
    std::fill(std::begin(arrHasGenD), std::end(arrHasGenD), false);
    std::fill(std::begin(vGenDMatchingIdx), std::end(vGenDMatchingIdx), -1);
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
      // if (debug) std::cout << "jEntry: " << jEntry << "\n";
      if (TMath::Abs(ptrGenPart_pdgId[ptrGenPart_genPartIdxMother[ig]]) ==
              pdgX2 &&
          TMath::Abs(ptrGenPart_pdgId[ig]) == pdgDown) {
        Byte_t iD = (static_cast<Byte_t>(signbit(
                         ptrGenPart_pdgId[ptrGenPart_genPartIdxMother[ig]]))
                     << 1) +
                    (signbit(ptrGenPart_pdgId[ig]));
        if (!arrHasGenD[iD]) {
          // if (debug)
          //   std::cout << "found d-quark no. " << static_cast<Int_t>(iD) << "
          //   at "
          //             << ig << "\n";
          arrHasGenD[iD] = true;
          vGenDMatchingIdx[iD] = ig;
        }
      }
    }
    for (Byte_t iD = 0; iD < 4; iD++) {
      if (!arrHasGenD[iD]) {
        haveAllGenDMatching = false;
        if (debug)
          std::cout << "d quark " << static_cast<Int_t>(iD) << "isn't found."
                    << "\n";
        break;
      }
    }
    // if (debug) std::cout << std::flush;
    Bool_t isGenLongLivedEvent = hasGenLepton && haveAllGenDMatching;
    {
      // if (BranchMounterHelper::RunForEachNIdxSorted<
      //         BranchMounterIterableForUntuplizer,
      //         std::vector<Int_t>::iterator, Int_t>(mounterGenDMatching,
      //         vGenDMatchingIdx.begin(), 4) < 4 &&
      //     debug) {
      //   Error("BranchMounterHelper::RunForEachIdxSorted",
      //         "Not all d indices found!");
      // }
      std::vector<Byte_t> idxSortVGenDMatchingIdx(4, 0);
      std::iota(idxSortVGenDMatchingIdx.begin(), idxSortVGenDMatchingIdx.end(),
                0);
      std::sort(idxSortVGenDMatchingIdx.begin(), idxSortVGenDMatchingIdx.end(),
                [&vGenDMatchingIdx](Byte_t i, Byte_t j) -> Bool_t {
                  return vGenDMatchingIdx[i] < vGenDMatchingIdx[j];
                });
      Int_t iii = 0;
      BranchMounterHelper::RunOnFunAt(
          mounterGenDMatching,
          static_cast<std::function<Bool_t(size_t)>>([
            &vGenDMatchingIdx, &idxSortVGenDMatchingIdx, &iii
          ](size_t idx)
              ->Bool_t {
                if (idx ==
                    static_cast<size_t>(
                        vGenDMatchingIdx[idxSortVGenDMatchingIdx[iii]])) {
                  iii++;
                  return true;
                }
                return false;
              }),
          static_cast<size_t>(nGenPart));
      if (debug && outputFileTree.Contains("_ctau-") && iii != 4) {
        Fatal("xAna_monoZ_preselect",
              "mounterGenDMatching failed to collect four values. iii = %d",
              iii);
      }
      for (Byte_t j1 = 0, j2 = 0; j1 < 4; j2 = ++j1) {
        while (j2 != idxSortVGenDMatchingIdx[j2]) {
          mounterGenDMatching.SwapAt(
              static_cast<size_t>(j2),
              static_cast<size_t>(idxSortVGenDMatchingIdx[j2]));
          std::swap(j2, idxSortVGenDMatchingIdx[j2]);
        }
      }
      if (debug) {
        for (Byte_t i = 0; i < 4; i++) {
          if (idxSortVGenDMatchingIdx[i] != i) {
            Fatal("xAna_monoZ_preselect",
                  "permutation appling algorithom has a bug.");
          }
        }
      }
    }


    // Float_t* ptrElectron_pt = data.GetPtrFloat("Electron_pt");
    // Float_t* ptrElectron_phi = data.GetPtrFloat("Electron_phi");
    // Float_t* ptrElectron_eta = data.GetPtrFloat("Electron_eta");
    // Float_t* ptrElectron_eCorr = data.GetPtrFloat("Electron_eCorr");

    // mounterElectronDynamics.PrepareForPushing(nElectron);
    // for (Byte_t i = 0; i < nElectron; i++) mounterElectronDynamics.PushOrSkip(true);

    Int_t nElectron = data.GetInt("nElectron");
    Int_t nMuon = data.GetInt("nMuon");
    for (Byte_t i=0; i<2; i++) {
      Int_t nLepton = i ? nMuon : nElectron;
      vMounterLeptonChained[i].PrepareForPushing(nLepton);
      if (debug) std::cout << (i ? "nMuon: " : "nElectron: ") << (i ? nMuon : nElectron) << std::endl;
      for (Byte_t j=0; j<nLepton; j++) {
        // Pt > 20 and |Eta| < 4.5
        vMounterLeptonChained[i].PushOrSkip(
            vMounterLeptonDynamics[i].Peek(0) > 20 &&
            TMath::Abs(vMounterLeptonDynamics[i].Peek(1)) < 4.5);
      }
    }
    Int_t nElectronPassedPtEta = vMounterLeptonDynamics[0].vvE.back().size();
    Int_t nMuonPassedPtEta = vMounterLeptonDynamics[1].vvE.back().size();
    if (debug) std::cout << "nElectronPassedPtEta: " << nElectronPassedPtEta << "\n"
        "nMuonPassedPtEta: " << nMuonPassedPtEta << std::endl;

    std::vector<Int_t> arrVLeptonIdxPassedPtEtaNumCorrect[2] = {{-1, -1}, {-1, -1}};

    Bool_t* ptrElectron_mvaFall17V2Iso_WPL =
        data.GetPtrBool("Electron_mvaFall17V2Iso_WPL");
    Bool_t* ptrElectron_mvaFall17V2Iso_WP90 =
        data.GetPtrBool("Electron_mvaFall17V2Iso_WP90");
    Bool_t* ptrElectron_mvaFall17V2Iso_WP80 =
        data.GetPtrBool("Electron_mvaFall17V2Iso_WP80");
    Bool_t* ptrElectron_isSoft = ptrElectron_mvaFall17V2Iso_WP80;
    Bool_t* ptrElectron_isTight = ptrElectron_mvaFall17V2Iso_WP90;
    Bool_t* ptrMuon_softId = data.GetPtrBool("Muon_softId");
    Bool_t* ptrMuon_tightId = data.GetPtrBool("Muon_tightId");
    Bool_t* ptrMuon_isSoft = ptrMuon_softId;
    Bool_t* ptrMuon_isTight = ptrMuon_tightId;

    Bool_t isNumElectronCorrect = true;
    for (Int_t iElectronPassed = 0, iiElectronNumCorrect = 0; iElectronPassed < nElectronPassedPtEta; iElectronPassed++) {
      const Int_t iElectron = vMounterLeptonIdxPassedPtEta[0].vvE.back()[iElectronPassed];
      if (ptrElectron_isSoft[iElectron]) {
        if (iiElectronNumCorrect == 2) {
          isNumElectronCorrect = false;
          std::fill(arrVLeptonIdxNumCorrect[0].begin(), arrVLeptonIdxNumCorrect[0].end(), -1);
          break;
        }
        arrVLeptonIdxNumCorrect[0][iiElectronNumCorrect] = iElectron;
        arrVLeptonIdxPassedPtEtaNumCorrect[0][iiElectronNumCorrect] = iElectronPassed;
        iiElectronNumCorrect++;
      }
    }
    if (arrVLeptonIdxNumCorrect[0].back() < 0) {
      isNumElectronCorrect = false;
      std::fill(arrVLeptonIdxNumCorrect[0].begin(), arrVLeptonIdxNumCorrect[0].end(), -1);
    }
    // nElectronPair = isNumElectronCorrect;

    Bool_t isNumMuonCorrect = true;
    for (Int_t iMuonPassed = 0, iiMuonNumCorrect = 0; iMuonPassed < nMuonPassedPtEta; iMuonPassed++) {
      const Int_t iMuon = vMounterLeptonIdxPassedPtEta[1].vvE.back()[iMuonPassed];
      if (ptrMuon_isSoft[iMuon]) {
        if (iiMuonNumCorrect == 2) {
          isNumMuonCorrect = false;
          std::fill(arrVLeptonIdxNumCorrect[1].begin(), arrVLeptonIdxNumCorrect[1].end(), -1);
          break;
        }
        arrVLeptonIdxNumCorrect[1][iiMuonNumCorrect] = iMuon;
        arrVLeptonIdxPassedPtEtaNumCorrect[1][iiMuonNumCorrect] = iMuonPassed;
        iiMuonNumCorrect++;
      }
    }
    if (arrVLeptonIdxNumCorrect[1].back() < 0) {
      isNumMuonCorrect = false;
      std::fill(arrVLeptonIdxNumCorrect[1].begin(), arrVLeptonIdxNumCorrect[1].end(), -1);
    }

    Bool_t isNumCorrect = (isNumElectronCorrect != isNumMuonCorrect);
    arrRecoIsLeptonNumCorrect[0] = isNumCorrect && isNumElectronCorrect;
    arrRecoIsLeptonNumCorrect[1] = isNumCorrect && isNumMuonCorrect;
    if (debug) std::cout << "arrRecoIsLeptonNumCorrect = {" << arrRecoIsLeptonNumCorrect[0] << ", " << arrRecoIsLeptonNumCorrect[1] << "}" << std::endl;

    // if (!isNumCorrect) {
    //   lambdaFillTheTrees();
    //   continue;
    // }
    if (debug)
      std::cout << "elepair, mupair: " << (Int_t)isNumElectronCorrect
                << (Int_t)isNumMuonCorrect << std::endl;

    for (Byte_t i = 0; i < 2; i++) arrIsPassingZMassCut[i] = false;


    TLorentzVector* ppElectronP4NumberCorrect[2];
    TLorentzVector* ptrElectronP4NumberCorrectSum;

    if (arrRecoIsLeptonNumCorrect[0]) {
      ptrElectronP4NumberCorrectSum = new TLorentzVector();
      for (Int_t iiElectronNumCorrect = 0; iiElectronNumCorrect < 2; iiElectronNumCorrect++) {
        const auto& indexElectronPassedPtEtaNumCorrect = arrVLeptonIdxPassedPtEtaNumCorrect[0][iiElectronNumCorrect];
        ppElectronP4NumberCorrect[iiElectronNumCorrect] = new TLorentzVector();
        // ppElectronP4NumberCorrect[i]->SetPtEtaPhiE(
        //     ptrElectron_pt[i], ptrElectron_eta[i], ptrElectron_phi[i],
        //     ptrElectron_eCorr[i]);
        ppElectronP4NumberCorrect[iiElectronNumCorrect]->SetPtEtaPhiM(
            vMounterLeptonDynamics[0].vvE[0][indexElectronPassedPtEtaNumCorrect],
            vMounterLeptonDynamics[0].vvE[1][indexElectronPassedPtEtaNumCorrect],
            vMounterLeptonDynamics[0].vvE[2][indexElectronPassedPtEtaNumCorrect],
            massElectron);
        *ptrElectronP4NumberCorrectSum += *ppElectronP4NumberCorrect[iiElectronNumCorrect];
      }
      // Double_t electronPair_mass = ptrElectronP4NumberCorrectSum->M();
      // vElectronPair_massNumCorrect.push_back(electronPair_mass);
      if (arrRecoIsLeptonNumCorrect[0]) {
        arrLeptonPair_pt[0] = ptrElectronP4NumberCorrectSum->Pt();
        arrLeptonPair_eta[0] = ptrElectronP4NumberCorrectSum->Eta();
        arrLeptonPair_phi[0] = ptrElectronP4NumberCorrectSum->Phi();
        arrLeptonPair_mass[0] = ptrElectronP4NumberCorrectSum->M();
      }

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

    // Float_t* ptrMuon_pt = data.GetPtrFloat("Muon_pt");
    // Float_t* ptrMuon_phi = data.GetPtrFloat("Muon_phi");
    // Float_t* ptrMuon_eta = data.GetPtrFloat("Muon_eta");
    // Float_t* ptrMuon_mass = data.GetPtrFloat("Muon_mass");

    TLorentzVector* ppMuonP4NumberCorrect[2];
    TLorentzVector* ptrMuonP4NumberCorrectSum;
    if (arrRecoIsLeptonNumCorrect[1]) {
      ptrMuonP4NumberCorrectSum = new TLorentzVector();
      for (Int_t iiMuonNumCorrect = 0; iiMuonNumCorrect < 2; iiMuonNumCorrect++) {
        const auto& indexMuonPassedPtEtaNumCorrect = arrVLeptonIdxPassedPtEtaNumCorrect[1][iiMuonNumCorrect];
        ppMuonP4NumberCorrect[iiMuonNumCorrect] = new TLorentzVector();
        ppMuonP4NumberCorrect[iiMuonNumCorrect]->SetPtEtaPhiM(
            vMounterLeptonDynamics[1].vvE[0][indexMuonPassedPtEtaNumCorrect],
            vMounterLeptonDynamics[1].vvE[1][indexMuonPassedPtEtaNumCorrect],
            vMounterLeptonDynamics[1].vvE[2][indexMuonPassedPtEtaNumCorrect],
            vMounterLeptonDynamics[1].vvE[3][indexMuonPassedPtEtaNumCorrect]);
        *ptrMuonP4NumberCorrectSum += *ppMuonP4NumberCorrect[iiMuonNumCorrect];
      }
      // Double_t muonPair_mass = ptrMuonP4NumberCorrectSum->M();
      // vMuonPair_massNumCorrect.push_back(muonPair_mass);
      if (arrRecoIsLeptonNumCorrect[1]) {
        arrLeptonPair_pt[1] = ptrMuonP4NumberCorrectSum->Pt();
        arrLeptonPair_eta[1] = ptrMuonP4NumberCorrectSum->Eta();
        arrLeptonPair_phi[1] = ptrMuonP4NumberCorrectSum->Phi();
        arrLeptonPair_mass[1] = ptrMuonP4NumberCorrectSum->M();
      }

      // leptonPair_massNumCorrect = muonPair_mass;
      Double_t massZCutUpper = massZ + 20;
      Double_t massZCutLower = massZ - 20;
      // Z mass cut for muon pairs
      arrIsPassingZMassCut[1] = (massZCutLower < arrLeptonPair_mass[1] &&
                                 arrLeptonPair_mass[1] < massZCutUpper);
      if (debug)
        std::cout << "massLPair: " << arrLeptonPair_mass[1] << std::endl;
    }

    // if (!arrIsPassingZMassCut[0] && !arrIsPassingZMassCut[1] &&
    //     !isGenLongLivedEvent) {
    //   lambdaFillTheTrees();
    //   continue;  // TODO
    // }

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

    if (debug) std::printf("Matching jets");
    // Jet_*
    vIdxJetPassed.clear();
    vIdxJetPassed.reserve(nJet);
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
    vIdxJetPassed.shrink_to_fit();
    vIdxFatJetPassed.clear();
    vIdxFatJetPassed.reserve(nFatJet);
    Float_t* ptrFatJet_pt_original = data.GetPtrFloat("FatJet_pt");
    Float_t* ptrFatJet_eta_original = data.GetPtrFloat("FatJet_eta");
    Float_t* ptrFatJet_phi_original = data.GetPtrFloat("FatJet_phi");
    Float_t* ptrFatJet_mass_original = data.GetPtrFloat("FatJet_mass");
    mounterFatJet.PrepareForPushing(nFatJet);
    {
      auto ptrFatJet_pt_original_cloned = ptrFatJet_pt_original;
      auto ptrFatJet_eta_original_cloned = ptrFatJet_eta_original;
      // auto ptrFatJet_phi_original_cloned = ptrFatJet_phi_original;
      // auto ptrFatJet_mass_original_cloned = ptrFatJet_mass_original;
      for (UInt_t idxFatJet = 0; idxFatJet < nFatJet; idxFatJet++) {
        Bool_t result = *ptrFatJet_pt_original_cloned > 200 &&
                        TMath::Abs(*ptrFatJet_eta_original_cloned) < 4.5;
        if (result) {
          vIdxFatJetPassed.push_back(idxFatJet);
        }
        mounterFatJet.PushOrSkip(result);
        ptrFatJet_pt_original_cloned++;
        ptrFatJet_eta_original_cloned++;
      }
    }
    nFatJetPassed = vIdxFatJetPassed.size();
    vIdxFatJetPassed.shrink_to_fit();
    mounterHT.Push();
    mounterMET.Push();
    mounterNMETJet.Push();
    {
      const Int_t &nMETJet = mounterNMETJet.vE.back();
      mounterMETJet.PrepareForPushing(nMETJet);
      for (Int_t i=0; i<nMETJet; i++) {
        mounterMETJet.PushOrSkip(true);
      }
    }

    if (!isGenLongLivedEvent) {
      lambdaFillTheTrees();
      continue;
    }

    Float_t* ptrGenPart_pt = data.GetPtrFloat("GenPart_pt");
    Float_t* ptrGenPart_eta = data.GetPtrFloat("GenPart_eta");
    Float_t* ptrGenPart_phi = data.GetPtrFloat("GenPart_phi");
    TLorentzVector* p4DMatching[4];
    TLorentzVector* p4DJetClosest[4];
    for (Byte_t iDMatching = 0; iDMatching < 4; iDMatching++) {
      if (debug) std::printf("iDMatching: %d\n", iDMatching);
      p4DJetClosest[iDMatching] = nullptr;
      vDeltaRDJet[iDMatching] = __FLT_MAX__;
      p4DMatching[iDMatching] = new TLorentzVector;
      vGenDMatchingPt[iDMatching] = ptrGenPart_pt[vGenDMatchingIdx[iDMatching]];
      vGenDMatchingEta[iDMatching] =
          ptrGenPart_eta[vGenDMatchingIdx[iDMatching]];
      vGenDMatchingPhi[iDMatching] =
          ptrGenPart_phi[vGenDMatchingIdx[iDMatching]];
      vGenDMatchingIdPassed[iDMatching] =
          vGenDMatchingPt[iDMatching] > 30 &&
          TMath::Abs(vGenDMatchingEta[iDMatching]) < 3;
      if (vGenDMatchingIdPassed[iDMatching]) {
        nGenDMatchingPassed++;
      }
      p4DMatching[iDMatching]->SetPtEtaPhiM(
          vGenDMatchingPt[iDMatching], vGenDMatchingEta[iDMatching],
          vGenDMatchingPhi[iDMatching], massDown);
      if (nJetPassed < 2) {
        continue;
      }
      for (Byte_t rankJetPassed = 0; rankJetPassed < nJetPassed;
           rankJetPassed++) {
        if (debug) std::printf("rankJetPassed: %d\n", rankJetPassed);
        TLorentzVector* p4Jet = new TLorentzVector;
        p4Jet->SetPtEtaPhiM(ptrJet_pt_original[vIdxJetPassed[rankJetPassed]],
                            ptrJet_eta_original[vIdxJetPassed[rankJetPassed]],
                            ptrJet_phi_original[vIdxJetPassed[rankJetPassed]],
                            ptrJet_mass_original[vIdxJetPassed[rankJetPassed]]);
        Double_t deltaRDJet = p4DMatching[iDMatching]->DeltaR(*p4Jet);
        if (!iDMatching || deltaRDJet < vDeltaRDJet[iDMatching]) {
          vDeltaRDJet[iDMatching] = deltaRDJet;
          std::swap(p4DJetClosest[iDMatching], p4Jet);
          vDRankJetPassedPt[iDMatching] = rankJetPassed;
          vDGenIdxJetClosest[iDMatching] = vIdxJetPassed[rankJetPassed];
        } else {
        }
        if (iDMatching) {
          delete p4Jet;
        }
      }
      if (debug) std::printf("Calculating jet deltaR\n");
      if (vDeltaRDJet[iDMatching] < 0.4) {
        vGenDMatchedId[iDMatching] = true;
        vJetMatchedRankJetPassedPt.push_back(vDRankJetPassedPt[iDMatching]);
        nGenDMatched++;
      } else {
        vGenDMatchedId[iDMatching] = false;
      }
    }
    if (nJetPassed >= 2) {
      if (debug && vJetMatchedRankJetPassedPt.size() > nDQuarksExpected) {
        Fatal(Form("xAna_monoZ_preselected, line %d", __LINE__), "Size of vJetMatchedRankJetPassedPt (%zu) is greater than %d", vJetMatchedRankJetPassedPt.size(), nDQuarksExpected);
      }
      std::sort(vJetMatchedRankJetPassedPt.begin(),
                vJetMatchedRankJetPassedPt.end());
      nJetMatched =
          std::distance(vJetMatchedRankJetPassedPt.begin(),
                        std::unique(vJetMatchedRankJetPassedPt.begin(),
                                    vJetMatchedRankJetPassedPt.end()));
      vJetMatchedRankJetPassedPt.resize(nJetMatched);
      if (nJetMatched) {
        areAllJetsMatchedLeading =
            vJetMatchedRankJetPassedPt.back() == nJetMatched - 1;
      } else {
        nJetMatched = 0;
        areAllJetsMatchedLeading = false;
      }

      areAllGenDMatchingPassed = nGenDMatchingPassed == 4;
      areAllGenDMatched = nGenDMatched == 4;
      if (debug)
        std::printf("areAllGenDMatchingPassed: %d\n", areAllGenDMatchingPassed);
      if (debug) std::printf("areAllGenDMatched: %d\n", areAllGenDMatched);
      // if (debug) std::cout << "line: " << __LINE__ << " jEntry: " << jEntry << std::endl;
      if (debug) for (Byte_t i=0; i<4; i++) {
        std::cout
        << "p4DJetClosest[" << static_cast<UShort_t>(i) << "] ("<< p4DJetClosest[i] << ") (px, py, pz, e) ("
        << p4DJetClosest[i]->Px() << ", "
        << p4DJetClosest[i]->Py() << ", "
        << p4DJetClosest[i]->Pz() << ", "
        << p4DJetClosest[i]->E()
        << ")" << std::endl;
      }
      for (Byte_t i = 0; i < 2; i++) {
        vDeltaRJetPair[i] =
            p4DJetClosest[i << 1]->DeltaR(*p4DJetClosest[(i << 1) + 1]);
      }
      vDeltaRBetweenJetPairs[0] =
          (*p4DJetClosest[0] + *p4DJetClosest[1])
              .DeltaR(*p4DJetClosest[2] + *p4DJetClosest[3]);
    }
    TLorentzVector* p4DPairMatching[2];
    TLorentzVector* p4DPairFatJetClosest[2];
    for (UInt_t iDPairMatching = 0; iDPairMatching < 2; iDPairMatching++) {
      // std::cout << vGenDPairMatchingMass.size() << '\n';
      // std::cout << (iDPairMatching << 1) << ((iDPairMatching << 1) + 1) <<
      // "\n"; std::cout << p4DMatching[(iDPairMatching << 1)+1]->Pt() <<'\n';
      // if (debug) std::cout << "line: " << __LINE__ << " jEntry: " << jEntry << std::endl;
      p4DPairMatching[iDPairMatching] = new TLorentzVector;
      *(p4DPairMatching[iDPairMatching]) =
          *(p4DMatching[iDPairMatching << 1]) +
          *(p4DMatching[(iDPairMatching << 1) + 1]);
      // if (debug) std::cout << "line: " << __LINE__ << " jEntry: " << jEntry << std::endl;
      // std::cout << "testtag" << std::endl;
      vGenDPairMatchingPt[iDPairMatching] =
          p4DPairMatching[iDPairMatching]->Pt();
      vGenDPairMatchingEta[iDPairMatching] =
          p4DPairMatching[iDPairMatching]->Eta();
      vGenDPairMatchingPhi[iDPairMatching] =
          p4DPairMatching[iDPairMatching]->Phi();
      vGenDPairMatchingMass[iDPairMatching] =
          p4DPairMatching[iDPairMatching]->M();
      vGenDPairMatchingIdPassed[iDPairMatching] =
          vGenDPairMatchingPt[iDPairMatching] > 200 &&
          TMath::Abs(vGenDPairMatchingEta[iDPairMatching]) < 4.5;
      if (vGenDPairMatchingIdPassed[iDPairMatching]) {
        nGenDPairMatchingPassed++;
      }
      if (nFatJetPassed < 2) {
        continue;
      }
      // Fat jet matching
      for (Byte_t rankFatJetPassed = 0; rankFatJetPassed < nFatJetPassed;
           rankFatJetPassed++) {
        if (debug) std::printf("rankFatJetPassed: %d\n", rankFatJetPassed);
        TLorentzVector* p4FatJet = new TLorentzVector;
        p4FatJet->SetPtEtaPhiM(
            ptrFatJet_pt_original[vIdxFatJetPassed[rankFatJetPassed]],
            ptrFatJet_eta_original[vIdxFatJetPassed[rankFatJetPassed]],
            ptrFatJet_phi_original[vIdxFatJetPassed[rankFatJetPassed]],
            ptrFatJet_mass_original[vIdxFatJetPassed[rankFatJetPassed]]);
        Double_t deltaRDPairFatJet =
            p4DPairMatching[iDPairMatching]->DeltaR(*p4FatJet);
        if (!rankFatJetPassed ||
            deltaRDPairFatJet < vDeltaRDPairFatJet[iDPairMatching]) {
          vDeltaRDPairFatJet[iDPairMatching] = deltaRDPairFatJet;
          vDPairRankFatJetPassedPt[iDPairMatching] = rankFatJetPassed;
          std::swap(p4DPairFatJetClosest[iDPairMatching], p4FatJet);
          vGenDPairIdxFatJetClosest[iDPairMatching] =
              vIdxFatJetPassed[rankFatJetPassed];
        } else {
        }
        if (rankFatJetPassed) {
          delete p4FatJet;
        }
      }
      if (vDeltaRDPairFatJet[iDPairMatching] < 0.8) {
        vGenDPairMatchedId[iDPairMatching] = true;
        // vDPairMatchedRankActual.emplace_back(
        vFatJetMatchedRankFatJetPassedPt.emplace_back(
            vDPairRankFatJetPassedPt[iDPairMatching]);
        nGenDPairMatched++;
      } else {
        vGenDPairMatchedId[iDPairMatching] = false;
      }
    }
    areAllGenDPairMatchingPassed = nGenDPairMatchingPassed == 2;
    areAllGenDPairMatched = nGenDPairMatched == 2;
    if (debug)
      std::printf(
          "areAllGenDPairMatchingPassed: %d, areAllGenDPairMatched: %d\n",
          areAllGenDPairMatchingPassed, areAllGenDPairMatched);
    if (nFatJetPassed >= 2) {
      std::sort(vFatJetMatchedRankFatJetPassedPt.begin(),
                vFatJetMatchedRankFatJetPassedPt.end());
      auto uniqueEnd = std::unique(vFatJetMatchedRankFatJetPassedPt.begin(),
                                   vFatJetMatchedRankFatJetPassedPt.end());
      nFatJetMatched =
          std::distance(vFatJetMatchedRankFatJetPassedPt.begin(), uniqueEnd);
      vFatJetMatchedRankFatJetPassedPt.resize(nFatJetMatched);
      if (nFatJetMatched) {
        areAllFatJetsMatchedLeading =
            vFatJetMatchedRankFatJetPassedPt.back() == nFatJetMatched - 1;
      } else {
        areAllFatJetsMatchedLeading = false;
      }
      if (debug) std::printf("nFatJetMatched: %d\n", nFatJetMatched);
    }

    lambdaFillTheTrees();
    // ttZMassCutted->Fill();
  }

  // ttNumCorrect->GetCurrentFile()->Write();
  // ttZMassCutted->GetCurrentFile()->Write();
  outFileTree->Write();

  // // TString outputFileHist = (TString) "output_" + outputFileHead + "_" +
  // //                          outputFileVar + "_" + outputFileTail +
  // "_hist"
  // + ".root";

  // // TFile *outFileHist = new TFile(outputFileHist.Data(),
  // toRecreateOutFile ? "recreate" : "update");

  // TString outImageDir =
  //     (TString) "../out_images/output_" + outputFileHead + "_" +
  //     outputFileTail;

  // // TCanvas* c1 = new TCanvas;
  // gStyle->SetOptStat(111111);
  // auto lambdaPrintHistograms = [/*&c1,*/ /*&outFileHist,*/
  // toRecreateOutFile, outImageDir, outputFileHead, outputFileVar,
  // outputFileTail](TTree *tt)-> void {
  //   std::cout << "Printing " << tt->GetName() << std::endl;
  //   TString nameTree = tt->GetName();
  //   TString outImageNameHead = (TString)outputFileHead + "_" +
  //   outputFileVar
  //   + "_" + nameTree; TString outImageCommonPath = outImageDir + "/" +
  //   outImageNameHead + "_"; TString outputFileHist = (TString) "output_" +
  //   outputFileHead + "_" +
  //                          outputFileVar + "_" + outputFileTail + "_hist_"
  //                          + nameTree + ".root";
  //   TFile *outFileHist = new TFile(outputFileHist.Data(), toRecreateOutFile
  //   ? "recreate" : "update"); for (TObject* leafObject :
  //   *(tt->GetListOfLeaves())) {
  //     TLeaf* leaf = (TLeaf*)leafObject;
  //     TString nameLeaf = leaf->GetName();
  //     if (nameLeaf.First("[") >= 0) {
  //       nameLeaf.Resize(nameLeaf.First("["));
  //     }
  //     TString nameHist = "h" + nameLeaf;
  //     tt->Draw(nameLeaf + ">>" + nameHist);
  //     TH1 *hist = (TH1 *) gDirectory->Get(nameHist);
  //     hist->SetTitle((TString)leaf->GetBranch()->GetTitle() + " (" +
  //     nameTree
  //     + ")"); TString outImagePath = outImageCommonPath + nameLeaf +
  //     ".svg";
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
