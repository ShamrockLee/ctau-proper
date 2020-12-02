#include <Rtypes.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStyle.h>
#include <TVectorT.h>
#include <TVectorFfwd.h>

#include <cstring>  // std::memcpy
#include <iostream>
#include <vector>

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
    std::string inputFile,       //< Path of input file
    TString outputFileTree,
    // std::string outputFileHead,  //< Output filename head ex: "preselected"
    // std::string outputFileVar,   //< Variables of the data ex:
    //                              //"Mx2-150_Mv-500_Mx1-1_ctau-1"
    // std::string outputFileTail,  //< Output filename tail (date) ex: "20200730"
    bool toRecreateOutFile = true, bool debug = false) {
  const Double_t massZ = 91.1876;  //< static mass of Z (constant)
  const Double_t massElectron =
      0.0005109989461;          //< static mass of electron (constant)
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
  Bool_t arrGenHasLepton[3];
  for (Byte_t i = 0; i < 2; i++) {
    TString name = "Gen_has" + namesLepton[i] + "Pair";
    TString title = "Whether this is a Z-" + namesLeptonNota[i] +
                    namesLeptonNota[i] + " event";
    arrTTGen[i]
        ->Branch(name, &(arrGenHasLepton[i]), name + "/O", 1)
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
    arrTTGen[i]->Branch(name, &(arrIsPassingZMassCut[i]), name + "/O", 1)->SetTitle(title);
    arrTTNumCorrect[i]->Branch(name, &(arrIsPassingZMassCut[i]), name + "/O", 1)->SetTitle(title);
  }
  auto lambdaFillTheTrees = [&]() -> void {
    for (Byte_t i = 0; i < 3; i++) {
      if (arrGenHasLepton[i]) {
        arrTTGen[i]->Fill();
      }
    }
    for (Byte_t i = 0; i < 2; i++) {
      if (arrRecoIsLeptonNumCorrect[i]) {
        arrTTNumCorrect[i]->Fill();
      }
      if (arrIsPassingZMassCut[i]) {
        arrTTZMassCutted[i]->Fill();
      }
    }
  };

  /// Number of electron pairs (1 or 0), bounded to the trees
  // UInt_t nElectronPair;
  // ttNumCorrect->Branch("nElectronPair", &nElectronPair, "nElectronPair/I",
  //                      __SIZEOF_INT__);
  // ttZMassCutted->Branch("nElectronPair", &nElectronPair, "nElectronPair/I",
  //                       __SIZEOF_INT__);

  /// Electron pair mass, bounded to the trees
  // std::vector<Float_t> vElectronPair_massNumCorrect;
  // vElectronPair_massNumCorrect.clear();
  // ttNumCorrect->Branch("ElectronPair_mass", &vElectronPair_massNumCorrect);
  // ttZMassCutted->Branch("ElectronPair_mass", &vElectronPair_massNumCorrect);

  /// Number of muon pairs (1 or 0), bounded to the trees
  // UInt_t nMuonPair;
  // ttNumCorrect->Branch("nMuonPair", &nMuonPair, "nMuonPair/I",
  // __SIZEOF_INT__); ttZMassCutted->Branch("nMuonPair", &nMuonPair,
  // "nMuonPair/I", __SIZEOF_INT__);

  /// Muon pair mass, bounded to the treess
  // std::vector<Float_t> vMuonPair_massNumCorrect;
  // ttNumCorrect->Branch("MuonPair_mass", &vMuonPair_massNumCorrect);
  // ttZMassCutted->Branch("MuonPair_mass", &vMuonPair_massNumCorrect);

  // /// lepton pair mass, bounded to the trees
  // Float_t leptonPair_massNumCorrect;
  // ttNumCorrect->Branch("LeptonPair_massNumCorrect",
  // &leptonPair_massNumCorrect, "LeptonPair_massNumCorrect/F",
  // __SIZEOF_FLOAT__);

  typedef struct {
    TString name;
    TString title;
  } ObjectDescription;

  UInt_t nJet, nJetPassed;
  for (Byte_t i=0; i<3; i++) {
    TString name = "nJet";
    TString title = "Number of jets";
    arrTTGen[i]->Branch(name, &nJet, name + "/I", __SIZEOF_INT__)->SetTitle(title);
    arrTTGen[i]->Branch(name + "Passed", &nJetPassed, name + "Passed" + "/I", __SIZEOF_INT__)->SetTitle(title + " passing preselections");
  }
  for (Byte_t i=0; i<2; i++) {
    TString name = "nJet";
    TString title = "Number of jets";
    arrTTZMassCutted[i]->Branch(name, &nJet, name + "/I", __SIZEOF_INT__)->SetTitle(title);
    arrTTZMassCutted[i]->Branch(name + "Passed", &nJetPassed, name + "Passed" + "/I", __SIZEOF_INT__)->SetTitle(title + " passing preselections");
  }

  UInt_t nLeafOriginal;
  std::vector<ObjectDescription> vLeafDescriptionJetBool,
      vLeafDescriptionJetUInt, vLeafDescriptionJetInt, vLeafDescriptionJetFloat;
  std::vector<ObjectDescription> vLeafDescriptionHTFloat;
  {
    TFile* tfIn = TFile::Open(inputFile.data(), "READ");
    TTree* ttIn = (TTree*)tfIn->Get(nameTreeIn);
    TObjArray* tarrLeafOriginal = ttIn->GetListOfLeaves();

    nLeafOriginal = tarrLeafOriginal->GetSize();

    std::vector<Bool_t> arrVIdxLeafLeptonID[2];
    std::vector<Int_t> vIdxLeafJet;
    vIdxLeafJet.clear();
    for (UInt_t i = 0; i < nLeafOriginal; i++) {
      TLeaf* tlfCurrent = (TLeaf*)tarrLeafOriginal->At(i);
      ObjectDescription descriptionCurrent;
      descriptionCurrent.name = (TString)tlfCurrent->GetName();
      descriptionCurrent.title = (TString)tlfCurrent->GetBranch()->GetTitle();
      TString typeNameCurrent = tlfCurrent->GetTypeName();
      if (descriptionCurrent.name.SubString("Jet_").Start() == 0) {
        if (typeNameCurrent == "Bool_t") {
          vLeafDescriptionJetBool.push_back(descriptionCurrent);
        } else if (typeNameCurrent == "Int_t") {
          vLeafDescriptionJetInt.push_back(descriptionCurrent);
        } else if (typeNameCurrent == "UInt_t") {
          vLeafDescriptionJetUInt.push_back(descriptionCurrent);
        } else if (typeNameCurrent == "Float_t") {
          vLeafDescriptionJetFloat.push_back(descriptionCurrent);
        } else if (debug) {
          std::cout << (TString) "LeafJet missed: " + descriptionCurrent.name +
                           "\t" + typeNameCurrent + "\t" +
                           descriptionCurrent.title
                    << std::endl;
        }
      }
      if (descriptionCurrent.name.SubString("SoftActivityJetHT").Start() == 0) {
        if (typeNameCurrent == "Float_t") {
          if (debug) std::cout << "Found " << descriptionCurrent.name << " with type " << typeNameCurrent << std::endl;
          vLeafDescriptionHTFloat.push_back(descriptionCurrent);
        }
      }
    }
  }

  UInt_t nLeafJetBool = vLeafDescriptionJetBool.size();
  UInt_t nLeafJetInt = vLeafDescriptionJetInt.size();
  UInt_t nLeafJetUInt = vLeafDescriptionJetUInt.size();
  UInt_t nLeafJetFloat = vLeafDescriptionJetFloat.size();
  UInt_t nLeafHTFloat = vLeafDescriptionHTFloat.size();

  std::vector<Float_t> arrVLeafJetFloat[nLeafJetFloat],
      arrVLeafJetUInt[nLeafJetUInt], arrVLeafJetInt[nLeafJetInt],
      arrVLeafJetBool[nLeafJetBool];
  // Histograms for jets
  // TH1F *arrHLeafJetBool[nLeafJetBool],
  // *arrHLeafJetInt[nLeafJetInt],
  // *arrHLeafJetUInt[nLeafJetUInt],
  // *arrHLeafJetFloat[nLeafJetFloat];
  for (UInt_t i = 0; i < nLeafJetFloat; i++) {
    TString name = vLeafDescriptionJetFloat[i].name;
    TString title = vLeafDescriptionJetFloat[i].title;
    // ttZMassCutted->Branch(name, &(arrVLeafJetFloat[i]))->SetTitle(title);
    for (Byte_t j = 0; j < 3; j++) {
      arrTTGen[j]->Branch(name, &(arrVLeafJetFloat[i]))->SetTitle(title);
    }
    for (Byte_t j = 0; j < 2; j++) {
      arrTTZMassCutted[j]
          ->Branch(name, &(arrVLeafJetFloat[i]))
          ->SetTitle(title);
    }
    // arrHLeafJetFloat[i] = new TH1F((TString)"h" +
    // vLeafDescriptionJetFloat[i].name, vLeafDescriptionJetFloat[i].title,
    // 20000, -100, 1900);
  }
  for (UInt_t i = 0; i < nLeafJetUInt; i++) {
    TString name = vLeafDescriptionJetUInt[i].name;
    TString title = vLeafDescriptionJetUInt[i].title;
    // ttZMassCutted->Branch(name, &(arrVLeafJetUInt[i]))->SetTitle(title);
    for (Byte_t j = 0; j < 3; j++) {
      arrTTGen[j]->Branch(name, &(arrVLeafJetUInt[i]))->SetTitle(title);
    }
    for (Byte_t j = 0; j < 2; j++) {
      arrTTZMassCutted[j]->Branch(name, &(arrVLeafJetUInt[i]))->SetTitle(title);
    }
    // arrHLeafJetUInt[i] = new TH1F((TString)"h" +
    // vLeafDescriptionJetUInt[i].name, vLeafDescriptionJetUInt[i].title, 30,
    // -0.5, 0.5);
  }
  for (UInt_t i = 0; i < nLeafJetInt; i++) {
    TString name = vLeafDescriptionJetInt[i].name;
    TString title = vLeafDescriptionJetInt[i].title;
    // ttZMassCutted->Branch(name, &(arrVLeafJetInt[i]))->SetTitle(title);
    for (Byte_t j = 0; j < 3; j++) {
      arrTTGen[j]->Branch(name, &(arrVLeafJetInt[i]))->SetTitle(title);
    }
    for (Byte_t j = 0; j < 2; j++) {
      arrTTZMassCutted[j]->Branch(name, &(arrVLeafJetInt[i]))->SetTitle(title);
    }
    // arrHLeafJetInt[i] = new TH1F((TString)"h" +
    // vLeafDescriptionJetInt[i].name, vLeafDescriptionJetInt[i].title, 30,
    // -0.5, 30-0.5);
  }
  for (UInt_t i = 0; i < nLeafJetBool; i++) {
    TString name = vLeafDescriptionJetBool[i].name;
    TString title = vLeafDescriptionJetBool[i].title;
    // ttZMassCutted->Branch(name, &(arrVLeafJetBool[i]))->SetTitle(title);
    for (Byte_t j = 0; j < 3; j++) {
      arrTTGen[j]->Branch(name, &(arrVLeafJetBool[i]))->SetTitle(title);
    }
    for (Byte_t j = 0; j < 2; j++) {
      arrTTZMassCutted[j]->Branch(name, &(arrVLeafJetBool[i]))->SetTitle(title);
    }
    // arrHLeafJetBool[i] = new TH1F((TString)"h" +
    // vLeafDescriptionJetBool[i].name, vLeafDescriptionJetBool[i].title, 2,
    // -0.5, 0.5);
  }

  Float_t arrHTFloat[nLeafHTFloat];
  for(UInt_t i=0; i<nLeafHTFloat; i++) {
    TString name = vLeafDescriptionHTFloat[i].name;
    TString title = vLeafDescriptionHTFloat[i].title;
    for (Byte_t j=0; j<3; j++) {
      arrTTGen[j]->Branch(name, &(arrHTFloat[i]), __SIZEOF_FLOAT__)->SetTitle(title);
    }
    for (Byte_t j=0; j<2; j++) {
      arrTTZMassCutted[j]->Branch(name, &(arrHTFloat[i]), __SIZEOF_FLOAT__);
    }
  }

  for (jEntry = 0; jEntry < nEntry; jEntry++) {
    data.GetEntry(jEntry);

    // vElectronPair_massNumCorrect.clear();
    // vMuonPair_massNumCorrect.clear();

    for (UInt_t i = 0; i < nLeafJetFloat; i++) {
      arrVLeafJetFloat[i].clear();
    }
    for (UInt_t i = 0; i < nLeafJetUInt; i++) {
      arrVLeafJetUInt[i].clear();
    }
    for (UInt_t i = 0; i < nLeafJetInt; i++) {
      arrVLeafJetInt[i].clear();
    }
    for (UInt_t i = 0; i < nLeafJetBool; i++) {
      arrVLeafJetBool[i].clear();
    }

    nJet = data.GetInt("nJet");
    // if (debug) std::cout << "nJet: " << nJet << std::endl;

    Bool_t Gen_isRightEvent = false;
    UInt_t nGenPart = data.GetInt("nGenPart");
    Int_t* ptrGenPart_pdgId = data.GetPtrInt("GenPart_pdgId");
    Int_t* ptrGenPart_genPartIdxMother =
    data.GetPtrInt("GenPart_genPartIdxMother");
    for (Byte_t i=0; i<3; i++) arrGenHasLepton[i] = false;
    for (UInt_t ig=0; ig<nGenPart; ig++) {
      if (ptrGenPart_pdgId[ptrGenPart_genPartIdxMother[ig]] == pdgZ) {
        Int_t genparIdAbs = TMath::Abs(ptrGenPart_pdgId[ig]);
        if (genparIdAbs == pdgElectron) {
          arrGenHasLepton[0] = true;
          Gen_isRightEvent = true;
          // break;
        } else if (genparIdAbs == pdgMuon) {
          arrGenHasLepton[1] = true;
          Gen_isRightEvent = true;
          // break;
        } else if (genparIdAbs == pdgTau) {
          arrGenHasLepton[2] = true;
          Gen_isRightEvent = true;
          // break;
        }
      }
    }

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

    for (Byte_t i=0; i<2; i++) arrIsPassingZMassCut[i] = 0;

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

    if (!arrIsPassingZMassCut[0] && !arrIsPassingZMassCut[1] && !Gen_isRightEvent) {
      lambdaFillTheTrees();
      continue; //TODO
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
    std::vector<UInt_t> vJJetPassed(nJet);
    vJJetPassed.clear();
    Float_t *ptrJet_pt_original = data.GetPtrFloat("Jet_pt");
    Float_t *ptrJet_eta_original = data.GetPtrFloat("Jet_eta");
    for (UInt_t jJet=0; jJet<nJet; jJet++) {
      if (ptrJet_pt_original[jJet] < 30) continue;
      if (TMath::Abs(ptrJet_eta_original[jJet]) > 3) continue;
      vJJetPassed.push_back(jJet);
    }
    nJetPassed = vJJetPassed.size();
    for (UInt_t i = 0; i < nLeafJetFloat; i++) {
      TString name = vLeafDescriptionJetFloat[i].name;
      Float_t* ptrFloat = data.GetPtrFloat(name);
      arrVLeafJetFloat[i].clear();
      for (UInt_t j: vJJetPassed) {
        // arrHLeafJetFloat[i]->Fill(ptrFloat[j]);
        arrVLeafJetFloat[i].push_back(ptrFloat[j]);
        // if (debug) std::cout << ptrFloat[j] << ", ";
      }
      // if (debug) std::cout << std::endl;
    }
    for (UInt_t i = 0; i < nLeafJetUInt; i++) {
      TString name = vLeafDescriptionJetUInt[i].name;
      Int_t* ptrUInt = data.GetPtrInt(name);
      arrVLeafJetUInt[i].clear();
      for (UInt_t j: vJJetPassed) {
        // arrHLeafJetUInt[i]->Fill(ptrUInt[j]);
        arrVLeafJetUInt[i].push_back(ptrUInt[j]);
      }
    }
    for (UInt_t i = 0; i < nLeafJetInt; i++) {
      TString name = vLeafDescriptionJetInt[i].name;
      Int_t* ptrInt = data.GetPtrInt(name);
      arrVLeafJetInt[i].clear();
      for (UInt_t j: vJJetPassed) {
        // arrHLeafJetInt[i]->Fill(ptrInt[j]);
        arrVLeafJetInt[i].push_back(ptrInt[j]);
      }
    }
    for (UInt_t i = 0; i < nLeafJetBool; i++) {
      TString name = vLeafDescriptionJetBool[i].name;
      Bool_t* ptrBool = data.GetPtrBool(name);
      for (UInt_t j: vJJetPassed) {
        // arrHLeafJetBool[i]->Fill(ptrBool[j]);
        arrVLeafJetBool[i].push_back(ptrBool[j]);
      }
    }
    for (UInt_t i=0; i<nLeafHTFloat; i++) {
      TString name = vLeafDescriptionHTFloat[i].name;
      arrHTFloat[i] = data.GetFloat(name);
    }

    lambdaFillTheTrees();
    // ttZMassCutted->Fill();
  }

  // ttNumCorrect->GetCurrentFile()->Write();
  // ttZMassCutted->GetCurrentFile()->Write();
  outFileTree->Write();

  
  // // TString outputFileHist = (TString) "output_" + outputFileHead + "_" +
  // //                          outputFileVar + "_" + outputFileTail + "_hist" + ".root";

  // // TFile *outFileHist = new TFile(outputFileHist.Data(), toRecreateOutFile ? "recreate" : "update");

  // TString outImageDir =
  //     (TString) "../out_images/output_" + outputFileHead + "_" + outputFileTail;
  
  // // TCanvas* c1 = new TCanvas;
  // gStyle->SetOptStat(111111);
  // auto lambdaPrintHistograms = [/*&c1,*/ /*&outFileHist,*/ toRecreateOutFile, outImageDir, outputFileHead, outputFileVar, outputFileTail](TTree *tt)-> void {
  //   std::cout << "Printing " << tt->GetName() << std::endl;
  //   TString nameTree = tt->GetName();
  //   TString outImageNameHead = (TString)outputFileHead + "_" + outputFileVar + "_" + nameTree;
  //   TString outImageCommonPath = outImageDir + "/" + outImageNameHead + "_";
  //   TString outputFileHist = (TString) "output_" + outputFileHead + "_" +
  //                          outputFileVar + "_" + outputFileTail + "_hist_" + nameTree + ".root";
  //   TFile *outFileHist = new TFile(outputFileHist.Data(), toRecreateOutFile ? "recreate" : "update");
  //   for (TObject* leafObject : *(tt->GetListOfLeaves())) {
  //     TLeaf* leaf = (TLeaf*)leafObject;
  //     TString nameLeaf = leaf->GetName();
  //     if (nameLeaf.First("[") >= 0) {
  //       nameLeaf.Resize(nameLeaf.First("["));
  //     }
  //     TString nameHist = "h" + nameLeaf;
  //     tt->Draw(nameLeaf + ">>" + nameHist);
  //     TH1 *hist = (TH1 *) gDirectory->Get(nameHist);
  //     hist->SetTitle((TString)leaf->GetBranch()->GetTitle() + " (" + nameTree + ")");
  //     TString outImagePath = outImageCommonPath + nameLeaf + ".svg";
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
