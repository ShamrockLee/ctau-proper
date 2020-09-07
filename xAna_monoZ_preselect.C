#include <TClonesArray.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>


#include <iostream>
#include <vector>
#include <cstring> // std::memcpy

#include "untuplizer.h"

// template <typename TypeData>
// TBranch* branchAdder(TTree* tree, const char* name, TypeData* address, const char* typestring=nullptr, Int_t buffsize=32000) {
//   if (typestring == nullptr) {
//     return tree->Branch(name, address);
//   } else {
//     return tree->Branch(name, address, ((TString)name + "/" + typestring), buffsize);
//   }
// }

typedef std::vector<Bool_t>::iterator IteratorVBool_t;
typedef std::vector<Int_t>::iterator IteratorVInt_t;
typedef std::vector<UInt_t>::iterator IteratorVUInt_t;
typedef std::vector<Float_t>::iterator IteratorVFloat_t;


void xAna_monoZ_preselect(std::string inputFile,
                          std::string outputFileHead,  // "jets-newdata"
                          std::string outputFileVar,   // "Mx2-150_Mx1-1_ctau-1"
                          std::string outputFileTail,  // "20200730"
                          bool toRecreateOutFile = true, bool debug = false) {

  const Double_t massZ = 91.1876;
  const Int_t pdgZ = 23;
  const Int_t pdgZp = 55;
  const Int_t pdgX2 = 18;
  const Int_t pdgX1 = 5000522;
  const Int_t pdgDown = 1;

  TString outputFileTree = (TString) "output_" + outputFileHead + "_" +
                       outputFileVar + "_" + outputFileTail + "_tree" + ".root";
  TFile* outFileTree =
      new TFile(outputFileTree.Data(), toRecreateOutFile ? "recreate" : "update");

  const char* nameTreeIn = "Events";
  TreeReader data(inputFile.data(), nameTreeIn);
  const Long64_t nEntry = data.GetEntriesFast();
  // if (debug) std::cout<< "nEntry: " << nEntry << std::endl;

  TTree* ttNumCorrect = new TTree("NumCorrect", "Variables of entries with 2 electrons/2 muons");
  TTree* ttPreselected = new TTree("Preselected", "Preselected variables");

  Long64_t jEntry;
  ttNumCorrect->Branch("jEntryOriginal", &jEntry, "jEntryOriginal/l", 64);
  ttPreselected->Branch("jEntryOriginal", &jEntry, "jEntryOriginal/l", 64);

  UInt_t nElectronPair;
  ttNumCorrect->Branch("nElectronPair", &nElectronPair, "nElectronPair/I", __SIZEOF_INT__);
  ttPreselected->Branch("nElectronPair", &nElectronPair, "nElectronPair/I", __SIZEOF_INT__);

  std::vector<Float_t> vElectronPair_massNumCorrect;
  vElectronPair_massNumCorrect.clear();
  ttNumCorrect->Branch("ElectronPair_mass", &vElectronPair_massNumCorrect, "ElectronPair_mass/F", __SIZEOF_FLOAT__);
  ttPreselected->Branch("ElectronPair_mass", &vElectronPair_massNumCorrect, "ElectronPair_mass/F", __SIZEOF_FLOAT__);

  UInt_t nMuonPair;
  ttNumCorrect->Branch("nMuonPair", &nMuonPair, "nMuonPair/I", __SIZEOF_INT__);
  ttPreselected->Branch("nMuonPair", &nMuonPair, "nMuonPair/I", __SIZEOF_INT__);

  std::vector<Float_t> vMuonPair_massNumCorrect;
  ttNumCorrect->Branch("MuonPair_mass", &vMuonPair_massNumCorrect, "MuonPair_mass/F", __SIZEOF_FLOAT__);
  ttPreselected->Branch("MuonPair_mass", &vMuonPair_massNumCorrect, "MuonPair_mass/F", __SIZEOF_FLOAT__);

  Float_t leptonPair_massNumCorrect;
  ttNumCorrect->Branch("LeptonPair_massNumCorrect", &leptonPair_massNumCorrect, "LeptonPair_massNumCorrect/F", __SIZEOF_FLOAT__);

  typedef struct {
    TString name;
    TString title;
  } ObjectDescription;

  UInt_t nJet;
  UInt_t nLeafOriginal;
  std::vector<ObjectDescription>
  vLeafDescriptionJetBool,
  vLeafDescriptionJetUInt,
  vLeafDescriptionJetInt,
  vLeafDescriptionJetFloat;
  {
    TFile* tfIn = new TFile(inputFile.data(), "READ");
    TTree* ttIn = (TTree*)tfIn->Get(nameTreeIn);
    TObjArray* tarrLeafOriginal = ttIn->GetListOfLeaves();

    nLeafOriginal = tarrLeafOriginal->GetSize();

    std::vector<Int_t> vIdxLeafJet;
    vIdxLeafJet.clear();
    for (UInt_t i = 0; i < nLeafOriginal; i++) {
      TLeaf* tlfCurrent = (TLeaf*)tarrLeafOriginal->At(i);
      ObjectDescription descriptionCurrent;
      descriptionCurrent.name = (TString)tlfCurrent->GetName();
      descriptionCurrent.title = (TString)tlfCurrent->GetTitle();
      TString typeNameCurrent = tlfCurrent->GetTypeName();
      if (descriptionCurrent.name.First("Jet_") == 0) {
        if (typeNameCurrent == "Bool_t") {
          vLeafDescriptionJetBool.push_back(descriptionCurrent);
        } else if (typeNameCurrent == "Int_t") {
          vLeafDescriptionJetInt.push_back(descriptionCurrent);
        } else if (typeNameCurrent == "UInt_t") {
          vLeafDescriptionJetUInt.push_back(descriptionCurrent);
        } else if (typeNameCurrent == "Float_t") {
          vLeafDescriptionJetFloat.push_back(descriptionCurrent);
        } else if (debug) {
          std::cout << (TString)"LeafJet missed: " + descriptionCurrent.name
          + "\t" + typeNameCurrent + "\t" + descriptionCurrent.title << std::endl;
        }
      }
    }
  }

  UInt_t nLeafJetBool = vLeafDescriptionJetBool.size();
  UInt_t nLeafJetInt = vLeafDescriptionJetInt.size();
  UInt_t nLeafJetUInt = vLeafDescriptionJetUInt.size();
  UInt_t nLeafJetFloat = vLeafDescriptionJetFloat.size();

  std::vector<Float_t> arrVLeafJetFloat[nLeafJetFloat],
  arrVLeafJetUInt[nLeafJetUInt],
  arrVLeafJetInt[nLeafJetInt],
  arrVLeafJetBool[nLeafJetBool];
  // Histograms for jets
  // TH1F *arrHLeafJetBool[nLeafJetBool],
  // *arrHLeafJetInt[nLeafJetInt],
  // *arrHLeafJetUInt[nLeafJetUInt],
  // *arrHLeafJetFloat[nLeafJetFloat];
  for (UInt_t i=0; i<nLeafJetFloat; i++) {
    TString name = vLeafDescriptionJetFloat[i].name;
    TString title = vLeafDescriptionJetFloat[i].title;
    ttPreselected->Branch(name, &(arrVLeafJetFloat[i]), name + "/F", __SIZEOF_FLOAT__ * 100)
    ->SetTitle(title);
    // arrHLeafJetFloat[i] = new TH1F((TString)"h" + vLeafDescriptionJetFloat[i].name,
    // vLeafDescriptionJetFloat[i].title,
    // 20000, -100, 1900);
  }
  for (UInt_t i=0; i<nLeafJetUInt; i++) {
    TString name = vLeafDescriptionJetUInt[i].name;
    TString title = vLeafDescriptionJetUInt[i].title;
    ttPreselected->Branch(name, &(arrVLeafJetUInt[i]), name + "/I", __SIZEOF_INT__ * 100)
    ->SetTitle(title);
    // arrHLeafJetUInt[i] = new TH1F((TString)"h" + vLeafDescriptionJetUInt[i].name,
    // vLeafDescriptionJetUInt[i].title,
    // 30, -0.5, 0.5);
  }
  for (UInt_t i=0; i<nLeafJetInt; i++) {
    TString name = vLeafDescriptionJetInt[i].name;
    TString title = vLeafDescriptionJetInt[i].title;
    ttPreselected->Branch(name, &(arrVLeafJetInt[i]), name + "/i", __SIZEOF_INT__ * 100)
    ->SetTitle(title);
    // arrHLeafJetInt[i] = new TH1F((TString)"h" + vLeafDescriptionJetInt[i].name,
    // vLeafDescriptionJetInt[i].title,
    // 30, -0.5, 30-0.5);
  }
  for (UInt_t i=0; i<nLeafJetBool; i++) {
    TString name = vLeafDescriptionJetBool[i].name;
    TString title = vLeafDescriptionJetBool[i].title;
    ttPreselected->Branch(name, &(arrVLeafJetBool[i]), name + "/O", __SIZEOF_INT__ * 100)
    ->SetTitle(title);
    // arrHLeafJetBool[i] = new TH1F((TString)"h" + vLeafDescriptionJetBool[i].name,
    // vLeafDescriptionJetBool[i].title,
    // 2, -0.5, 0.5);
  }

  for (jEntry = 0; jEntry < nEntry; jEntry++) {

    vElectronPair_massNumCorrect.clear();
    vMuonPair_massNumCorrect.clear();

    data.GetEntry(jEntry);

    nJet = data.GetInt("nJet");
    // if (debug) std::cout << "nJet: " << nJet << std::endl;

    // Int_t nGenPart = data.GetInt("nGenPart");
    // Int_t* ptrGenPart_pdgId = data.GetPtrInt("GenPart_pdgId");
    // Int_t* ptrGenPart_genPartIdxMother = data.GetPtrInt("GenPart_genPartIdxMother");

    Int_t nElectron = data.GetInt("nElectron");
    Bool_t* ptrElectron_mvaFall17V2Iso_WPL =
        data.GetPtrBool("Electron_mvaFall17V2Iso_WPL");
    // Bool_t* ptrElectron_mvaFall17V2Iso_WP90 =
    //     data.GetPtrBool("Electron_mvaFall17V2Iso_WP90"); // Probably too tight
    // Bool_t* ptrElectron_mvaFall17V2Iso_WP80 =
    //     data.GetPtrBool("Electron_mvaFall17V2Iso_WP80");
    Bool_t* ptrElectron_isSoft = ptrElectron_mvaFall17V2Iso_WPL;
    Bool_t* ptrElectron_isTight = ptrElectron_mvaFall17V2Iso_WPL;
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
    nElectronPair = isNumElectronCorrect;

    Int_t nMuonSoft, nMuonTight;
    {
      Int_t i;
      for (i = 0, nMuonSoft = 0, nMuonTight = 0; i < nElectron;
           nMuonSoft += ptrMuon_isSoft[i], nMuonTight += ptrMuon_isTight[i],
          i++)
        ;
    }
    Bool_t isNumMuonCorrect = (nMuonSoft == 2 && nMuonTight == 2);
    nMuonPair = isNumMuonCorrect;

    if (isNumElectronCorrect == isNumMuonCorrect) {
      continue;
    }
    if (debug) std::cout << "elepair, mupair: "
    << (Int_t)isNumElectronCorrect
    << (Int_t)isNumMuonCorrect << std::endl;

    Float_t* ptrElectron_pt = data.GetPtrFloat("Electron_pt");
    Float_t* ptrElectron_phi = data.GetPtrFloat("Electron_phi");
    Float_t* ptrElectron_eta = data.GetPtrFloat("Electron_eta");
    Float_t* ptrElectron_mass = data.GetPtrFloat("Electron_mass");

    TLorentzVector* ppElectronP4NumberCorrect[nElectronTight];
    TLorentzVector* ptrElectronP4NumberCorrectSum;
    if (isNumElectronCorrect) {
      ptrElectronP4NumberCorrectSum = new TLorentzVector();
      for (Int_t i = 0; i < nElectronTight; i++) {
        ppElectronP4NumberCorrect[i] = new TLorentzVector();
        ppElectronP4NumberCorrect[i]->SetPtEtaPhiM(
            ptrElectron_pt[i], ptrElectron_eta[i], ptrElectron_phi[i],
            ptrElectron_mass[i]);
        *ptrElectronP4NumberCorrectSum += *ppElectronP4NumberCorrect[i];
      }
      Double_t electronPair_mass = ptrElectronP4NumberCorrectSum->M();
      vElectronPair_massNumCorrect.push_back(electronPair_mass);
      leptonPair_massNumCorrect = electronPair_mass;
      Double_t massZCutUpper = massZ + 20;
      Double_t massZCutLower = massZ - 20;
      // Z mass cut for electron pairs
      if (massZCutLower > electronPair_mass || electronPair_mass > massZCutUpper) {
        ttNumCorrect->Fill();
        continue;
      }
      if (debug) std::cout << "massLPair: " << electronPair_mass << std::endl;
    }

    Float_t* ptrMuon_pt = data.GetPtrFloat("Muon_pt");
    Float_t* ptrMuon_phi = data.GetPtrFloat("Muon_phi");
    Float_t* ptrMuon_eta = data.GetPtrFloat("Muon_eta");
    Float_t* ptrMuon_mass = data.GetPtrFloat("Muon_mass");

    TLorentzVector* ppMuonP4NumberCorrect[nMuonTight];
    TLorentzVector* ptrMuonP4NumberCorrectSum;
    if (isNumMuonCorrect) {
      ptrMuonP4NumberCorrectSum = new TLorentzVector();
      for (Int_t i = 0; i < nMuonTight; i++) {
        ppMuonP4NumberCorrect[i] = new TLorentzVector();
        ppMuonP4NumberCorrect[i]->SetPtEtaPhiM(ptrMuon_pt[i], ptrMuon_eta[i],
                                               ptrMuon_phi[i], ptrMuon_mass[i]);
        *ptrMuonP4NumberCorrectSum += *ppMuonP4NumberCorrect[i];
      }
      Double_t muonPair_mass = ptrMuonP4NumberCorrectSum->M();
      vMuonPair_massNumCorrect.push_back(muonPair_mass);
      leptonPair_massNumCorrect = muonPair_mass;
      Double_t massZCutUpper = massZ + 20;
      Double_t massZCutLower = massZ - 20;
      // Z mass cut for muon pairs
      if (massZCutLower > muonPair_mass || muonPair_mass > massZCutUpper) {
        ttNumCorrect->Fill();
        continue;
      }
      if (debug) std::cout << "massLPair: " << muonPair_mass << std::endl;
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

    ttNumCorrect->Fill();
    // Jet_*
    for (UInt_t i=0; i<nLeafJetFloat; i++) {
      TString name = vLeafDescriptionJetFloat[i].name;
      Float_t *ptrFloat = data.GetPtrFloat(name);
      arrVLeafJetFloat[i].clear();
        for (UInt_t j = 0; j < nJet; j++) {
          // arrHLeafJetFloat[i]->Fill(ptrFloat[j]);
          arrVLeafJetFloat[i].push_back(ptrFloat[j]);
          // if (debug) std::cout << ptrFloat[j] << ", ";
        }
        // if (debug) std::cout << std::endl;
    }
    for (UInt_t i=0; i<nLeafJetUInt; i++) {
      TString name  = vLeafDescriptionJetUInt[i].name;
      Int_t *ptrUInt = data.GetPtrInt(name);
      arrVLeafJetUInt[i].clear();
      for (UInt_t j=0; j<nJet; j++) {
        // arrHLeafJetUInt[i]->Fill(ptrUInt[j]);
        arrVLeafJetUInt[i].push_back(ptrUInt[j]);
      }
    }
    for (UInt_t i=0; i<nLeafJetInt; i++) {
      TString name = vLeafDescriptionJetInt[i].name;
      Int_t *ptrInt = data.GetPtrInt(name);
      arrVLeafJetInt[i].clear();
      for (UInt_t j=0; j<nJet; j++) {
        // arrHLeafJetInt[i]->Fill(ptrInt[j]);
        arrVLeafJetInt[i].push_back(ptrInt[j]);
      }
    }
    for (UInt_t i=0; i<nLeafJetBool; i++) {
      TString name = vLeafDescriptionJetBool[i].name;
      Bool_t *ptrBool = data.GetPtrBool(name);
      for (UInt_t j=0; j<nJet; j++) {
        // arrHLeafJetBool[i]->Fill(ptrBool[j]);
        arrVLeafJetBool[i].push_back(ptrBool[j]);
      }
    }
    ttPreselected->Fill();
  }

  // ttNumCorrect->GetCurrentFile()->Write();
  ttPreselected->GetCurrentFile()->Write();

  TString outImageDir = (TString)"../out_images/output_" + outputFileHead + "_" + outputFileTail;
  TString outImageNameHead = (TString)outputFileHead + "_" + outputFileVar;
  TString outImageCommonPath = outImageDir + "/" + outImageNameHead + "_";
  TCanvas *c1 = new TCanvas;
  gStyle->SetOptStat(111111);

  for (TObject *leafObject: *(ttPreselected->GetListOfLeaves())) {
    TLeaf *leaf = (TLeaf *) leafObject;
    TString nameLeaf = leaf->GetName();
    if (nameLeaf.First("[") >=0 ) {
      nameLeaf.Resize(nameLeaf.First("["));
    }
    TString outImagePath = outImageCommonPath + nameLeaf + ".svg";
    c1->Clear();
    ttPreselected->Draw(nameLeaf);
    c1->Print(outImagePath);
  }
  c1->Close();

  outFileTree->Close();
  delete outFileTree;
  outFileTree = new TFile(outputFileTree.Data(), "update");
  for (UInt_t i=0; i<nLeafJetFloat; i++) {
    // arrHLeafJetFloat[i]->Write();
    // outFileTree->WriteObject(arrHLeafJetFloat[i], arrHLeafJetFloat[i]->GetName());
  }
  for (UInt_t i=0; i<nLeafJetUInt; i++) {
    // arrHLeafJetUInt[i]->Write();
    // outFileTree->WriteObject(arrHLeafJetUInt[i], arrHLeafJetUInt[i]->GetName());
  }
  for (UInt_t i=0; i<nLeafJetInt; i++) {
    // arrHLeafJetInt[i]->Write();
    // outFileTree->WriteObject(arrHLeafJetInt[i], arrHLeafJetInt[i]->GetName());
  }
  for (UInt_t i=0; i<nLeafJetBool; i++) {
    // arrHLeafJetBool[i]->Write();
    // outFileTree->WriteObject(arrHLeafJetBool[i], arrHLeafJetBool[i]->GetName());
  }

  outFileTree->Close();
}