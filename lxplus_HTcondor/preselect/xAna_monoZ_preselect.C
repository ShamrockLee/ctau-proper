#include <Rtypes.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TMath.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

template<typename V = ROOT::RDFDetail::RInferredType>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumn(ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> &df, const TString nameColumn) {
  const TString typenameColumn = df.GetColumnType(nameColumn.Data());
  Int_t nBins;
  Double_t lowerLimit, upperLimit;
  if (typenameColumn.Contains("Bool") || typenameColumn.Contains("bool")) {
    nBins = 2;
    lowerLimit = 0.;
    upperLimit = 2.;
  } else if (typenameColumn.Contains("int") || typenameColumn.Contains("Int")) {
    nBins = 20001;
    lowerLimit = - static_cast<Double_t>(nBins / 2) - 0.5;
    upperLimit = lowerLimit + nBins;
  // } else if (typenameColumn.Contains("float") || typenameColumn.Contains("Float") || typenameColumn.Contains("double") || typenameColumn.Contains("Double")) {
  } else {
    nBins = 20000;
    lowerLimit = -static_cast<Double_t>(nBins) / 2;
    upperLimit = lowerLimit + nBins;
  }
  return df.Histo1D<V>({"h" + nameColumn, nameColumn, nBins, lowerLimit, upperLimit}, nameColumn.Data());
}

void xAna_monoZ_preselect(const TString fileIn, const TString fileOut, const Bool_t recreate=true, const Bool_t debug=false) {
  ROOT::EnableImplicitMT();
  constexpr Double_t massZ = 91.1876;  //< static mass of Z (constant)
  constexpr Double_t massElectron =
      0.0005109989461;  //< static mass of electron (constant)
  constexpr Double_t massDown = 0.0048;
  constexpr Int_t pdgZ = 23;        //< pdgid of Z
  constexpr Int_t pdgZp = 55;       //< pdgid of Z'
  constexpr Int_t pdgX2 = 18;       //< pdgid of x2
  constexpr Int_t pdgX1 = 5000522;  //< pdgid of x1
  constexpr Int_t pdgDown = 1;      //< pdgid of down quark
  constexpr Int_t pdgElectron = 11;
  constexpr Int_t pdgMuon = 13;
  constexpr Int_t pdgTau = 15;
  const Bool_t isSignal = fileOut.Contains("ignal");
  // const TString namesLepton[] = {"Electron", "Muon",
  //                                "Tau"};  //< the name of the leptons (Xxx)
  // const TString namesLeptonLower[] = {"electron", "muon",
  //                                     "tau"};  //< the name of the leptons (xxx)
  // const TString namesLeptonNota[] = {"e", "mu", "tau"};
  
  const ROOT::RDataFrame dfIn("tree/treeMaker", {fileIn});
  ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> &&dfOriginalTemp = dfIn;
  for (const TString prefNameCol: {"ele", "mu", "THINjet", "FATjet"}) {
    const TString nameColP4 = prefNameCol + "P4";
    dfOriginalTemp = dfOriginalTemp
    .Define(prefNameCol + "Pt", [](const TClonesArray tcaP4)->std::vector<Float_t>{
      std::vector<Float_t> vResult(tcaP4.GetSize(), 0.);
      for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
        vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->Pt();
      }
      return vResult;
    }, {nameColP4})
    .Define(prefNameCol + "Eta", [](const TClonesArray tcaP4)->std::vector<Float_t>{
      std::vector<Float_t> vResult(tcaP4.GetSize(), 0.);
      for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
        vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->Eta();
      }
      return vResult;
    }, {nameColP4})
    .Define(prefNameCol + "Phi", [](const TClonesArray tcaP4)->std::vector<Float_t>{
      std::vector<Float_t> vResult(tcaP4.GetSize(), 0.);
      for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
        vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->Phi();
      }
      return vResult;
    }, {nameColP4})
    .Define(prefNameCol + "M", [](const TClonesArray tcaP4)->std::vector<Float_t>{
      std::vector<Float_t> vResult(tcaP4.GetSize(), 0.);
      for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
        vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->M();
      }
      return vResult;
    }, {nameColP4});
  }
  dfOriginalTemp = dfOriginalTemp
  .Define("THINjetIdxPassed", [](const Int_t THINnJet, const std::vector<Float_t> THINjetPt, const std::vector<Float_t> THINjetEta)->std::vector<Int_t>{
    std::vector<Int_t> THINjetIdxPassed(THINnJet);
    THINjetIdxPassed.clear();
    for (Int_t iJet = 0; iJet < THINnJet; ++iJet) {
      if (THINjetPt[iJet] > 30. && TMath::Abs(THINjetEta[iJet]) < 2.5) {
        THINjetIdxPassed.emplace_back(iJet);
      }
    }
    return THINjetIdxPassed;
  }, {"THINnJet", "THINjetPt", "THINjetEta"});
  auto dfOriginal = dfOriginalTemp
  std::vector<ROOT::RDF::RResultPtr<TH1D>> vHistViewOriginal(200);
  vHistViewOriginal.clear();
  vHistViewOriginal.emplace_back(GetHistFromColumn<Int_t>(dfOriginal, "THINnJet"));
  vHistViewOriginal.emplace_back(GetHistFromColumn(dfOriginal, "THINjetPt"));
  vHistViewOriginal.emplace_back(GetHistFromColumn<std::vector<Int_t>>(dfOriginal, "THINjetIdxPassed"));
  TFile* tfOut = TFile::Open(fileOut, recreate ? "recreate" : "update");
  tfOut->mkdir("Original", "The distributions of unfiltered entries");
  tfOut->cd("Original");
  for (auto &&histView: vHistViewOriginal) {
    histView->Write();
  }
  tfOut->Close();
}
