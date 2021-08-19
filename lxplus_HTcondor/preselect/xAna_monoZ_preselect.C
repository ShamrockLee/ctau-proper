#include <Rtypes.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <ROOT/RVec.hxx>
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

#include <string>
#include <string_view>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <regex>

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

Bool_t RefgetE2D(std::string &typenameE, const std::string typenameCol) {
  std::regex r ("^ROOT::VecOps::RVec<vector<(.+)>>$");
  std::smatch m;
  const Bool_t result = (std::regex_match(typenameCol, m, r));
  if (result) {
    typenameE = m[1];
  }
  return result;
}

Bool_t RefgetE1D(std::string &typenameE, const std::string typenameCol) {
  std::regex r ("^ROOT::VecOps::RVec<(.+)>$");
  std::smatch m;
  const Bool_t result = (std::regex_match(typenameCol, m, r));
  if (result) {
    typenameE = m[1];
  }
  return result;
}
void xAna_monoZ_preselect(const std::string fileIn, const std::string fileOut, const Bool_t recreate=true, const Bool_t debug=false) {
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
  const Bool_t isSignal = TString(fileOut).Contains("ignal");
  // const TString namesLepton[] = {"Electron", "Muon",
  //                                "Tau"};  //< the name of the leptons (Xxx)
  // const TString namesLeptonLower[] = {"electron", "muon",
  //                                     "tau"};  //< the name of the leptons (xxx)
  // const TString namesLeptonNota[] = {"e", "mu", "tau"};
  TFile *tfIn = TFile::Open(fileIn.c_str());
  TTree *ttIn = tfIn->Get<TTree>("tree/treeMaker");
  ROOT::RDataFrame dfIn(*ttIn);
  std::vector<std::string>
    vNameColOriginal = {};
  std::vector<std::vector<std::string>> vvNameColPreselectedTHINjet = { {}, {} };
  ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> &&dfOriginalTemp = std::move(dfIn);
  {
    auto dfOriginalCurrent = dfOriginalTemp;
    auto vNamesCol = dfOriginalCurrent.GetColumnNames();
    std::vector<std::string> vTypenamesCol(vNamesCol.size(), "");
    for (size_t iCol = 0; iCol < vNamesCol.size(); ++iCol) {
      vTypenamesCol[iCol] = dfOriginalCurrent.GetColumnType(vNamesCol[iCol]);
    }
    dfOriginalTemp = std::move(dfOriginalCurrent);
    for (const std::string nameCol: vNamesCol) {
      const TString tstrNameCol(nameCol);
      // if (debug) std::cerr << tstrNameCol << "\tLength: " << tstrNameCol.Length() << std::endl;
      if (tstrNameCol.EndsWith("P4")) {
        const TString tstrPrefNameCol = tstrNameCol(0, tstrNameCol.Length() - 2);
        const std::string prefNameCol (tstrPrefNameCol.Data());
        if (debug) std::cerr << "prefNameCol: " << prefNameCol << std::endl;
        dfOriginalTemp = dfOriginalTemp
        .Define(prefNameCol + "Pt", [](const TClonesArray tcaP4) {
          ROOT::RVec<Float_t> vResult(tcaP4.GetSize(), 0.);
          for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
            vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->Pt();
          }
          return vResult;
        }, {nameCol})
        .Define(prefNameCol + "Eta", [](const TClonesArray tcaP4) {
          ROOT::RVec<Float_t> vResult(tcaP4.GetSize(), 0.);
          for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
            vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->Eta();
          }
          return vResult;
        }, {nameCol})
        .Define(prefNameCol + "Phi", [](const TClonesArray tcaP4) {
          ROOT::RVec<Float_t> vResult(tcaP4.GetSize(), 0.);
          for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
            vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->Phi();
          }
          return vResult;
        }, {nameCol})
        .Define(prefNameCol + "E", [](const TClonesArray tcaP4) {
          ROOT::RVec<Float_t> vResult(tcaP4.GetSize(), 0.);
          for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
            vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->E();
          }
          return vResult;
        }, {nameCol})
        .Define(prefNameCol + "M", [](const TClonesArray tcaP4) {
          ROOT::RVec<Float_t> vResult(tcaP4.GetSize(), 0.);
          for (Int_t iJet = 0; iJet < tcaP4.GetSize(); ++iJet) {
            vResult[iJet] = static_cast<TLorentzVector *>(tcaP4[iJet])->M();
          }
          return vResult;
        }, {nameCol});
        continue;
      }
    }
  }
  dfOriginalTemp = dfOriginalTemp
  .Define("eleIdxPassPtEta", [](const Int_t nEle, const ROOT::RVec<Float_t> elePt, const ROOT::RVec<Float_t> eleEta){
    ROOT::RVec<size_t> vResult(nEle);
    vResult.clear();
    for (Int_t iEle = 0; iEle < nEle; ++iEle) {
      if (elePt[iEle] > 20 && TMath::Abs(eleEta[iEle]) < 2.5) {
        vResult.emplace_back(iEle);
      }
    }
    return vResult;
  }, {"nEle", "elePt", "eleEta"})
  .Define("nElePassPtEta", "static_cast<Int_t>(eleIdxPassPtEta.size())")
  .Redefine("nEle", "nElePassPtEta")
  .Define("muIdxPassPtEta", [](const Int_t nMu, const ROOT::RVec<Float_t> muPt, const ROOT::RVec<Float_t> muEta) {
    ROOT::RVec<size_t> vResult(nMu);
    vResult.clear();
    for (Int_t iMu = 0; iMu < nMu; ++iMu) {
      if (muPt[iMu] > 20 && TMath::Abs(muEta[iMu]) < 2.5) {
        vResult.emplace_back(iMu);
      }
    }
    return vResult;
  }, {"nMu", "muPt", "muEta"})
  .Define("nMuPassPtEta", "static_cast<Int_t>(muIdxPassPtEta.size())")
  .Redefine("nMu", "nMuPassPtEta")
  .Define("THINjetIdxPassPtEta", [](const Int_t THINnJet, const ROOT::RVec<Float_t> THINjetPt, const ROOT::RVec<Float_t> THINjetEta) {
    ROOT::RVec<size_t> vResult(THINnJet);
    vResult.clear();
    for (Int_t iJet = 0; iJet < THINnJet; ++iJet) {
      if (THINjetPt[iJet] > 30. && TMath::Abs(THINjetEta[iJet]) < 2.5) {
        vResult.emplace_back(iJet);
      }
    }
    return vResult;
  }, {"THINnJet", "THINjetPt", "THINjetEta"})
  .Define("THINnJetPassPtEta", "static_cast<Int_t>(THINjetIdxPassPtEta.size())")
  .Redefine("THINnJet", "THINnJetPassPtEta")
  .Define("FATjetIdxPassPtEta", [](const Int_t FATnJet, const ROOT::RVec<Float_t> FATjetPt, const ROOT::RVec<Float_t> FATjetEta) {
    ROOT::RVec<size_t> vResult(FATnJet);
    vResult.clear();
    for (Int_t iJet = 0; iJet < FATnJet; ++iJet) {
      if (FATjetPt[iJet] > 200. && TMath::Abs(FATjetEta[iJet]) < 2.5) {
        vResult.emplace_back(iJet);
      }
    }
    return vResult;
  }, {"FATnJet", "FATjetPt", "FATjetEta"})
  .Define("FATnJetPassPtEta", "static_cast<Int_t>(FATjetIdxPassPtEta.size())")
  .Redefine("FATnJet", "FATnJetPassPtEta");
  for (std::string nameCol: {"nEle", "nMu", "THINnJet", "FATnJet"}) {
    vNameColOriginal.emplace_back(nameCol + "PassPtEta");
  }
  {
    auto dfOriginalCurrent = dfOriginalTemp;
    auto vNamesCol = dfOriginalCurrent.GetColumnNames();
    std::vector<std::string> vTypenamesCol(vNamesCol.size());
    vTypenamesCol.clear();
    for (const std::string nameCol: vNamesCol) {
      vTypenamesCol.emplace_back(dfOriginalCurrent.GetColumnType(nameCol));
    }
    dfOriginalTemp = std::move(dfOriginalCurrent);
    for (std::string &nameCol: vNamesCol) {
      std::regex r("^is(.+)Muon$");
      std::smatch m;
      if (std::regex_match(nameCol, m, r)) {
        const std::string nameAlias = "muIsPass" + std::string(m[1]);
        if (debug) std::cerr << "nameAlias: " << nameAlias << std::endl;
        dfOriginalTemp = dfOriginalTemp.Define(nameAlias, nameCol);
        nameCol = nameAlias;
      }
    }
    for (size_t iCol = 0; iCol < vNamesCol.size(); ++iCol) {
      const std::string &nameCol = vNamesCol[iCol];
      const std::string &typenameCol = vTypenamesCol[iCol];
      const TString tstrNameCol = nameCol;
      const TString tstrTypenameCol = typenameCol;
      for (const std::string pref: {"ele", "mu", "THINjet", "FATjet"}) {
        if (tstrNameCol.BeginsWith(pref) &&
          ! tstrNameCol.EndsWith("PassPtEta") && (
          tstrTypenameCol.Contains("vector") ||
          tstrTypenameCol.Contains("RVec") ||
          tstrTypenameCol.Contains("array")
        ) && (
          tstrTypenameCol.Contains("ouble") ||
          tstrTypenameCol.Contains("loat") ||
          tstrTypenameCol.Contains("Long") ||
          tstrTypenameCol.Contains("long") ||
          tstrTypenameCol.Contains("Int") ||
          tstrTypenameCol.Contains("int")||
          tstrTypenameCol.Contains("Bool") ||
          tstrTypenameCol.Contains("bool")
        ) && ! (
          tstrTypenameCol.Contains("TObjArray") ||
          tstrTypenameCol.Contains("TClonesArray")
        )) {
          if (debug) std::cerr << nameCol << ", " << typenameCol << std::endl;
          // vNameColOriginal.emplace_back(nameCol);
          std::string typenameE = "";
          if (RefgetE2D(typenameE, typenameCol)) {
            if (debug) std::cerr << "Element type of the 2D vector: " << typenameE << std::endl;
            // // Doesn't work
            // dfOriginalTemp = dfOriginalTemp
            // .Redefine(nameCol, Form("ROOT::VecOps::Map(ROOT::VecOps::Take(%s, %sIdxPassPtEta), [](const std::vector<%s> &vE){ROOT::VecOps::RVec<%s> result(vE); return result;})", nameCol.c_str(), pref.c_str(), typenameE.c_str(), typenameE.c_str()));
            dfOriginalTemp = dfOriginalTemp
            .Redefine(nameCol,
              "ROOT::VecOps::RVec<ROOT::VecOps::RVec<" + typenameE + ">> vvResult(" + nameCol + ".size());"
              "vvResult.clear();"
              "for (size_t iVE: " + pref + "IdxPassPtEta) {"
              "  const auto& vE = " + nameCol + "[iVE];"
              "  ROOT::VecOps::RVec<" + typenameE + "> resultSub(vE.size());"
              "  resultSub.clear();"
              "  for (auto iterElement = vE.cbegin(); iterElement != vE.cend(); ++iterElement) {"
              "    resultSub.push_back(*iterElement);"
              "  } "
              "  vvResult.push_back(resultSub);"
              "} "
              "return vvResult");
          } else {
            if (!(tstrTypenameCol.Contains("Bool") || tstrTypenameCol.Contains("bool"))){
              vNameColOriginal.emplace_back(nameCol);
            }
            // See ROOT issue 8857 and pr 8859
            if (debug) std::cerr << "Redefining .." << std::endl;
            dfOriginalTemp = dfOriginalTemp
            .Redefine(nameCol, Form("ROOT::VecOps::Take(%s, %sIdxPassPtEta)", nameCol.c_str(), pref.c_str()));
          }
          break;
        }
      }
    }
  }
  // vNameColOriginal.emplace_back("eleIsPassLoose");
  dfOriginalTemp = dfOriginalTemp
  .Define("leptonPairFlavor", [](const ROOT::RVec<Bool_t> eleIsLoose, const ROOT::RVec<Bool_t> eleIsMedium, const ROOT::RVec<Bool_t> eleIsTight, const ROOT::RVec<Bool_t> muIsLoose, const ROOT::RVec<Bool_t> muIsMedium, const ROOT::RVec<Bool_t> muIsTight)->Int_t{
    const size_t nEleLoose = ROOT::VecOps::Sum(eleIsLoose);
    const size_t nEleMedium = ROOT::VecOps::Sum(eleIsMedium);
    const size_t nEleTight = ROOT::VecOps::Sum(eleIsTight);
    const size_t nMuLoose = ROOT::VecOps::Sum(muIsLoose);
    const size_t nMuMedium = ROOT::VecOps::Sum(muIsMedium);
    const size_t nMuTight = ROOT::VecOps::Sum(muIsTight);
    if (nEleTight >= 4 || nMuTight >= 4) {
      return -1;
    }
    const Bool_t isEleNumCorrect = nEleMedium >= 1 && nEleLoose >= 2;
    const Bool_t isMuNumCorrect = nMuMedium >= 1 && nMuLoose >= 2;
    if (isEleNumCorrect && isMuNumCorrect) {
      return -1;
    } else if (isEleNumCorrect) {
      return 0;
    } else if (isMuNumCorrect) {
      return 1;
    }
    return -1;
  },{"eleIsPassLoose", "eleIsPassMedium", "eleIsPassTight", "muIsPassLoose", "muIsPassMedium", "muIsPassTight"});
  auto dfOriginal = dfOriginalTemp;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> vHistViewOriginal(vNameColOriginal.size());
  vHistViewOriginal.clear();
  for (const TString &nameCol: vNameColOriginal) {
    vHistViewOriginal.emplace_back(GetHistFromColumn(dfOriginal, nameCol));
  }
  std::array<ROOT::RDF::RNode, 2> aDfNumCorrect = {dfOriginal, dfOriginal};
  {
    const std::array<std::string, 2> aPref{"Ele", "Mu"};
    const std::array<std::string, 2> aPrefLower{"ele", "mu"};
    for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
      aDfNumCorrect[iLepFlav] = aDfNumCorrect[iLepFlav]
      .Filter(Form("leptonPairFlavor==%zu", iLepFlav))
      .Define("idxLeptonNumCorrect", [](const Int_t nLep, const ROOT::RVec<Bool_t> isLoose, const ROOT::RVec<Bool_t> isMedium, const ROOT::RVec<Bool_t> isTight)->ROOT::RVec<Int_t>{
        ROOT::RVec<Int_t> vResult(2, -1);
        if (nLep < 2) {
          return vResult;
        }
        for (vResult[0] = 0; vResult[0] < nLep - 1; ++(vResult[0])) {
          for (vResult[1] = vResult[0] + 1; vResult[1] < nLep; ++(vResult[1])) {
            if (ROOT::VecOps::All(ROOT::VecOps::Take(isLoose, vResult)) && ROOT::VecOps::Sum(ROOT::VecOps::Take(isMedium, vResult)) >= 1) {
              return vResult;
            }
          }
        }
        std::fill(vResult.begin(), vResult.end(), -1);
        return vResult;
      }, {"n" + aPref[iLepFlav], aPrefLower[iLepFlav] + "IsPassLoose", aPrefLower[iLepFlav] + "IsPassMedium", aPrefLower[iLepFlav] + "IsPassTight"});
    }
  }
  TFile* tfOut = TFile::Open(fileOut.c_str(), recreate ? "recreate" : "update");
  tfOut->mkdir("allEventsCounter", tfIn->Get<TDirectory>("allEventsCounter")->GetTitle());
  tfOut->cd("allEventsCounter");
  tfIn->Get<TH1>("allEventsCounter/totalEvents")->Write();
  tfOut->cd("/");
  tfOut->mkdir("Original", "The distributions of unfiltered entries");
  tfOut->cd("Original");
  for (auto &&histView: vHistViewOriginal) {
    histView->Write();
  }
  if (debug) std::cerr << "Completed!" << std::endl;
  tfOut->Close();
  tfIn->Close();
}

#if false
int main(int argc, char** argv) {
  xAna_monoZ_preselect(argv[0], argv[1], true, true);
}
#endif
