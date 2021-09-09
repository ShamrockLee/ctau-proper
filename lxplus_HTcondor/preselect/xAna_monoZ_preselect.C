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
#include <TError.h>

#include <string>
#include <string_view>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <regex>

template<typename V = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumn(
    D &df,
    const std::string nameColumn,
    const std::string typenameColumn,
    const std::string expression,
    const Int_t binDensityOrder,
    const Double_t alignment,
    const Bool_t isLowerAssigned,
    const Int_t lowerLimitBins,
    const Bool_t isUpperAssigned,
    const Int_t upperLimitBins) {
  const Int_t nBins = TMath::Nint((upperLimitBins - lowerLimitBins) * TMath::Power(10., binDensityOrder));
  const Double_t binWidth = TMath::Power(10., -binDensityOrder);
  const Double_t lowerLimit = alignment + binWidth * lowerLimitBins;
  const Double_t upperLimit = alignment + binWidth * upperLimitBins;
  return df.template Histo1D<V>({("h" + nameColumn).c_str(),
      Form("{name: %s, typename: %s, expression: %s, binDensityOrder: %d, alignment: %F, isLowerAssigned: %s, lowerLimitBins: %d, isUpperAssigned: %s, upperLimitBins: %d}",
          nameColumn.c_str(), typenameColumn.c_str(), expression.c_str(), binDensityOrder, alignment, isLowerAssigned ? "true" : "false", lowerLimitBins, isUpperAssigned ? "true" : "false", upperLimitBins),
      nBins, lowerLimit, upperLimit}, expression);
}

template<typename V = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumn(D &df, const std::string nameColumn, const std::string nameColumnStripped) {
  const std::string typenameColumn = df.GetColumnType(nameColumn);
  const TString tstrTypenameColumn = typenameColumn;
  const TString tstrNameColumnStripped = nameColumnStripped;
  Int_t binDensityOrder = 0;
  Double_t alignment = 0.;
  Bool_t isLowerAssigned = false;
  Bool_t isUpperAssigned = false;
  Int_t halfWidth = 10000;
  Int_t lowerLimitBins = -halfWidth;
  Int_t upperLimitBins = halfWidth;
  Int_t nBins;
  Double_t lowerLimit, upperLimit;
  if (tstrTypenameColumn.Contains("Bool") || tstrTypenameColumn.Contains("bool")) {
    binDensityOrder = 0;
    alignment = 0.;
    isLowerAssigned = true;
    lowerLimitBins = 0.;
    isUpperAssigned = true;
    upperLimitBins = 2.;
  } else if (tstrTypenameColumn.Contains("Int") || tstrTypenameColumn.Contains("int")) {
    binDensityOrder = 0;
    alignment = -0.5;
    if (tstrNameColumnStripped.Contains("Idx")
      || tstrNameColumnStripped.Contains("Rank")) {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    } else if (tstrNameColumnStripped.BeginsWith("n")
      || tstrNameColumnStripped.Contains("nJet")
      || tstrNameColumnStripped.Contains("nEle")
      || tstrNameColumnStripped.Contains("nMu")
      || tstrNameColumnStripped.Contains("nPho")
      || tstrNameColumnStripped.Contains("_n")) {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    }
  } else {
    alignment = 0.;
    if (tstrNameColumnStripped.EndsWith("Pt")) {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    } else if (tstrNameColumnStripped.EndsWith("Eta")) {
      binDensityOrder = 1;
    } else if (tstrNameColumnStripped.EndsWith("Phi")) {
      binDensityOrder = 2;
      isUpperAssigned = true;
      upperLimitBins = TMath::CeilNint(TMath::Pi() * TMath::Power(10., binDensityOrder));
      isLowerAssigned = true;
      lowerLimitBins = -upperLimitBins;
    }  else if (tstrNameColumnStripped.EndsWith("E")) {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    } else if (tstrNameColumnStripped.EndsWith("M")
      || tstrNameColumnStripped.EndsWith("M2")
      || tstrNameColumnStripped.EndsWith("Mt")
      || tstrNameColumnStripped.EndsWith("Mt2")) {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    } else if (tstrNameColumnStripped.EndsWith("Sig")) {
      binDensityOrder = 2;
      isLowerAssigned = true;
      lowerLimitBins = 0;
    } else if (tstrNameColumnStripped.EndsWith("jetArea")) {
      binDensityOrder = 2;
      isLowerAssigned = true;
      lowerLimitBins = 0;
    } else if (tstrNameColumnStripped.EndsWith("EF")) {
      binDensityOrder = 2;
      isLowerAssigned = 0;
      lowerLimitBins = 0;
      isUpperAssigned = 0;
      upperLimitBins = TMath::Power(10, binDensityOrder);
    } else if (tstrNameColumnStripped(0, tstrNameColumnStripped.Length() - 2).String().EndsWith("Cov")) {
      binDensityOrder = 2;
      isLowerAssigned = 0;
      lowerLimitBins = 0;
      isUpperAssigned = 0;
      upperLimitBins = TMath::Power(10, binDensityOrder);
    } else if (tstrNameColumnStripped.Contains("Btag")) {
      binDensityOrder = 2;
    } else if (tstrNameColumnStripped.Contains("Ctag")) {
      binDensityOrder = 2;
    }
  }
  return GetHistFromColumn<V, D>(
      df, nameColumn, typenameColumn, nameColumn,
      binDensityOrder, alignment,
      isLowerAssigned, lowerLimitBins,
      isUpperAssigned, upperLimitBins);
}

template<typename V = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumn(D &df, const std::string nameColumn) {
  return GetHistFromColumn<V, D>(df, nameColumn, nameColumn);
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

template<typename E>
E SumRVecWithInit(const ROOT::RVec<E> v, const E init) {
  return std::accumulate(v.begin(), v.end(), init);
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

  const std::array<std::string, 2> aPrefLepFlav{"Ele", "Mu"};
  const std::array<std::string, 2> aPrefLepFlavLower{"ele", "mu"};
  const std::array<std::string, 2> aPrefAKShort{"THIN", "FAT"};

  TFile *tfIn = TFile::Open(fileIn.c_str());
  TTree *ttIn = tfIn->Get<TTree>("tree/treeMaker");
  ROOT::RDataFrame dfIn(*ttIn);
  std::vector<std::string> vNameColOriginal = {};
  std::array<std::vector<std::string>, 2> avNameColNumCorrect, avNameColZMassCutted;
  for (auto pav: {&avNameColNumCorrect, &avNameColZMassCutted}) {
    for (auto &v: *pav) {
      v.clear();
    }
  }
  std::array<std::array<std::vector<std::string>, 2>, 2> aavNameColPreselected, aavNameColMatching, aavNameColAllMatched;
  for (auto paav: {&aavNameColPreselected, &aavNameColMatching, &aavNameColAllMatched}) {
    for (auto &av: *paav) {
      for (auto &v: av) {
        v.clear();
      }
    }
  }

  std::vector<ROOT::RDF::RResultPtr<TH1D>> vHistViewOriginal = {};
  std::array<std::vector<ROOT::RDF::RResultPtr<TH1D>>, 2> avHistViewNumCorrect, avHistViewZMassCutted;
  for (auto pav: {&avHistViewNumCorrect, &avHistViewZMassCutted}) {
    for (auto &v: *pav) {
      v.clear();
    }
  }
  std::array<std::array<std::vector<ROOT::RDF::RResultPtr<TH1D>>, 2>, 2> aavHistViewPreselected, aavHistViewMatching, aavHistViewAllMatched;
  for (auto paav: {&aavHistViewPreselected, &aavHistViewMatching, &aavHistViewAllMatched}) {
    for (auto &av: *paav) {
      for (auto &v: av) {
        v.clear();
      }
    }
  }

  ROOT::RDF::RNode dfOriginal(dfIn);
  {
    auto vNamesCol = dfOriginal.GetColumnNames();
    for (const std::string nameCol: vNamesCol) {
      const std::string typenameCol(dfOriginal.GetColumnType(nameCol));
      const TString tstrTypenameCol(typenameCol.c_str());
      std::string typenameE;
      if (!RefgetE1D(typenameE, typenameCol)) {
        if (!tstrTypenameCol.Contains("TClonesArray")
          && !tstrTypenameCol.Contains("TObjArray")) {
          // if (debug) std::cerr << "Working on " << nameCol << ", " << typenameCol << std::endl;
          if (tstrTypenameCol.Contains("Double")
            || tstrTypenameCol.Contains("double")
            || tstrTypenameCol.Contains("Float")
            || tstrTypenameCol.Contains("float")
            || tstrTypenameCol.Contains("Long")
            || tstrTypenameCol.Contains("long")
            || tstrTypenameCol.Contains("Short")
            || tstrTypenameCol.Contains("short")
            || tstrTypenameCol.Contains("Int")
            || tstrTypenameCol.Contains("int")
            || tstrTypenameCol.Contains("Char")
            || tstrTypenameCol.Contains("Byte")
            || tstrTypenameCol.Contains("char")
            || tstrTypenameCol.Contains("Bool")
            || tstrTypenameCol.Contains("bool")
          ) {
            vNameColOriginal.emplace_back(nameCol);
            for (auto paav: {&aavNameColPreselected, &aavNameColMatching, &aavNameColAllMatched}) {
              for (auto &av: *paav) {
                for (auto &v: av) {
                  v.emplace_back(nameCol);
                }
              }
            }
          }
        }
      }
    }
    for (const std::string nameCol: vNamesCol) {
      const TString tstrNameCol(nameCol);
      // if (debug) std::cerr << tstrNameCol << "\tLength: " << tstrNameCol.Length() << std::endl;
      if (tstrNameCol.EndsWith("P4")) {
        const TString tstrPrefNameCol = tstrNameCol(0, tstrNameCol.Length() - 2);
        const std::string prefNameCol (tstrPrefNameCol.Data());
        if (debug) std::cerr << "prefNameCol: " << prefNameCol << std::endl;
        dfOriginal = dfOriginal
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
  dfOriginal = dfOriginal
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
    auto vNamesCol = dfOriginal.GetColumnNames();
    std::vector<std::string> vTypenamesCol(vNamesCol.size());
    vTypenamesCol.clear();
    for (const auto nameCol: vNamesCol){
      vTypenamesCol.emplace_back(dfOriginal.GetColumnType(nameCol));
      // if (debug) std::cerr << vTypenamesCol.back() << ",\tsize: " << vTypenamesCol.size() << std::endl;
    }
    for (std::string &nameCol: vNamesCol) {
      std::regex r("^is(.+)Muon$");
      std::smatch m;
      if (std::regex_match(nameCol, m, r)) {
        const std::string nameAlias = "muIsPass" + std::string(m[1]);
        if (debug) std::cerr << "nameAlias: " << nameAlias << std::endl;
        dfOriginal = dfOriginal.Define(nameAlias, nameCol);
        nameCol = nameAlias;
      }
    }
    for (size_t iCol = 0; iCol < vNamesCol.size(); ++iCol) {
      const std::string &nameCol = vNamesCol[iCol];
      const std::string &typenameCol = vTypenamesCol[iCol];
      // if (debug) std::cerr << "Working on " << nameCol << ",\t" << typenameCol << std::endl;
      const TString tstrNameCol = nameCol;
      const TString tstrTypenameCol = typenameCol;
      std::string typenameE = "";
      const Int_t nDims = RefgetE2D(typenameE, typenameCol) ? 2 : (RefgetE1D(typenameE, typenameCol) ? 1 : 0);
      // if (debug) std::cerr << "nDims: " << nDims << std::endl;
      const TString tstrTypenameE = typenameE;
      for (const std::string pref: {"ele", "mu", "THINjet", "FATjet"}) {
        if (tstrNameCol.BeginsWith(pref)) {
          if (nDims == 2) {
            if (debug) std::cerr << "Redefining 2D: " << nameCol << std::endl;
            dfOriginal = dfOriginal
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
          } else if (nDims == 1) {
            if (debug) std::cerr << "Redefining 1D: " << nameCol << std::endl;
            dfOriginal = dfOriginal
            .Redefine(nameCol, Form("ROOT::VecOps::Take(%s, %sIdxPassPtEta)", nameCol.c_str(), pref.c_str()));
            if (tstrTypenameE.Contains("double")
              || tstrTypenameE.Contains("float")
              || tstrTypenameE.Contains("long")
              || tstrTypenameE.Contains("int")
              || tstrTypenameE.Contains("short")
              || tstrTypenameE.Contains("char")
              // || tstrTypenameE.Contains("bool")
            ) {
              if (debug) std::cerr << "Pick " << nameCol << std::endl;
              vNameColOriginal.emplace_back(nameCol);
              for (auto paav: {&aavNameColPreselected, &aavNameColMatching, &aavNameColAllMatched}) {
                for (auto &av: *paav) {
                  for (auto &v: av) {
                    v.emplace_back(nameCol);
                  }
                }
              }
            }
          }
          break;
        } else {
          continue;
        }
      }
    }
    for (std::size_t iCol = 0; iCol < vNamesCol.size(); ++iCol) {
      const std::string &nameCol = vNamesCol[iCol];
      const std::string &typenameCol = vTypenamesCol[iCol];
      const TString tstrNameCol = nameCol;
      const TString tstrTypenameCol = typenameCol;
      for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
        if (tstrNameCol.BeginsWith(aPrefLepFlavLower[iLepFlav])) {
          if (tstrNameCol.BeginsWith((aPrefLepFlavLower[iLepFlav] + "Is"))) {
            std::string nameColNew = "n" + aPrefLepFlav[iLepFlav] + tstrNameCol(aPrefLepFlavLower[iLepFlav].length() + 2, nameCol.length()).Data();
            if (debug) std::cerr << "Defining " << nameColNew << std::endl;
            dfOriginal = dfOriginal.Define(nameColNew,
              [](const ROOT::RVec<Bool_t> rvId)->Int_t{
                return ROOT::VecOps::Sum(rvId, 0);
              }, {nameCol});
            vNameColOriginal.emplace_back(nameColNew);
            for (auto paav: {&aavNameColPreselected, &aavNameColMatching, &aavNameColAllMatched}) {
              for (auto &av: *paav) {
                for (auto &v: av) {
                  v.emplace_back(nameCol);
                }
              }
            }
          }
          break;
        } else {
          continue;
        }
      }
    }
  }
  // vNameColOriginal.emplace_back("eleIsPassLoose");
  dfOriginal = dfOriginal
  .Define("leptonPairFlavor", [&debug](const Int_t nElePassLoose, const Int_t nElePassMedium, const Int_t nElePassTight, const Int_t nMuPassLoose, const Int_t nMuPassMedium, const Int_t nMuPassTight)->Int_t{
    // if (debug) Info("lambda:leptonPairFlavor", "%d, %d, %d, %d, %d, %d", nElePassLoose, nElePassMedium, nElePassTight, nMuPassLoose, nMuPassMedium, nMuPassTight);
    if (nElePassTight >= 4 || nMuPassTight >= 4) {
      return -1;
    }
    const Bool_t isEleNumCorrect = nElePassMedium >= 1 && nElePassLoose >= 2;
    const Bool_t isMuNumCorrect = nMuPassMedium >= 1 && nMuPassLoose >= 2;
    if (isEleNumCorrect && isMuNumCorrect) {
      return -1;
    } else if (isEleNumCorrect) {
      return 0;
    } else if (isMuNumCorrect) {
      return 1;
    }
    return -1;
  },{"nElePassLoose", "nElePassMedium", "nElePassTight", "nMuPassLoose", "nMuPassMedium", "nMuPassTight"});  std::array<ROOT::RDF::RNode, 2> aDfNumCorrect = {dfOriginal, dfOriginal};
  vNameColOriginal.emplace_back("leptonPairFlavor");
  {
    for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
      aDfNumCorrect[iLepFlav] = aDfNumCorrect[iLepFlav]
      .Filter(Form("leptonPairFlavor==%zu", iLepFlav))
      .Define(aPrefLepFlavLower[iLepFlav] + "PairedIdx", [&debug](const Int_t nLep, const ROOT::RVec<Bool_t> isLoose, const ROOT::RVec<Bool_t> isMedium, const ROOT::RVec<Bool_t> isTight)->ROOT::RVec<Int_t>{
        // if (debug) Info("lambda:PairedIdx", "nLep: %d, (%zu, %zu, %zu)", nLep, isLoose.size(), isMedium.size(), isTight.size());
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
      }, {"n" + aPrefLepFlav[iLepFlav], aPrefLepFlavLower[iLepFlav] + "IsPassLoose", aPrefLepFlavLower[iLepFlav] + "IsPassMedium", aPrefLepFlavLower[iLepFlav] + "IsPassTight"});
      // if (debug) {
      //   std::cerr << "Columns in dfNumCorrect[" << iLepFlav << "]: " << std::flush;
      //   for (const auto nameCol: aDfNumCorrect[iLepFlav].GetColumnNames()) {
      //     std::cerr << nameCol << ", " << std::flush;
      //   }
      //   std::cerr << "}" << std::endl;
      // }
      for (const auto nameCol : aDfNumCorrect[iLepFlav].GetColumnNames()) {
        const TString tstrNameCol(nameCol);
        const TString tstrTypenameCol(aDfNumCorrect[iLepFlav].GetColumnType(nameCol).c_str());
        if (tstrNameCol.BeginsWith(aPrefLepFlavLower[iLepFlav].c_str()) && tstrTypenameCol.BeginsWith("ROOT::VecOps::RVec") && ! tstrNameCol.BeginsWith(aPrefLepFlavLower[iLepFlav] + "PairedIdx")) {
          const std::string nameColNew = aPrefLepFlavLower[iLepFlav] + "Paired" + tstrNameCol(aPrefLepFlavLower[iLepFlav].length(), tstrNameCol.Length()).Data();
          if (aDfNumCorrect[iLepFlav].HasColumn(nameColNew)) continue;
          if (debug) std::cerr << "Defining " << nameColNew << std::endl;
          avNameColNumCorrect[iLepFlav].emplace_back(nameColNew);
          for (auto paav: {&aavNameColPreselected, &aavNameColMatching, &aavNameColAllMatched}) {
            for (auto &v: (*paav)[iLepFlav]) {
              v.emplace_back(nameColNew);
            }
          }
          aDfNumCorrect[iLepFlav] = aDfNumCorrect[iLepFlav]
          .Define(nameColNew, "ROOT::VecOps::Take(" + nameCol + ", " + aPrefLepFlavLower[iLepFlav] + "PairedIdx)");
        }
      }
      aDfNumCorrect[iLepFlav] = aDfNumCorrect[iLepFlav]
      .Define(
        aPrefLepFlavLower[iLepFlav] + "PairM",
        "static_cast<Float_t>(ROOT::VecOps::Sum(ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiEVector>("
        "  " + aPrefLepFlavLower[iLepFlav] + "PairedPt,"
        "  " + aPrefLepFlavLower[iLepFlav] + "PairedEta,"
        "  " + aPrefLepFlavLower[iLepFlav] + "PairedPhi,"
        "  " + aPrefLepFlavLower[iLepFlav] + "PairedE"
        "), ROOT::Math::PtEtaPhiEVector()).M())")
      .Define(aPrefLepFlavLower[iLepFlav] + "PairIsPassZ",
        "TMath::Abs(" + aPrefLepFlavLower[iLepFlav] + "PairM" + " - " + std::to_string(massZ) + ") < 20.");
      avNameColNumCorrect[iLepFlav].emplace_back(aPrefLepFlavLower[iLepFlav] + "PairM");
      avNameColNumCorrect[iLepFlav].emplace_back(aPrefLepFlavLower[iLepFlav] + "PairIsPassZ");
      for (auto paav: {&aavNameColPreselected, &aavNameColMatching, &aavNameColAllMatched}) {
        for (auto &v: (*paav)[iLepFlav]) {
          v.emplace_back(aPrefLepFlavLower[iLepFlav] + "PairM");
        }
      }
    }
  }  std::array<ROOT::RDF::RNode, 2> aDfZMassCutted = aDfNumCorrect;
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    aDfZMassCutted[iLepFlav] = aDfZMassCutted[iLepFlav].Filter(aPrefLepFlavLower[iLepFlav] + "PairIsPassZ");
  }
  // https://stackoverflow.com/questions/28541488/nested-aggregate-initialization-of-stdarray
  std::array<std::array<ROOT::RDF::RNode, 2>, 2> aaDfPreselected {{ {aDfZMassCutted[0], aDfZMassCutted[0]}, {aDfZMassCutted[1], aDfZMassCutted[1]} }};
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      aaDfPreselected[iLepFlav][iAK] = aaDfPreselected[iLepFlav][iAK].Filter(aPrefAKShort[iAK] + "nJet" + " > 2");
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    const size_t nCol = avNameColNumCorrect[iLepFlav].size();
    avHistViewNumCorrect[iLepFlav].clear();
    avHistViewNumCorrect[iLepFlav].reserve(nCol);
    for (const std::string &nameCol: avNameColNumCorrect[iLepFlav]) {
      avHistViewNumCorrect[iLepFlav].emplace_back(GetHistFromColumn(aDfNumCorrect[iLepFlav],nameCol));
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    const size_t nCol = avNameColZMassCutted[iLepFlav].size();
    avHistViewZMassCutted[iLepFlav].clear();
    avHistViewZMassCutted[iLepFlav].reserve(nCol);
    for (const std::string &nameCol: avNameColZMassCutted[iLepFlav]) {
      avHistViewZMassCutted[iLepFlav].emplace_back(GetHistFromColumn(aDfZMassCutted[iLepFlav], nameCol));
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      const size_t nCol = aavHistViewPreselected[iLepFlav][iAK].size();
      aavHistViewPreselected[iLepFlav][iAK].clear();
      aavHistViewPreselected[iLepFlav][iAK].reserve(nCol);
      for (const std::string &nameCol: aavNameColPreselected[iLepFlav][iAK]) {
        aavHistViewPreselected[iLepFlav][iAK].emplace_back(GetHistFromColumn(aaDfPreselected[iLepFlav][iAK], nameCol));
      }
    }
  }
  TFile* tfOut = TFile::Open(fileOut.c_str(), recreate ? "recreate" : "update");
  tfOut->mkdir("allEventsCounter", tfIn->Get<TDirectory>("allEventsCounter")->GetTitle())->cd();
  tfIn->Get<TH1>("allEventsCounter/totalEvents")->Write();
  {
    const size_t nCol = vNameColOriginal.size();
    vHistViewOriginal.clear();
    vHistViewOriginal.reserve(nCol);
    for (const std::string &nameCol: vNameColOriginal) {
      vHistViewOriginal.emplace_back(GetHistFromColumn(dfOriginal, nameCol));
    }
  }
  tfOut->cd("/");
  tfOut->mkdir("Original", "Unfiltered entries")->cd();
  for (auto &&histView: vHistViewOriginal) {
    histView->Write();
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->cd("/");
    tfOut->mkdir(("NumCorrect" + aPrefLepFlav[iLepFlav]).c_str(), ("Entries with correct " + aPrefLepFlavLower[iLepFlav] + " numbers").c_str())->cd();
    for (auto &&histView: avHistViewNumCorrect[iLepFlav]) {
      histView->Write();
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->cd("/");
    tfOut->mkdir(("ZMassCutted" + aPrefLepFlav[iLepFlav]).c_str(), "Entries passing the Z mass cut")->cd();
    for (auto &&histView: avHistViewZMassCutted[iLepFlav]) {
      histView->Write();
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      tfOut->cd("/");
      tfOut->mkdir(("Preselected" + aPrefLepFlav[iLepFlav] + aPrefAKShort[iAK] + "jet").c_str(), "Preselected entries")->cd();
      for (auto &&histView: aavHistViewPreselected[iLepFlav][iAK]) {
        histView->Write();
      }
    }
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
