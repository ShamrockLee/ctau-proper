#include <Rtypes.h>
#include <TROOT.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDF/HistoModels.hxx>
#include <Math/Vector4D.h>
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

/**Calculate deltaR2 of two ROOT::Math::LorentzVector
 * \\(\\Delta R ^2= \\Delta \\eta ^2 + \\Delta \\phi ^2\\)
 */
template<typename TypeVector1, typename TypeVector2>
inline Double_t GetDeltaR2(TypeVector1 lv1, TypeVector2 lv2) {
  return 
    (static_cast<Double_t>(lv1.Eta()) - static_cast<Double_t>(lv2.Eta())) *
    (static_cast<Double_t>(lv1.Eta()) - static_cast<Double_t>(lv2.Eta())) +
    (static_cast<Double_t>(lv1.Phi()) - static_cast<Double_t>(lv2.Phi())) *
    (static_cast<Double_t>(lv1.Phi()) - static_cast<Double_t>(lv2.Phi()));
}

/**Generate a histogram of the specified expression lazily
 * from an RDataFram
 * by applying Histo1D() with specified binning information
 * 
 * This overload does the actual application
 * and store the binning specification
 * in the histogram title in JSON format
 */
template<typename V = ROOT::RDFDetail::RInferredType, typename W = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumnCustom(
    D &df,
    const std::string nameColumn,
    const std::string typenameColumn,
    const Int_t binDensityOrder,
    const Double_t alignment,
    const Bool_t isLowerAssigned,
    const Int_t lowerLimitBins,
    const Bool_t isUpperAssigned,
    const Int_t upperLimitBins,
    const std::string expression = "",
    const std::string exprWeight = "") {
  const Int_t nBins = TMath::Nint((upperLimitBins - lowerLimitBins) * TMath::Power(10., binDensityOrder));
  const Double_t binWidth = TMath::Power(10., -binDensityOrder);
  const Double_t lowerLimit = alignment + binWidth * lowerLimitBins;
  const Double_t upperLimit = alignment + binWidth * upperLimitBins;
  const char*  &&cstrJSON = Form("{name:%s,typename:%s,expression:%s,%sbinDensityOrder:%d,alignment:%F,isLowerAssigned:%s,lowerLimitBins:%d,isUpperAssigned:%s,upperLimitBins:%d}",
          nameColumn.c_str(), typenameColumn.c_str(), expression.c_str(), exprWeight.length() ? ("exprWeight:" + exprWeight + ",").c_str() : "",
          binDensityOrder, alignment, isLowerAssigned ? "true" : "false", lowerLimitBins, isUpperAssigned ? "true" : "false", upperLimitBins);
  const ROOT::RDF::TH1DModel model {("h" + nameColumn).c_str(), cstrJSON,
      nBins, lowerLimit, upperLimit};
  return exprWeight.length()
    ? df.template Histo1D<V, W>(model, expression, exprWeight)
    : df.template Histo1D<V>(model, expression);
}

/**Generate a histogram of the specified expression lazily
 * from an RDataFram by applying Histo1D()
 * 
 * The binning information (an ROOT::RDF::TH1DModel insteance)
 * is calculated from the column name or another string specified.
 * 
 * This overload decide the bin number and lower and upper limits
 * from the specified string (`nameColumnStripped`)
 */
template<typename V = ROOT::RDFDetail::RInferredType, typename W = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumnCustom(D &df, const std::string nameColumn, const std::string typenameColumn, const std::string expression, std::string exprWeight = "", const std::string nameColumnStripped = "") {
  if (!nameColumnStripped.length()) {
    return GetHistFromColumnCustom<V, W, D>(df, nameColumn, typenameColumn, expression, exprWeight, nameColumn);
  }
  const TString tstrTypenameColumn = typenameColumn;
  const TString tstrNameColumnStripped = nameColumnStripped;
  Int_t binDensityOrder = 0;
  Double_t alignment = 0.;
  Bool_t isLowerAssigned = false;
  Bool_t isUpperAssigned = false;
  Int_t lowerLimitBins = -1000;
  Int_t upperLimitBins = 5000;
  Int_t nBins;
  Double_t lowerLimit, upperLimit;
  if (nameColumn == "Counter") {
    binDensityOrder = 0;
    alignment = 0;
    isLowerAssigned = true;
    lowerLimitBins = 0;
    isUpperAssigned = true;
    upperLimitBins = 1;
  } else if (tstrTypenameColumn.Contains("Bool") || tstrTypenameColumn.Contains("bool")) {
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
    } else if (tstrNameColumnStripped.EndsWith("Sgn")) {
      isLowerAssigned = true;
      lowerLimitBins = -1;
      isUpperAssigned = true;
      upperLimitBins = 2;
    }
  } else {
    alignment = 0.;
    if (nameColumn == "mcWeight") {
      binDensityOrder = 2;
    } else if (tstrNameColumnStripped.EndsWith("Pt")) {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    } else if (tstrNameColumnStripped.EndsWith("Rho")) {
      binDensityOrder = 1;
      isLowerAssigned = true;
      lowerLimitBins = 0;
    } else if (tstrNameColumnStripped.EndsWith("Eta")) {
      binDensityOrder = 1;
      lowerLimitBins = -upperLimitBins;
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
      isLowerAssigned = true;
      lowerLimitBins = 0;
      isUpperAssigned = true;
      upperLimitBins = TMath::Power(10, binDensityOrder);
    } else if (tstrNameColumnStripped(0, tstrNameColumnStripped.Length() - 2).String().EndsWith("Cov")) {
      binDensityOrder = 2;
      isLowerAssigned = true;
      lowerLimitBins = 0;
      isUpperAssigned = true;
      upperLimitBins = TMath::Power(10, binDensityOrder);
    } else if (tstrNameColumnStripped.Contains("Btag")) {
      binDensityOrder = 2;
    } else if (tstrNameColumnStripped.Contains("Ctag")) {
      binDensityOrder = 2;
    }
  }
  return GetHistFromColumnCustom<V, W, D>(
      df, nameColumn, typenameColumn,
      binDensityOrder, alignment,
      isLowerAssigned, lowerLimitBins,
      isUpperAssigned, upperLimitBins,
      expression, exprWeight);
}

/**Generate a histogram of an RDataFrame column lazily
 * from an RDataFram by applying Histo1D()
 * 
 * The binning information (an ROOT::RDF::TH1DModel insteance)
 * is calculated from the column name or another string specified.
 */
template<typename V = ROOT::RDFDetail::RInferredType, typename W = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumn(D &df, const std::string nameColumn, std::string exprWeight = "", const std::string nameColumnStripped = "") {
  const std::string typenameColumn = df.GetColumnType(nameColumn);
  return GetHistFromColumnCustom<V, W, D>(df, nameColumn, typenameColumn, nameColumn, exprWeight, nameColumnStripped);
}

/**Determine if a type name string describes a 2D Rvec
 * and get the element type name string
 */
Bool_t RefgetE2D(std::string &typenameE, const std::string typenameCol) {
  std::regex r ("^ROOT::VecOps::RVec\\s*<\\s*(vector|ROOT::VecOps::RVec)\\s*<\\s*(.*\\S)\\s*>\\s*>\\s*$");
  std::smatch m;
  const Bool_t result = (std::regex_match(typenameCol, m, r));
  if (result) {
    typenameE = m[2];
  }
  return result;
}

/**Determine if a type name string describes a 1D Rvec
 * and get the element type name string
 */
Bool_t RefgetE1D(std::string &typenameE, const std::string typenameCol) {
  std::regex r ("^ROOT::VecOps::RVec\\s*<\\s*(.*\\S)\\s*>\\s*$");
  std::smatch m;
  const Bool_t result = (std::regex_match(typenameCol, m, r));
  if (result) {
    typenameE = m[1];
  }
  return result;
}

inline std::string GetExprTakeIdx2D(const std::string nameCol, const std::string nameIdx, const std::string typenameE, const Bool_t debug = false) {
  return
    // "if (true && " + nameCol + ".size() <= ROOT::VecOps::Max(" + nameIdx + ")) {"
    // "  Fatal(\"lambda:RedefineWithIdx\", \"" + nameIdx + " maximum (%d) exceeds the " + nameCol + " size (%zu)\", ROOT::VecOps::Max(" + nameIdx + "), " + nameCol + ".size());"
    // "}"
    "ROOT::VecOps::RVec<ROOT::VecOps::RVec<" + typenameE + ">> vvResult(" + nameCol + ".size());"
    "vvResult.clear();"
    "for (size_t iVE: " + nameIdx + ") {"
    "  const auto& vE = " + nameCol + "[iVE];"
    "  ROOT::VecOps::RVec<" + typenameE + "> resultSub(vE.size());"
    "  resultSub.clear();"
    "  for (auto iterElement = vE.cbegin(); iterElement != vE.cend(); ++iterElement) {"
    "    resultSub.push_back(*iterElement);"
    "  } "
    "  vvResult.push_back(resultSub);"
    "} "
    "return vvResult"
    // // Doesn't work
    // "ROOT::VecOps::Map("
    // + nameCol +
    // ", [" + nameIdx + "](ROOT::RVec<" + typenameE + "> vSub){"
    // "  return ROOT::VecOps::Take(vSub, " + nameIdx + ");"
    // "})"
  ;
}

inline std::string GetExprTakeIdx1D(const std::string nameCol, const std::string nameIdx, const std::string typenameE, const Bool_t debug = false) {
  return Form("ROOT::VecOps::Take(%s, %s)", nameCol.c_str(), nameIdx.c_str());
}

template <class D = ROOT::RDF::RNode>
void RedefinePrefWithIdx(D &df, const std::string pref, const std::vector<std::string> vPrefExcl, const std::string nameIdx, const std::function<void(const std::string nameCol, const Int_t nDims, const std::string typenameE)> funPickNameCol = nullptr, const Bool_t debug = false) {
  for (const std::string &nameCol: df.GetColumnNames()) {
    if (std::find(vPrefExcl.cbegin(), vPrefExcl.cend(), pref) != vPrefExcl.cend()) continue;
    std::string typenameCol = df.GetColumnType(nameCol);
    const TString tstrNameCol = nameCol;
    std::string typenameE = "";
    const Int_t nDims = RefgetE2D(typenameE, typenameCol) ? 2 : (RefgetE1D(typenameE, typenameCol) ? 1 : 0);
    const TString tstrTypenameE = typenameE;
    if (tstrNameCol.BeginsWith(pref)) {
      if (nDims == 2) {
        if (debug) std::cerr << "Redefining 2D: " << nameCol << " (" << typenameE << ") (" << typenameCol << ")" << std::endl;
        df = df
        .Redefine(nameCol, GetExprTakeIdx2D(nameCol, nameIdx, typenameE, debug));
      } else if (nDims == 1) {
        if (debug) std::cerr << "Redefining 1D: " << nameCol << " (" << typenameE << ") (" << typenameCol << ")" << std::endl;
        df = df
        .Redefine(nameCol, GetExprTakeIdx1D(nameCol, nameIdx, typenameE, debug));
      }
      if (funPickNameCol && nDims > 0) {
        funPickNameCol(nameCol, nDims, typenameE);
      }
    }
  }
}

void xAna_monoZ_preselect(const std::string fileIn, const std::string fileOut, const Bool_t enableIMT=true, const Bool_t debug=false) {
  typedef ROOT::Math::PtEtaPhiMVector TypeLorentzVector;
  const std::string typenameLorentzVector = "ROOT::Math::PtEtaPhiMVector";
  if (enableIMT) ROOT::EnableImplicitMT();
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
  std::array<std::vector<std::string>, 2> avNameColGen, avNameColHasLPair, avNameColHasVtx, avNameColNoTau, avNameColLPairedPassPt, avNameColZMassCutted, avNameColNoExtraL;
  for (auto pav: {&avNameColGen, &avNameColHasLPair, &avNameColHasVtx, &avNameColNoTau, &avNameColLPairedPassPt, &avNameColZMassCutted, &avNameColNoExtraL}) {
    for (auto &v: *pav) {
      v.clear();
    }
  }
  std::array<std::array<std::vector<std::string>, 2>, 2> aavNameColHasJet, aavNameColLPairPassPt, aavNameColMatching, aavNameColAllMatched;
  for (auto paav: {&aavNameColHasJet, &aavNameColMatching, &aavNameColAllMatched}) {
    for (auto &av: *paav) {
      for (auto &v: av) {
        v.clear();
      }
    }
  }

  std::vector<ROOT::RDF::RResultPtr<TH1D>> vHistViewOriginal = {};
  std::array<std::vector<ROOT::RDF::RResultPtr<TH1D>>, 2> avHistViewGen, avHistViewHasLPair, avHistViewHasVtx, avHistViewNoTau, avHistViewLPairedPassPt, avHistViewZMassCutted, avHistViewNoExtraL;
  for (auto pav: {&avHistViewGen, &avHistViewHasLPair, &avHistViewHasVtx, &avHistViewNoTau, &avHistViewLPairedPassPt, &avHistViewZMassCutted, &avHistViewNoExtraL}) {
    for (auto &v: *pav) {
      v.clear();
    }
  }
  std::array<std::array<std::vector<ROOT::RDF::RResultPtr<TH1D>>, 2>, 2> aavHistViewHasJet, aavHistViewLPairPassPt, aavHistViewMatching, aavHistViewAllMatched;
  for (auto paav: {&aavHistViewHasJet, &aavHistViewLPairPassPt, &aavHistViewMatching, &aavHistViewAllMatched}) {
    for (auto &av: *paav) {
      for (auto &v: av) {
        v.clear();
      }
    }
  }
  // Begin the Original stage
  ROOT::RDF::RNode dfOriginal(dfIn);
  dfOriginal = dfOriginal
  .Define("Counter", [](){return 0.5;}, {})
  // // Doesn't work
  // // Error in <TClonesArray::operator=>: cannot copy TClonesArray's when classes are different
  // // Causes segfault
  // .Define("HPSTauP4", "HPSTau_4Momentum")
  // A simple implementation of Sgn function
  // https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c/1903975#1903975
  .Define("mcWeightSgn", "static_cast<Int_t>((mcWeight > 0) - (0 > mcWeight))");
  {
    auto vNamesCol = dfOriginal.GetColumnNames();
    for (const std::string nameCol: vNamesCol) {
      const TString tstrNameCol(nameCol);
      const std::string typenameCol(dfOriginal.GetColumnType(nameCol));
      const TString tstrTypenameCol(typenameCol);
      std::string typenameE;
      if (!RefgetE1D(typenameE, typenameCol) && !tstrNameCol.BeginsWith("mcWeight")) {
        if (!tstrTypenameCol.Contains("TClonesArray")
          && !tstrTypenameCol.Contains("TObjArray")) {
          if (debug) std::cerr << "Working on " << nameCol << ", " << typenameCol << std::endl;
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
            for (auto pav: {&avNameColHasLPair, &avNameColHasVtx, &avNameColNoTau, &avNameColLPairedPassPt, &avNameColZMassCutted, &avNameColNoExtraL}) {
              for (auto &v: *pav) {
                v.emplace_back(nameCol);
              }
            }
            for (auto paav: {&aavNameColHasJet, &aavNameColLPairPassPt, &aavNameColMatching, &aavNameColAllMatched}) {
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
    if (debug) std::cerr << "avNameColHasVtx[0].size(): " << avNameColHasVtx[0].size() << std::endl;
    std::vector<std::string> vPref;
    vPref.clear();
    for (std::string nameCol: vNamesCol) {
      TString tstrNameCol(nameCol);
      const std::string typenameCol = dfOriginal.GetColumnType(nameCol);
      // if (debug) std::cerr << tstrNameCol << "\tLength: " << tstrNameCol.Length() << std::endl;
      if (tstrNameCol.EndsWith("P4") || nameCol == "HPSTau_4Momentum") {
        if (typenameCol == "TClonesArray") {
          dfOriginal = dfOriginal
          .Redefine(nameCol, [/*prefNameCol*/](const TClonesArray &tcaP4){
            // Info("lambda:P4TcaTlzToRVecLz", "Working on %s, prefNameCol: %s", tcaP4.GetName(), prefNameCol.c_str());
            const size_t n = tcaP4.GetSize();
            ROOT::RVec<TypeLorentzVector> vResult(n, TypeLorentzVector());
            for (size_t i = 0; i < n; ++i) {
              const auto &tlvP4 = *static_cast<TLorentzVector *>(tcaP4[i]);
              vResult[i].SetPxPyPzE(tlvP4.Px(), tlvP4.Py(), tlvP4.Pz(), tlvP4.E());
            }
            return vResult;
          }, {{ nameCol }});
        }
        // // It does not work, but I don't know why
        // for (const std::string suff: {"Pt", "Eta", "Phi", "E", "M"}) {
        //   const std::string expression = "ROOT::VecOps::Map("
        //       + nameCol + ", "
        //       "[](const ROOT::Math::PtEtaPhiMVector p4)->Double_t{"
        //       "return static_cast<Double_t>(p4." + suff + "());"
        //       "})";
        //   if (debug) std::cerr << expression << std::endl;
        //   dfOriginal = dfOriginal
        //   .Define(
        //       prefNameCol + suff,
        //       expression);
        // }
        TString tstrPrefNameCol = tstrNameCol(0, tstrNameCol.Length() - 2);
        if (debug) std::cerr << "prefNameCol: " << tstrPrefNameCol << std::endl;
        if (nameCol == "HPSTau_4Momentum") {
          nameCol = "HPSTauP4";
          tstrNameCol = nameCol;
          tstrPrefNameCol = "HPSTau";
          dfOriginal = dfOriginal.Define(nameCol, "HPSTau_4Momentum");
        }
        vPref.emplace_back(tstrPrefNameCol.Data());
      }
    }
    // std::vector<size_t> vIdxPrefSorted(vPref.size(), 0);
    // std::iota(vIdxPrefSorted.begin(), vIdxPrefSorted.end());
    // std::sort(vIdxPrefSorted.begin(), vIdxPrefSorted.end(), [&vPref](size_t i, size_t j){ return vPref[i] > vPref[j]; });
    for (std::string prefNameCol: { "ele", "mu", "THINjet", "FATjet", "HPSTau" }) {
      std::string nameCol = prefNameCol + "P4";
      TString tstrNameCol(nameCol);
      TString tstrPrefNameCol(prefNameCol);
      dfOriginal = dfOriginal
      .Define(prefNameCol + "Rank", [](const ROOT::RVec<TypeLorentzVector> &vP4) {
        ROOT::RVec<Int_t> vResult(vP4.size());
        std::iota(vResult.begin(), vResult.end(), 0);
        std::stable_sort(vResult.begin(), vResult.end(), [ &vP4 ](const Int_t ia, const Int_t ib)->Bool_t{ return vP4[ia].Pt() > vP4[ib].Pt(); });
        return vResult;
      }, {{ nameCol }})
      .Define("is" + prefNameCol + "Presorted", [](const ROOT::RVec<Int_t> vRank)->Bool_t{
        return std::is_sorted(vRank.begin(), vRank.end());
      }, {{ prefNameCol + "Rank" }});
      vNameColOriginal.emplace_back("is" + prefNameCol + "Presorted");
      if (!tstrPrefNameCol.Contains("gen") && !tstrPrefNameCol.Contains("Gen") && !tstrPrefNameCol.Contains("GEN")) {
        RedefinePrefWithIdx(dfOriginal, prefNameCol, { prefNameCol + "jet", prefNameCol + "Jet" + prefNameCol + "GenJet" }, prefNameCol + "Rank", nullptr, debug);
      }
      dfOriginal = dfOriginal
      .Define(prefNameCol + "Pt", [](const ROOT::RVec<TypeLorentzVector> &vP4) {
        return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4) {
          return p4.Pt();
        });
      }, {{ nameCol }})
      .Define(prefNameCol + "Eta", [](const ROOT::RVec<TypeLorentzVector> &vP4) {
        return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4) {
          return p4.Eta();
        });
      }, {{ nameCol }})
      .Define(prefNameCol + "Phi", [](const ROOT::RVec<TypeLorentzVector> &vP4) {
        return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4) {
          return p4.Phi();
        });
      }, {{ nameCol }})
      .Define(prefNameCol + "E", [](const ROOT::RVec<TypeLorentzVector> &vP4) {
        return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4) {
          return p4.E();
        });
      }, {{ nameCol }})
      .Define(prefNameCol + "M", [](const ROOT::RVec<TypeLorentzVector> &vP4) {
        return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4) {
          return p4.M();
        });
      }, {{ nameCol }});
      continue;
    }
  }
  dfOriginal = dfOriginal
  .Define("eleIdxPassBasicCuts", [](const Int_t nEle, const ROOT::RVec<TypeLorentzVector> &eleP4){
    ROOT::RVec<Int_t> vResult(nEle);
    vResult.clear();
    for (Int_t iEle = 0; iEle < nEle; ++iEle) {
      if (eleP4[iEle].Pt() > 20. && TMath::Abs(eleP4[iEle].Eta()) < 2.4) {
        vResult.emplace_back(iEle);
      }
    }
    return vResult;
  }, {"nEle", "eleP4"})
  .Define("nElePassBasicCuts", "static_cast<Int_t>(eleIdxPassBasicCuts.size())")
  .Redefine("nEle", "nElePassBasicCuts")
  // // JIT compilation doesn't pass without an obvious reason
  // .Define("muRelIso", "(muChHadIso + ROOT::VecOps::Map(muNeHadIso + muGamIso + 0.5 * muPUPt, [](const Double_t x){ return x > 0 ? x : 0.; })) / muPt")
  .Define("muRelIso", [](
      const ROOT::RVec<Float_t> &muChHadIso, const ROOT::RVec<Float_t> &muNeHadIso, const ROOT::RVec<Float_t> &muGamIso, const ROOT::RVec<Float_t> &muPUPt, const ROOT::RVec<Double_t> &muPt) {
    return (muChHadIso + ROOT::VecOps::Map(muNeHadIso + muGamIso + 0.5 * muPUPt, [](const Double_t x){ return x > 0 ? x : 0.; })) / muPt;
  }, { "muChHadIso", "muNeHadIso", "muGamIso", "muPUPt", "muPt" })
  .Define("muIdxPassBasicCuts", [](const Int_t nMu, const ROOT::RVec<TypeLorentzVector> &muP4, const ROOT::RVec<Double_t> &muRelIso, const ROOT::RVec<Int_t> muTrkLayers) {
    ROOT::RVec<Int_t> vResult(nMu);
    vResult.clear();
    for (Int_t iMu = 0; iMu < nMu; ++iMu) {
      if (muP4[iMu].Pt() > 20. && TMath::Abs(muP4[iMu].Eta()) < 2.4 && muRelIso[iMu] < 0.15 && muTrkLayers[iMu] >= 5) {
        vResult.emplace_back(iMu);
      }
    }
    return vResult;
  }, {"nMu", "muP4", "muRelIso", "muTrkLayers"})
  .Define("nMuPassBasicCuts", "static_cast<Int_t>(muIdxPassBasicCuts.size())")
  .Redefine("nMu", "nMuPassBasicCuts")
  .Define("THINjetIdxPassBasicCuts", [](const Int_t THINnJet, const ROOT::RVec<TypeLorentzVector> THINjetP4
      , const ROOT::RVec<TypeLorentzVector> eleP4, const ROOT::RVec<TypeLorentzVector> muP4) {
    const auto lepP4 = ROOT::VecOps::Concatenate(eleP4, muP4);
    ROOT::RVec<Int_t> vResult(THINnJet);
    vResult.clear();
    for (Int_t iJet = 0; iJet < THINnJet; ++iJet) {
      const auto &p4THIN = THINjetP4[iJet];
      if (
          p4THIN.Pt() > 30. && TMath::Abs(p4THIN.Eta()) < 2.5
          && ROOT::VecOps::All(ROOT::VecOps::Map(
              lepP4,
              [&p4THIN](const TypeLorentzVector p4Lep) {
                return GetDeltaR2(p4THIN, p4Lep) > 0.4 * 0.4;
          }))
      ) {
        vResult.emplace_back(iJet);
      }
    }
    return vResult;
  }, {"THINnJet", "THINjetP4", "eleP4", "muP4"})
  .Define("THINnJetPassBasicCuts", "static_cast<Int_t>(THINjetIdxPassBasicCuts.size())")
  .Redefine("THINnJet", "THINnJetPassBasicCuts")
  .Define("FATjetIdxPassBasicCuts", [](const Int_t FATnJet, const ROOT::RVec<TypeLorentzVector> FATjetP4
      , const ROOT::RVec<TypeLorentzVector> eleP4, const ROOT::RVec<TypeLorentzVector> muP4) {
    const auto lepP4 = ROOT::VecOps::Concatenate(eleP4, muP4);
    ROOT::RVec<Int_t> vResult(FATnJet);
    vResult.clear();
    for (Int_t iJet = 0; iJet < FATnJet; ++iJet) {
      const auto &p4FAT = FATjetP4[iJet];
      if (p4FAT.Pt() > 170. && TMath::Abs(p4FAT.Eta()) < 2.5
          && ROOT::VecOps::All(ROOT::VecOps::Map(
              lepP4,
              [&p4FAT](const TypeLorentzVector p4Lep) {
                return GetDeltaR2(p4FAT, p4Lep) > 0.8 * 0.8;
          }))
      ) {
        vResult.emplace_back(iJet);
      }
    }
    return vResult;
  }, {"FATnJet", "FATjetP4", "eleP4", "muP4"})
  .Define("FATnJetPassBasicCuts", "static_cast<Int_t>(FATjetIdxPassBasicCuts.size())")
  .Redefine("FATnJet", "FATnJetPassBasicCuts")
  ;
  for (std::string nameCol: {"nEle", "nMu", "THINnJet", "FATnJet"}) {
    vNameColOriginal.emplace_back(nameCol + "PassBasicCuts");
  }
  {
    auto vNamesCol = dfOriginal.GetColumnNames();
    std::vector<std::string> vTypenamesCol(vNamesCol.size());
    vTypenamesCol.clear();
    for (const auto nameCol: vNamesCol){
      vTypenamesCol.emplace_back(dfOriginal.GetColumnType(nameCol));
      // if (debug) std::cerr << vTypenamesCol.back() << ",\tsize: " << vTypenamesCol.size() << std::endl;
    }
    for (size_t iCol = 0; iCol < vNamesCol.size(); ++iCol) {
      std::string &nameCol = vNamesCol[iCol];
      std::string &typenameCol = vTypenamesCol[iCol];
      std::regex r("^is(.+)Muon$");
      std::smatch m;
      if (std::regex_match(nameCol, m, r) && typenameCol == "ROOT::VecOps::RVec<bool>") {
        const std::string nameAlias = "muIsPass" + std::string(m[1]);
        if (debug) std::cerr << "nameAlias: " << nameAlias << std::endl;
        dfOriginal = dfOriginal.Define(nameAlias, nameCol);
        nameCol = nameAlias;
      }
    }
    for (const std::string pref: {"ele", "mu", "THINjet", "FATjet"}) {
      RedefinePrefWithIdx(dfOriginal, pref, {}, pref + "IdxPassBasicCuts",
        [debug, &vNameColOriginal, &aavNameColHasJet, &aavNameColLPairPassPt, &aavNameColMatching, &aavNameColAllMatched]
        (const std::string nameCol, const Int_t nDims, const std::string typenameE){
          const TString tstrTypenameE = typenameE;
          if (nDims == 1 && !tstrTypenameE.Contains("Vector") && !tstrTypenameE.Contains("4D") && !tstrTypenameE.Contains("3D") && !tstrTypenameE.Contains("2D")
            && ( tstrTypenameE.Contains("double")
              || tstrTypenameE.Contains("float")
              || tstrTypenameE.Contains("long")
              || tstrTypenameE.Contains("int")
              || tstrTypenameE.Contains("short")
              || tstrTypenameE.Contains("char")
              // || tstrTypenameE.Contains("bool")
            )
          ) {
            if (debug) std::cerr << "Pick " << nameCol << std::endl;
            vNameColOriginal.emplace_back(nameCol);
            for (auto paav: {&aavNameColHasJet, &aavNameColLPairPassPt, &aavNameColMatching, &aavNameColAllMatched}) {
              for (auto &av: *paav) {
                for (auto &v: av) {
                  v.emplace_back(nameCol);
                }
              }
            }
          }
        }, debug);
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
            const auto funGen = [](const ROOT::RVec<Bool_t> rvId)->Int_t{
             return ROOT::VecOps::Sum(rvId, 0);
            };
            if (dfOriginal.HasColumn(nameColNew)) {
              dfOriginal = dfOriginal
              .Redefine(nameColNew, funGen, {{ nameCol }});
            } else {
              dfOriginal = dfOriginal
              .Define(nameColNew, funGen, {{ nameCol }});
            }
            vNameColOriginal.emplace_back(nameColNew);
            for (auto paav: {&aavNameColHasJet, &aavNameColLPairPassPt, &aavNameColMatching, &aavNameColAllMatched}) {
              for (auto &av: *paav) {
                for (auto &v: av) {
                  v.emplace_back(nameColNew);
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
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    dfOriginal = dfOriginal
    .Define(aPrefLepFlavLower[iLepFlav] + "PairedIdx", [&debug](
        const Int_t nLep,
        const ROOT::RVec<Bool_t> isLoose, const ROOT::RVec<Bool_t> isMedium, const ROOT::RVec<Bool_t> isTight,
        const ROOT::RVec<TypeLorentzVector> &lepP4)->ROOT::RVec<Int_t>{
        // if (debug) Info("lambda:PairedIdx", "nLep: %d, (%zu, %zu, %zu)", nLep, isLoose.size(), isMedium.size(), isTight.size());
        ROOT::RVec<Int_t> vResult(2, -1);
        if (nLep < 2) {
          return vResult;
        }
        for (vResult[0] = 0; vResult[0] < nLep - 1; ++(vResult[0])) {
          for (vResult[1] = vResult[0] + 1; vResult[1] < nLep; ++(vResult[1])) {
            if (ROOT::VecOps::All(ROOT::VecOps::Take(isMedium, vResult))) {
              return vResult;
            }
          }
        }
        std::fill(vResult.begin(), vResult.end(), -1);
        return vResult;
      }, {"n" + aPrefLepFlav[iLepFlav],
          aPrefLepFlavLower[iLepFlav] + "IsPassLoose", aPrefLepFlavLower[iLepFlav] + "IsPassMedium", aPrefLepFlavLower[iLepFlav] + "IsPassTight",
          aPrefLepFlavLower[iLepFlav] + "P4"});
  }
  // Lazily register histogram action for the Original stage
  {
    const size_t nCol = vNameColOriginal.size();
    vHistViewOriginal.clear();
    vHistViewOriginal.reserve(nCol + 2);
    vHistViewOriginal.emplace_back(GetHistFromColumn(dfOriginal, "mcWeight"));
    vHistViewOriginal.emplace_back(GetHistFromColumn(dfOriginal, "mcWeightSgn"));
    for (const std::string &nameCol: vNameColOriginal) {
      vHistViewOriginal.emplace_back(GetHistFromColumn(dfOriginal, nameCol, "mcWeightSgn"));
    }
  }
  // Begin the Gen stages in case the input dataset is MC Signal
  std::array<ROOT::RDF::RNode, 2> aDfGen {dfOriginal, dfOriginal};
  if (isSignal) {
    for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
      aDfGen[iLepFlav] = aDfGen[iLepFlav]
      .Filter("ROOT::VecOps::Any("
        "genParId==" + std::to_string(iLepFlav ? pdgMuon : pdgElectron) +
        "&& genMomParId==" + std::to_string(pdgZ) +
        ")");
      avNameColGen[iLepFlav] = vNameColOriginal;
    }
    // Lazily register histogram action for the Gen stages
    for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
      const size_t nCol = avNameColGen[iLepFlav].size();
      avHistViewGen[iLepFlav].clear();
      avHistViewGen[iLepFlav].reserve(nCol + 2);
      avHistViewGen[iLepFlav].emplace_back(GetHistFromColumn(aDfGen[iLepFlav], "mcWeight"));
      avHistViewGen[iLepFlav].emplace_back(GetHistFromColumn(aDfGen[iLepFlav], "mcWeightSgn"));
      for (const std::string &nameCol: avNameColGen[iLepFlav]) {
        avHistViewGen[iLepFlav].emplace_back(GetHistFromColumn(aDfGen[iLepFlav], nameCol, "mcWeightSgn"));
      }
    }
  }
  // Begin the HasLPair stages
  std::array<ROOT::RDF::RNode, 2> aDfHasLPair = aDfGen;
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    aDfHasLPair[iLepFlav] = aDfHasLPair[iLepFlav]
    .Filter((aPrefLepFlavLower[iLepFlav] + "PairedIdx[1]!=-1").c_str());
    // if (debug) {
    //   std::cerr << "Columns in dfHasLPair[" << iLepFlav << "]: " << std::flush;
    //   for (const auto nameCol: aDfHasLPair[iLepFlav].GetColumnNames()) {
    //     std::cerr << nameCol << ", " << std::flush;
    //   }
    //   std::cerr << "}" << std::endl;
    // }
    for (const auto nameCol : aDfHasLPair[iLepFlav].GetColumnNames()) {
      const TString tstrNameCol(nameCol);
      const std::string typenameCol = aDfHasLPair[iLepFlav].GetColumnType(nameCol);
      const TString tstrTypenameCol = typenameCol;
      if (tstrNameCol.BeginsWith(aPrefLepFlavLower[iLepFlav].c_str()) && tstrTypenameCol.BeginsWith("ROOT::VecOps::RVec") && ! tstrNameCol.BeginsWith(aPrefLepFlavLower[iLepFlav] + "PairedIdx")) {
        const std::string nameColNew = aPrefLepFlavLower[iLepFlav] + "Paired" + tstrNameCol(aPrefLepFlavLower[iLepFlav].length(), tstrNameCol.Length()).Data();
        if (aDfHasLPair[iLepFlav].HasColumn(nameColNew)) continue;
        std::string typenameE = "";
        const Int_t nDims = RefgetE2D(typenameE, typenameCol) ? 2 : (RefgetE1D(typenameE, typenameCol) ? 1 : 0);
        const TString tstrTypenameE = typenameE;
        if (nDims) {
          if (debug) std::cerr << "Defining " << nameColNew << std::endl;
          aDfHasLPair[iLepFlav] = aDfHasLPair[iLepFlav]
          .Define(nameColNew, "ROOT::VecOps::Take(" + nameCol + ", " + aPrefLepFlavLower[iLepFlav] + "PairedIdx)");
        }
        if (nDims == 1
            && !tstrTypenameE.Contains("Vector") && !tstrTypenameE.Contains("4D") && !tstrTypenameE.Contains("3D") && !tstrTypenameE.Contains("2D")
            && (tstrTypenameE.Contains("double")
                || tstrTypenameE.Contains("float")
                || tstrTypenameE.Contains("long")
                || tstrTypenameE.Contains("short")
                || tstrTypenameE.Contains("int")
                || tstrTypenameE.Contains("char")
                // || tstrTypenameE.Contains("bool")
            )
        ) {
          avNameColHasLPair[iLepFlav].emplace_back(nameColNew);
          for (auto paav: {&aavNameColHasJet, &aavNameColLPairPassPt, &aavNameColMatching, &aavNameColAllMatched}) {
            for (auto &v: (*paav)[iLepFlav]) {
              v.emplace_back(nameColNew);
            }
          }
        }
      }
    }
    aDfHasLPair[iLepFlav] = aDfHasLPair[iLepFlav]
    .Define(
      aPrefLepFlavLower[iLepFlav] + "PairP4",
      "ROOT::VecOps::Sum("
      + aPrefLepFlavLower[iLepFlav] + "PairedP4"
      ", " + typenameLorentzVector + "())");
    for (const std::string suff: {"Pt", "Eta", "Phi", "E", "M", "Mt2"}) {
      const std::string nameColNew = aPrefLepFlavLower[iLepFlav] + "Pair" + suff;
      if (debug) std::cerr << "Defining " << nameColNew << std::endl;
      aDfHasLPair[iLepFlav] = aDfHasLPair[iLepFlav]
      .Define(
        nameColNew,
        aPrefLepFlavLower[iLepFlav] + "PairP4." + suff + "()");
      for (auto pav: {&avNameColHasLPair, &avNameColZMassCutted}) {
        (*pav)[iLepFlav].emplace_back(nameColNew);
      }
      for (auto paav: {&aavNameColHasJet, &aavNameColLPairPassPt, &aavNameColMatching, &aavNameColAllMatched}) {
        for (auto &v: (*paav)[iLepFlav]) {
          v.emplace_back(nameColNew);
        }
      }
    }
    aDfHasLPair[iLepFlav] = aDfHasLPair[iLepFlav]
    .Define(aPrefLepFlavLower[iLepFlav] + "PairIsPassZ",
      "TMath::Abs(" + aPrefLepFlavLower[iLepFlav] + "PairM" + " - " + std::to_string(massZ) + ") < 15."
      // " && " + aPrefLepFlavLower[iLepFlav] + "PairPt > 130."
    );
    for (auto pav: {&avNameColHasLPair, &avNameColZMassCutted}) {
      (*pav)[iLepFlav].emplace_back(aPrefLepFlavLower[iLepFlav] + "PairIsPassZ");
    }
    for (auto paav: {&aavNameColHasJet, &aavNameColLPairPassPt, &aavNameColMatching, &aavNameColAllMatched}) {
      for (auto &v: (*paav)[iLepFlav]) {
        v.emplace_back(aPrefLepFlavLower[iLepFlav] + "PairIsPassZ");
      }
    }
    aDfHasLPair[iLepFlav] = aDfHasLPair[iLepFlav]
    .Define("HPSTauIdxPassBasicCuts", [](
      const ROOT::RVec<TypeLorentzVector> &HPSTauP4, const ROOT::RVec<TypeLorentzVector> &lepPairedP4,
      const ROOT::RVec<Bool_t> disc_decayModeFindingNewDMs, const ROOT::RVec<Bool_t> disc_byVTightIsolationMVA3newDMwLT) {
      const Int_t nHPSTau = HPSTauP4.size();
      ROOT::RVec<size_t> vResult(nHPSTau);
      vResult.clear();
      for (Int_t iHPSTau = 0; iHPSTau < nHPSTau; ++iHPSTau) {
        if (HPSTauP4[iHPSTau].Pt() > 18. && TMath::Abs(HPSTauP4[iHPSTau].Eta()) < 2.5
          && disc_decayModeFindingNewDMs[iHPSTau] && disc_byVTightIsolationMVA3newDMwLT[iHPSTau]) {
          Bool_t isPass = true;
          for (const auto& p4Lep: lepPairedP4) {
            if (GetDeltaR2(HPSTauP4[iHPSTau], p4Lep) < 0.4 * 0.4) {
              isPass = false;
              break;
            }
          }
          if (isPass) vResult.emplace_back(iHPSTau);
        }
      }
      return vResult;
    }, {"HPSTauP4", aPrefLepFlavLower[iLepFlav] + "PairedP4", "disc_decayModeFindingNewDMs", "disc_byVTightIsolationMVA3newDMwLT"})
    .Define("nHPSTauPassBasicCuts", "static_cast<Int_t>(HPSTauIdxPassBasicCuts.size())")
    .Define("nHPSTau", "nHPSTauPassBasicCuts");
    RedefinePrefWithIdx(aDfHasLPair[iLepFlav], "HPSTau", {}, "HPSTauIdxPassBasicCuts",
      [debug, &avNameColHasLPair]
      (const std::string nameCol, const Int_t nDims, const std::string typenameE){
        const TString tstrTypenameE = typenameE;
        if (nDims == 1 && !tstrTypenameE.Contains("Vector") && !tstrTypenameE.Contains("4D") && !tstrTypenameE.Contains("3D") && !tstrTypenameE.Contains("2D")
          && ( tstrTypenameE.Contains("double")
            || tstrTypenameE.Contains("float")
            || tstrTypenameE.Contains("long")
            || tstrTypenameE.Contains("int")
            || tstrTypenameE.Contains("short")
            || tstrTypenameE.Contains("char")
            // || tstrTypenameE.Contains("bool")
          )
        ) {
          if (debug) std::cerr << "Pick " << nameCol << std::endl;
          for (auto pav: { &avNameColHasLPair }) {
            for (auto &v: *pav) {
              v.emplace_back(nameCol);
            }
          }
        }
      }, debug);
  }
  // Lazily register histogram action for the HasLPair stages
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    const size_t nCol = avNameColHasLPair[iLepFlav].size();
    avHistViewHasLPair[iLepFlav].clear();
    avHistViewHasLPair[iLepFlav].reserve(nCol + 2);
    avHistViewHasLPair[iLepFlav].emplace_back(GetHistFromColumn(aDfHasLPair[iLepFlav], "mcWeight"));
    avHistViewHasLPair[iLepFlav].emplace_back(GetHistFromColumn(aDfHasLPair[iLepFlav], "mcWeightSgn"));
    for (const std::string &nameCol: avNameColHasLPair[iLepFlav]) {
      avHistViewHasLPair[iLepFlav].emplace_back(GetHistFromColumn(aDfHasLPair[iLepFlav], nameCol, "mcWeightSgn"));
    }
  }
  // Begin the HasVtx stages
  std::array<ROOT::RDF::RNode, 2> aDfHasVtx = aDfHasLPair;
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    aDfHasVtx[iLepFlav] = aDfHasVtx[iLepFlav]
    .Filter("nVtx>0");
  }
  // Lazily register histogram action for the HasVtx stages
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    const size_t nCol = avNameColHasVtx[iLepFlav].size();
    avHistViewHasVtx[iLepFlav].clear();
    avHistViewHasVtx[iLepFlav].reserve(nCol + 2);

    if (debug) std::cerr << "Specifying hist for mcWeight and mcWeightSgn" << std::endl;
    avHistViewHasVtx[iLepFlav].emplace_back(GetHistFromColumn(aDfHasVtx[iLepFlav], "mcWeight"));
    avHistViewHasVtx[iLepFlav].emplace_back(GetHistFromColumn(aDfHasVtx[iLepFlav], "mcWeightSgn"));
    for (const std::string &nameCol: avNameColHasVtx[iLepFlav]) {
      if (debug) std::cerr << "Specifying hist for column " << nameCol << std::endl;
      avHistViewHasVtx[iLepFlav].emplace_back(GetHistFromColumn(aDfHasVtx[iLepFlav], nameCol, "mcWeightSgn"));
    }
  }
  // Begin the LPairPassPt stages
  std::array<ROOT::RDF::RNode, 2> aDfLPairedPassPt = aDfHasVtx;
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    aDfLPairedPassPt[iLepFlav] = aDfLPairedPassPt[iLepFlav]
    .Filter([](const ROOT::RVec<TypeLorentzVector> &lepPairedP4) {
      return lepPairedP4[0].Pt() > 25. && lepPairedP4[1].Pt() > 20.;
    }, {aPrefLepFlavLower[iLepFlav] + "PairedP4"});
  }
  // Lazily register histogram action for the LPairedPassPt stage
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    const size_t nCol = avNameColLPairedPassPt[iLepFlav].size();
    avHistViewLPairedPassPt[iLepFlav].clear();
    avHistViewLPairedPassPt[iLepFlav].reserve(nCol + 2);
    avHistViewLPairedPassPt[iLepFlav].emplace_back(GetHistFromColumn(aDfLPairedPassPt[iLepFlav], "mcWeight"));
    avHistViewLPairedPassPt[iLepFlav].emplace_back(GetHistFromColumn(aDfLPairedPassPt[iLepFlav], "mcWeightSgn"));
    for (const std::string &nameCol: avNameColLPairedPassPt[iLepFlav]) {
      avHistViewLPairedPassPt[iLepFlav].emplace_back(GetHistFromColumn(aDfLPairedPassPt[iLepFlav], nameCol, "mcWeightSgn"));
    }
  }
  // Begin the ZMassCutted stages
  std::array<ROOT::RDF::RNode, 2> aDfZMassCutted = aDfLPairedPassPt;
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    aDfZMassCutted[iLepFlav] = aDfZMassCutted[iLepFlav].Filter(aPrefLepFlavLower[iLepFlav] + "PairIsPassZ");
  }
  // Lazily register histogram action for the ZMassCutted stages
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    const size_t nCol = avNameColZMassCutted[iLepFlav].size();
    avHistViewZMassCutted[iLepFlav].clear();
    avHistViewZMassCutted[iLepFlav].reserve(nCol + 2);
    avHistViewZMassCutted[iLepFlav].emplace_back(GetHistFromColumn(aDfZMassCutted[iLepFlav], "mcWeight"));
    avHistViewZMassCutted[iLepFlav].emplace_back(GetHistFromColumn(aDfZMassCutted[iLepFlav], "mcWeightSgn"));
    for (const std::string &nameCol: avNameColZMassCutted[iLepFlav]) {
      avHistViewZMassCutted[iLepFlav].emplace_back(GetHistFromColumn(aDfZMassCutted[iLepFlav], nameCol, "mcWeightSgn"));
    }
  }
  // Begin the NoExtraL stage
  std::array<ROOT::RDF::RNode, 2> aDfNoExtraL = aDfZMassCutted;
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    aDfNoExtraL[iLepFlav] = aDfNoExtraL[iLepFlav]
    .Filter(Form("nElePassLoose<=%d&&nMuPassSoft<=%d", iLepFlav ? 0 : 2, iLepFlav ? 2 : 0));
  }
  // Lazily register histogram action for the NoExtraL stages
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    const size_t nCol = avNameColNoExtraL[iLepFlav].size();
    avHistViewNoExtraL[iLepFlav].clear();
    avHistViewNoExtraL[iLepFlav].reserve(nCol + 2);
    avHistViewZMassCutted[iLepFlav].emplace_back(GetHistFromColumn(aDfNoExtraL[iLepFlav], "mcWeight"));
    avHistViewZMassCutted[iLepFlav].emplace_back(GetHistFromColumn(aDfNoExtraL[iLepFlav], "mcWeightSgn"));
    for (const std::string &nameCol: avNameColNoExtraL[iLepFlav]) {
      avHistViewNoExtraL[iLepFlav].emplace_back(GetHistFromColumn(aDfNoExtraL[iLepFlav], nameCol, "mcWeightSgn"));
    }
  }
  // Begin the NoTau stage
  std::array<ROOT::RDF::RNode, 2> aDfNoTau = aDfNoExtraL;
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    aDfNoTau[iLepFlav] = aDfNoTau[iLepFlav]
    .Filter("nHPSTau==0");
  }
  // Lazily register histogram action for the NoTau stages
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    const size_t nCol = avNameColNoTau[iLepFlav].size();
    avHistViewNoTau[iLepFlav].clear();
    avHistViewNoTau[iLepFlav].reserve(nCol + 2);
    avHistViewNoTau[iLepFlav].emplace_back(GetHistFromColumn(aDfNoTau[iLepFlav], "mcWeight"));
    avHistViewNoTau[iLepFlav].emplace_back(GetHistFromColumn(aDfNoTau[iLepFlav], "mcWeightSgn"));
    for (const std::string &nameCol: avNameColNoTau[iLepFlav]) {
      avHistViewNoTau[iLepFlav].emplace_back(GetHistFromColumn(aDfNoTau[iLepFlav], nameCol, "mcWeightSgn"));
    }
  }
  // Begin the HasJet stage
  // https://stackoverflow.com/questions/28541488/nested-aggregate-initialization-of-stdarray
  std::array<std::array<ROOT::RDF::RNode, 2>, 2> aaDfHasJet {{ {aDfNoTau[0], aDfNoTau[0]}, {aDfNoTau[1], aDfNoTau[1]} }};
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      aaDfHasJet[iLepFlav][iAK] = aaDfHasJet[iLepFlav][iAK].Filter(aPrefAKShort[iAK] + "nJet" + " >= 2");
    }
  }
  // Lazily register histogram action for the HasJet stages
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      const size_t nCol = aavHistViewHasJet[iLepFlav][iAK].size();
      aavHistViewHasJet[iLepFlav][iAK].clear();
      aavHistViewHasJet[iLepFlav][iAK].reserve(nCol + 2);
      aavHistViewHasJet[iLepFlav][iAK].emplace_back(GetHistFromColumn(aaDfHasJet[iLepFlav][iAK], "mcWeight"));
      aavHistViewHasJet[iLepFlav][iAK].emplace_back(GetHistFromColumn(aaDfHasJet[iLepFlav][iAK], "mcWeightSgn"));
      for (const std::string &nameCol: aavNameColHasJet[iLepFlav][iAK]) {
        aavHistViewHasJet[iLepFlav][iAK].emplace_back(GetHistFromColumn(aaDfHasJet[iLepFlav][iAK], nameCol, "mcWeightSgn"));
      }
    }
  }
  // Begin the LPairPassPt stage
  std::array<std::array<ROOT::RDF::RNode, 2>, 2> aaDfLPairPassPt = aaDfHasJet;
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      aaDfLPairPassPt[iLepFlav][iAK] = aaDfLPairPassPt[iLepFlav][iAK].Filter(aPrefLepFlavLower[iLepFlav] + "PairP4.Pt() >= 50.");
    }
  }
  // Lazily register histogram action for the LPairPassPt stages
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      const size_t nCol = aavHistViewLPairPassPt[iLepFlav][iAK].size();
      aavHistViewLPairPassPt[iLepFlav][iAK].clear();
      aavHistViewLPairPassPt[iLepFlav][iAK].reserve(nCol + 2);
      aavHistViewLPairPassPt[iLepFlav][iAK].emplace_back(GetHistFromColumn(aaDfLPairPassPt[iLepFlav][iAK], "mcWeight"));
      aavHistViewLPairPassPt[iLepFlav][iAK].emplace_back(GetHistFromColumn(aaDfLPairPassPt[iLepFlav][iAK], "mcWeightSgn"));
      for (const std::string &nameCol: aavNameColLPairPassPt[iLepFlav][iAK]) {
        aavHistViewLPairPassPt[iLepFlav][iAK].emplace_back(GetHistFromColumn(aaDfLPairPassPt[iLepFlav][iAK], nameCol, "mcWeightSgn"));
      }
    }
  }
  // Realize the actions and write to output files
  TFile* tfOut = TFile::Open(fileOut.c_str(), "recreate");
  // tfOut->mkdir("allEventsCounter", tfIn->Get<TDirectory>("allEventsCounter")->GetTitle())->cd();
  // TH1 *histTotalEvents = tfIn->Get<TH1>("allEventsCounter/totalEvents");
  // if (isSignal) {
  //   histTotalEvents->Scale(1. / 3);
  // }
  // histTotalEvents->Write();
  tfOut->cd("/");
  tfOut->mkdir("Original", "Unfiltered entries");
  if (isSignal) {
    for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
      tfOut->mkdir(("Gen" + aPrefLepFlav[iLepFlav]).c_str(), ("GEN-level " + aPrefLepFlavLower[iLepFlav] + " events").c_str());
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->mkdir(("HasLPair" + aPrefLepFlav[iLepFlav]).c_str(), ("Entries with correct " + aPrefLepFlavLower[iLepFlav] + " numbers").c_str());
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->mkdir(("HasVtx" + aPrefLepFlav[iLepFlav]).c_str(), "Entries with at least 1 vertex");
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->mkdir(("LPairedPassPt" + aPrefLepFlav[iLepFlav]).c_str(), Form("Entries with paired %s passing pT cuts", aPrefLepFlav[iLepFlav].c_str()));
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->mkdir(("ZMassCutted" + aPrefLepFlav[iLepFlav]).c_str(), "Entries passing the Z mass cut");
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->mkdir(("NoExtraL" + aPrefLepFlav[iLepFlav]).c_str(), "Entries with no extra leptons");
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->mkdir(("NoTau" + aPrefLepFlav[iLepFlav]).c_str(), "Entries with no cutted HPSTau");
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      tfOut->mkdir(("HasJet" + aPrefLepFlav[iLepFlav] + aPrefAKShort[iAK] + "jet").c_str(), ("Entries with at least 1 " + aPrefAKShort[iAK] + " jet").c_str());
    }
  }
  // for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
  //   for (size_t iAK = 0; iAK < 2; ++iAK) {
  //     tfOut->mkdir(("LPairPassPt" + aPrefLepFlav[iLepFlav] + aPrefAKShort[iAK] + "jet").c_str(), ("Entries with " + aPrefLepFlav[iLepFlav] + " Pair pT < 50").c_str());
  //   }
  // }
  tfOut->Close();
  std::vector<ROOT::RDF::RResultPtr<ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager>>> vSn;
  vSn.clear();
  {
    const std::vector<std::string> vNameColJetCmp ({
      "nGenPar"//, "genParPt", "genParEta", "genParPhi", "genParE", "genParM", "genParId", "genMomParId", "genParIndex", "genParQ", "genParSt"
      // , "nEle", "elePt", "eleEta", "elePhi", "eleE", "eleM", "eleCharge", "eleChargeConsistent"
      , "elePairedIdx"// , "elePairedPt", "elePairedEta", "elePairedPhi", "elePairedE", "elePairedM", "elePairedCharge", "elePairedChargeConsistent"
      // , "nMu", "muPt", "muEta", "muPhi", "muE", "muM", "muCharge"
      , "muPairedIdx"// , "muPairedPt", "muPairedEta", "muPairedPhi", "muPairedE", "muPairedM", "muPairedCharge"
      , "THINnJet", "THINjetIdxPassBasicCuts"// , "THINjetPt", "THINjetEta", "THINjetPhi", "THINjetE", "THINjetM", "THINjetArea", "THINjetCharge"
      , "FATnJet", "FATjetIdxPassBasicCuts"// , "FATjetPt", "FATjetEta", "FATjetPhi", "FATjetE", "FATjetM", "FATjetArea", "FATjetCharge"
      // , "pfMetCorrPt", "pfMetCorrEta", "pfMetCorrPhi", "pfMetCorrSumEt", "pfMetCorrSig", "pfMetCorrUnc"
    });
    for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
      for (size_t iAK = 0; iAK < 2; ++iAK) {
        vSn.emplace_back(aaDfHasJet[iLepFlav][iAK].Snapshot("HasJet" + aPrefLepFlav[iLepFlav] + aPrefAKShort[iAK] + "jet/tree", fileOut, vNameColJetCmp, {"update", ROOT::kZLIB, 1, false, 99, true, true}));
        // vSn.emplace_back(aaDfHasJet[iLepFlav][iAK].Snapshot("HasJet" + aPrefLepFlav[iLepFlav] + aPrefAKShort[iAK] + "jet/tree", fileOut, aavNameColHasJet[iLepFlav][iAK], {"update", ROOT::kZLIB, 1, false, 99, true, true}));
      }
    }
  }
  for (auto &&sn: vSn) {
    sn.GetValue();
  }
  tfOut = TFile::Open(fileOut.c_str(), "update");
  tfOut->cd("/");
  tfOut->cd("Original");
  for (auto &&histView: vHistViewOriginal) {
    histView->Write();
  }
  if (isSignal) {
    for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
      tfOut->cd("/");
      tfOut->cd(("Gen" + aPrefLepFlav[iLepFlav]).c_str());
      for (auto &&histView: avHistViewGen[iLepFlav]) {
        histView->Write();
      }
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->cd("/");
    tfOut->cd(("HasLPair" + aPrefLepFlav[iLepFlav]).c_str());
    for (auto &&histView: avHistViewHasLPair[iLepFlav]) {
      histView->Write();
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->cd("/");
    tfOut->cd(("HasVtx" + aPrefLepFlav[iLepFlav]).c_str());
    for (auto &&histView: avHistViewHasVtx[iLepFlav]) {
      histView->Write();
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->cd("/");
    tfOut->cd(("LPairedPassPt" + aPrefLepFlav[iLepFlav]).c_str());
    for (auto &&histView: avHistViewLPairedPassPt[iLepFlav]) {
      histView->Write();
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->cd("/");
    tfOut->cd(("ZMassCutted" + aPrefLepFlav[iLepFlav]).c_str());
    for (auto &&histView: avHistViewZMassCutted[iLepFlav]) {
      histView->Write();
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->cd("/");
    tfOut->cd(("NoExtraL" + aPrefLepFlav[iLepFlav]).c_str());
    for (auto &&histView: avHistViewNoExtraL[iLepFlav]) {
      histView->Write();
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    tfOut->cd("/");
    tfOut->cd(("NoTau" + aPrefLepFlav[iLepFlav]).c_str());
    for (auto &&histView: avHistViewNoTau[iLepFlav]) {
      histView->Write();
    }
  }
  for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
    for (size_t iAK = 0; iAK < 2; ++iAK) {
      tfOut->cd("/");
      tfOut->cd(("HasJet" + aPrefLepFlav[iLepFlav] + aPrefAKShort[iAK] + "jet").c_str());
      for (auto &&histView: aavHistViewHasJet[iLepFlav][iAK]) {
        histView->Write();
      }
    }
  }
  // for (size_t iLepFlav = 0; iLepFlav < 2; ++iLepFlav) {
  //   for (size_t iAK = 0; iAK < 2; ++iAK) {
  //     tfOut->cd("/");
  //     tfOut->cd(("LPairPassPt" + aPrefLepFlav[iLepFlav] + aPrefAKShort[iAK] + "jet").c_str());
  //     for (auto &&histView: aavHistViewLPairPassPt[iLepFlav][iAK]) {
  //       histView->Write();
  //     }
  //   }
  // }
  tfOut->Close();
  if (debug) std::cerr << "Completed!" << std::endl;
  tfIn->Close();
}

#if true
int main(int argc, char** argv) {
  xAna_monoZ_preselect(argv[1], argv[2], true, true);
  return 0;
}
#endif
