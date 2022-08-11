#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDF/HistoModels.hxx>
#include <Math/Vector4D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>
#include <Math/RootFinder.h>
#include <Math/RootFinderAlgorithms.h>
#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TError.h>

#include <string>
#include <cstring>
#include <regex>

#define OPTPARSE_IMPLEMENTATION
#define OPTPARSE_API static
#include "skeeto_optparse.h"

template <class E>
ROOT::RVec<E> UniqueSort(const ROOT::RVec<E> &v)
{
  ROOT::RVec<E> r(v);
  std::sort(r.begin(), r.end());
  typename ROOT::RVec<E>::iterator &&rIterEnd = std::unique(r.begin(), r.end());
  r.resize(std::distance(r.begin(), rIterEnd));
  return r;
}

template <class E, class Compare, class Predicate>
ROOT::RVec<E> UniqueSort(const ROOT::RVec<E> &v, Compare &&c, Predicate &&p)
{
  ROOT::RVec<E> r(v);
  std::sort(r.begin(), r.end(), std::forward<Compare>(c));
  typename ROOT::RVec<E>::iterator &&rIterEnd = std::unique(r.begin(), r.end(), std::forward<Predicate>(p));
  r.resize(std::distance(r.begin(), rIterEnd));
  return r;
}

template <class E, class Compare>
ROOT::RVec<E> UniqueSort(const ROOT::RVec<E> &v, const Compare &c)
{
  ROOT::RVec<E> r(v);
  std::sort(r.begin(), r.end(), c);
  typename ROOT::RVec<E>::iterator &&rIterEnd = std::unique(r.begin(), r.end(), [&c](const E &a, const E &b)
                                                            { return !c(a, b) && !c(b, a); });
  r.resize(std::distance(r.begin(), rIterEnd));
  return r;
}

template <class E>
void sort_uniquely(std::vector<E> &v)
{
  std::sort(v.begin(), v.end());
  typename std::vector<E>::iterator &&vIterEnd = std::unique(v.begin(), v.end());
  v.resize(std::distance(v.begin(), vIterEnd));
}

template <class E>
ROOT::RVec<E> UniqueFirst(const ROOT::RVec<E> &v)
{
  using size_type = typename ROOT::RVec<E>::size_type;
  ROOT::RVec<size_type> vIdx(ROOT::VecOps::StableArgsort(v));
  typename ROOT::RVec<size_type>::iterator &&vIdxIterEnd = std::unique(
      vIdx.begin(), vIdx.end(),
      [&vIdx](const size_type ia, const size_type ib)
      {
        return vIdx[ia] == vIdx[ib];
      });
  vIdx.resize(std::distance(vIdx.begin(), vIdxIterEnd));
  std::sort(vIdx.begin(), vIdx.end());
  ROOT::RVec<E> vResult(vIdx.size());
  for (const size_type &&i : vIdx)
  {
    vResult.emplace_back(vIdx[i]);
  }
  return vResult;
}

// https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
bool ReplaceStringFirst(std::string &str, const std::string &from, const std::string &to)
{
  const size_t iStart = str.find(from);
  if (iStart == str.npos)
  {
    return false;
  }
  str.replace(iStart, from.length(), to);
  return true;
}

std::string GetBasename(const std::string &path)
{
  if (!path.size())
  {
    return "";
  }
  std::regex r("^(?:.*/)?([^/]+)/?");
  std::smatch m;
  if (std::regex_match(path, m, r))
  {
    return m[1];
  }
  return "";
}

template <class Vector1, class Vector2>
inline typename Vector1::Scalar DetXY(const Vector1 &v1, const Vector2 &v2)
{
  return v1.X() * v2.Y() - v1.Y() * v2.X();
}

/**Generate a histogram of the specified expression lazily
 * from an RDataFram
 * by applying Histo1D() with specified binning information
 *
 * This overload does the actual application
 * and store the binning specification
 * in the histogram title in JSON format
 */
template <typename V = ROOT::RDFDetail::RInferredType, typename W = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
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
    const std::string exprWeight = "")
{
  const Int_t nBins = (upperLimitBins - lowerLimitBins);
  const Double_t binWidth = TMath::Power(10., -binDensityOrder);
  const Double_t lowerLimit = alignment + binWidth * lowerLimitBins;
  const Double_t upperLimit = alignment + binWidth * upperLimitBins;
  const char *&&cstrJSON = Form("{\"name\":\"%s\",\"typename\":\"%s\",\"expression\":\"%s\",%s\"binDensityOrder\":%d,\"alignment\":%F,\"isLowerAssigned\":%s,\"lowerLimitBins\":%d,\"isUpperAssigned\":%s,\"upperLimitBins\":%d}",
                                nameColumn.c_str(), typenameColumn.c_str(), expression.c_str(), exprWeight.length() ? ("\"exprWeight\":\"" + exprWeight + "\",").c_str() : "",
                                binDensityOrder, alignment, isLowerAssigned ? "true" : "false", lowerLimitBins, isUpperAssigned ? "true" : "false", upperLimitBins);
  const ROOT::RDF::TH1DModel model{("h" + nameColumn).c_str(), cstrJSON,
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
template <typename V = ROOT::RDFDetail::RInferredType, typename W = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumnCustom(D &df, const std::string nameColumn, const std::string typenameColumn, const std::string expression, std::string exprWeight = "", const std::string nameColumnStripped = "")
{
  if (!nameColumnStripped.length())
  {
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
  if (nameColumn == "Counter")
  {
    binDensityOrder = 0;
    alignment = 0;
    isLowerAssigned = true;
    lowerLimitBins = 0;
    isUpperAssigned = true;
    upperLimitBins = 1;
  }
  else if (tstrTypenameColumn.Contains("Bool") || tstrTypenameColumn.Contains("bool"))
  {
    binDensityOrder = 0;
    alignment = 0.;
    isLowerAssigned = true;
    lowerLimitBins = 0;
    isUpperAssigned = true;
    upperLimitBins = 2;
  }
  else if (tstrTypenameColumn.Contains("Int") || tstrTypenameColumn.Contains("int") || ((tstrTypenameColumn.Contains("short") || tstrTypenameColumn.Contains("long")) && !(tstrTypenameColumn.Contains("float") || tstrTypenameColumn.Contains("double"))))
  {
    binDensityOrder = 0;
    alignment = -0.5;
    if (tstrNameColumnStripped.Contains("Idx") || tstrNameColumnStripped.Contains("Rank"))
    {
      isLowerAssigned = true;
      lowerLimitBins = 0;
      if (tstrNameColumnStripped.Contains("THINjet"))
      {
        isUpperAssigned = true;
        upperLimitBins = 10;
      }
      else if (tstrNameColumnStripped.Contains("FATjet"))
      {
        isUpperAssigned = true;
        upperLimitBins = 10;
      }
    }
    else if (
        (tstrNameColumnStripped.BeginsWith("n") && !tstrNameColumnStripped.Contains("PuppinJet")) || tstrNameColumnStripped.Contains("nJet") || tstrNameColumnStripped.Contains("nEle") || tstrNameColumnStripped.Contains("nMu") || tstrNameColumnStripped.Contains("nPho") || tstrNameColumnStripped.Contains("nGen") || tstrNameColumnStripped.Contains("_n"))
    {
      isLowerAssigned = true;
      lowerLimitBins = 0;
      // nGenDPair 0, 1, 2
      if (tstrNameColumnStripped.BeginsWith("nGenDPair") || (tstrNameColumnStripped.EndsWith("JetClosestToGenDPair")))
      {
        isUpperAssigned = true;
        upperLimitBins = 3;
        // nGenD 0, 1, 2, 3, 4
      }
      else if (tstrNameColumnStripped.BeginsWith("nGenD") || (tstrNameColumnStripped.EndsWith("JetClosestToGenD")))
      {
        isUpperAssigned = true;
        upperLimitBins = 5;
      }
    }
    else if (tstrNameColumnStripped.EndsWith("Sgn"))
    {
      isLowerAssigned = true;
      lowerLimitBins = -1;
      isUpperAssigned = true;
      upperLimitBins = 2;
    }
  }
  else
  {
    alignment = 0.;
    if (nameColumn == "mcWeight")
    {
      binDensityOrder = 2;
    }
    else if (tstrNameColumnStripped.EndsWith("LogDiff"))
    {
      binDensityOrder = 1;
      isLowerAssigned = true;
      lowerLimitBins = -60;
      isUpperAssigned = true;
      upperLimitBins = 60;
    }
    else if (tstrNameColumnStripped.EndsWith("Pt"))
    {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    }
    else if (tstrNameColumnStripped.EndsWith("Rho"))
    {
      binDensityOrder = 1;
      isLowerAssigned = true;
      lowerLimitBins = 0;
    }
    else if (tstrNameColumnStripped.EndsWith("Eta"))
    {
      binDensityOrder = 1;
      upperLimitBins = 100;
      lowerLimitBins = -upperLimitBins;
    }
    else if (tstrNameColumnStripped.EndsWith("Phi"))
    {
      binDensityOrder = 2;
      isUpperAssigned = true;
      upperLimitBins = TMath::CeilNint(TMath::Pi() * TMath::Power(10., binDensityOrder));
      isLowerAssigned = true;
      lowerLimitBins = -upperLimitBins;
    }
    else if (tstrNameColumnStripped.EndsWith("E") || tstrNameColumnStripped.EndsWith("Et"))
    {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    }
    else if (tstrNameColumnStripped.EndsWith("M") || tstrNameColumnStripped.EndsWith("M2") || tstrNameColumnStripped.EndsWith("MTTwo") || tstrNameColumnStripped.EndsWith("MTTwo2"))
    {
      isLowerAssigned = true;
      lowerLimitBins = 0;
    }
    else if (tstrNameColumnStripped.EndsWith("DeltaR") || tstrNameColumnStripped.EndsWith("DeltaR2"))
    {
      isLowerAssigned = true;
      lowerLimitBins = 0;
      binDensityOrder = 2;
    }
    else if (tstrNameColumnStripped.EndsWith("Sig"))
    {
      binDensityOrder = 2;
      isLowerAssigned = true;
      lowerLimitBins = 0;
    }
    else if (tstrNameColumnStripped.EndsWith("jetArea"))
    {
      binDensityOrder = 2;
      isLowerAssigned = true;
      lowerLimitBins = 0;
    }
    else if (tstrNameColumnStripped.EndsWith("EF"))
    {
      binDensityOrder = 2;
      isLowerAssigned = true;
      lowerLimitBins = 0;
      isUpperAssigned = true;
      upperLimitBins = TMath::Power(10, binDensityOrder);
    }
    else if (tstrNameColumnStripped(0, tstrNameColumnStripped.Length() - 2).String().EndsWith("Cov"))
    {
      binDensityOrder = 2;
      isLowerAssigned = true;
      lowerLimitBins = 0;
      isUpperAssigned = true;
      upperLimitBins = TMath::Power(10, binDensityOrder);
    }
    else if (tstrNameColumnStripped.Contains("Btag"))
    {
      binDensityOrder = 2;
    }
    else if (tstrNameColumnStripped.Contains("Ctag"))
    {
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
template <typename V = ROOT::RDFDetail::RInferredType, typename W = ROOT::RDFDetail::RInferredType, typename D = ROOT::RDF::RNode>
ROOT::RDF::RResultPtr<TH1D> GetHistFromColumn(D &df, const std::string nameColumn, std::string exprWeight = "", const std::string nameColumnStripped = "")
{
  const std::string typenameColumn = df.GetColumnType(nameColumn);
  return GetHistFromColumnCustom<V, W, D>(df, nameColumn, typenameColumn, nameColumn, exprWeight, nameColumnStripped);
}

/**Get the value of a parameter from the file name string
 *
 * E.g. when keypattern=="Mx1", filename=="Mx2-1_Mv-500_Mx1-0p1", val will be "0p1"
 */
Bool_t RefgetParamFilename(std::string &val, const std::string filename, const std::string keypattern)
{
  std::regex r("^(?:.*?[._/])?" + keypattern + "-([^._/]+).*$");
  std::smatch m;
  const Bool_t result = std::regex_match(filename, m, r);
  if (result)
  {
    val = m[1];
  }
  return result;
}

/** Get the value of a parameter from the file name string
 * with "p" substituted with "."
 *
 * E.g. when keypattern=="Mx1", filename=="Mx2-1_Mv-500_Mx1-0p1", val will be "0.1"
 */
Bool_t RefgetParamFilenameNum(std::string &val, const std::string filename, const std::string keypattern)
{
  const Bool_t result = RefgetParamFilename(val, filename, keypattern);
  if (result)
  {
    ReplaceStringFirst(val, "p", ".");
  }
  return result;
}

/**Determine if a type name string describes a 2D Rvec
 * and get the element type name string.
 *
 * E.g.
 * "ROOT::VecOps::RVec<ROOT::VecOps::RVec<Double_t> >" to "Double_t"
 * "ROOT::VecOps::RVec<vector<LorentzVector<PtEtaPhiM4D<double> > > >" to LorentzVector<PtEtaPhiM4D<double> >
 *
 * Regex pattern "(?:xxxx)" represents a non-capture group that will not be listed as a match result
 * while the normal "(xxxx)" will be listed
 * If matched, m[0] will be the whole string, and the group result starts from m[1]
 *
 * "\\s" matches to a white-space-like character such as " " and "\t",
 * and "\\S" matches to a character other than "\\s" does
 */
Bool_t RefgetE2D(std::string &typenameE, const std::string typenameCol)
{
  std::regex r("^ROOT::VecOps::RVec\\s*<\\s*(?:vector|ROOT::VecOps::RVec)\\s*<\\s*(.*\\S)\\s*>\\s*>\\s*$");
  std::smatch m;
  const Bool_t result = (std::regex_match(typenameCol, m, r));
  if (result)
  {
    typenameE = m[1];
  }
  return result;
}

/**Determine if a type name string describes a 1D Rvec
 * and get the element type name string.
 */
Bool_t RefgetE1D(std::string &typenameE, const std::string typenameCol)
{
  std::regex r("^ROOT::VecOps::RVec\\s*<\\s*(.*\\S)\\s*>\\s*$");
  std::smatch m;
  const Bool_t result = (std::regex_match(typenameCol, m, r));
  if (result)
  {
    typenameE = m[1];
  }
  return result;
}

inline std::string GetExprTakeIdx2D(const std::string nameCol, const std::string nameIdx, const std::string typenameE, const int debug = 0)
{
  return
      // "if (true && " + nameCol + ".size() <= ROOT::VecOps::Max(" + nameIdx + ")) {"
      // "  Fatal(\"lambda:RedefineWithIdx\", \"" + nameIdx + " maximum (%d) exceeds the " + nameCol + " size (%zu)\", ROOT::VecOps::Max(" + nameIdx + "), " + nameCol + ".size());"
      // "}"
      "ROOT::VecOps::RVec<ROOT::VecOps::RVec<" + typenameE + ">> vvResult(" + nameCol + ".size());"
                                                                                        "vvResult.clear();"
                                                                                        "for (size_t iVE: " +
      nameIdx + ") {"
                "  const auto& vE = " +
      nameCol + "[iVE];"
                "  ROOT::VecOps::RVec<" +
      typenameE + "> resultSub(vE.size());"
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

inline std::string GetExprTakeIdx1D(const std::string nameCol, const std::string nameIdx, const std::string typenameE, const Bool_t debug = 0)
{
  return Form("ROOT::VecOps::Take(%s, %s)", nameCol.c_str(), nameIdx.c_str());
}

template <class D = ROOT::RDF::RNode>
void RedefinePrefWithIdx(D &df, const std::string pref, const std::vector<std::string> vPrefExcl, const std::string nameIdx, const std::function<void(const std::string nameCol, const Int_t nDims, const std::string typenameE)> funPickNameCol = nullptr, const int debug = 0)
{
  for (const std::string &nameCol : df.GetColumnNames())
  {
    if (std::find(vPrefExcl.cbegin(), vPrefExcl.cend(), pref) != vPrefExcl.cend())
      continue;
    std::string typenameCol = df.GetColumnType(nameCol);
    const TString tstrNameCol = nameCol;
    std::string typenameE = "";
    const Int_t nDims = RefgetE2D(typenameE, typenameCol) ? 2 : (RefgetE1D(typenameE, typenameCol) ? 1 : 0);
    const TString tstrTypenameE = typenameE;
    if (tstrNameCol.BeginsWith(pref))
    {
      if (nDims == 2)
      {
        if (debug)
          std::cerr << "Redefining 2D: " << nameCol << " (" << typenameE << ") (" << typenameCol << ")" << std::endl;
        df = df
                 .Redefine(nameCol, GetExprTakeIdx2D(nameCol, nameIdx, typenameE, debug));
      }
      else if (nDims == 1)
      {
        if (debug)
          std::cerr << "Redefining 1D: " << nameCol << " (" << typenameE << ") (" << typenameCol << ")" << std::endl;
        df = df
                 .Redefine(nameCol, GetExprTakeIdx1D(nameCol, nameIdx, typenameE, debug));
      }
      if (funPickNameCol && nDims > 0)
      {
        funPickNameCol(nameCol, nDims, typenameE);
      }
    }
  }
}

void DefineP4Component(ROOT::RDF::RNode &df, const std::string prefNameCol)
{
  const std::string nameCol = prefNameCol + "P4";
  for (const auto suf : {"Pt", "Eta", "Phi", "E", "Et", "M"})
  {
    df = df.Define(prefNameCol + suf, nameCol + "." + suf + "()");
  }
}

typedef ROOT::Math::PtEtaPhiMVector TypeLorentzVector;

void DefineP4Components1D(ROOT::RDF::RNode &df, const std::string prefNameCol)
{
  const std::string nameCol = prefNameCol + "P4";
  df = df
           .Define(prefNameCol + "Pt", [](const ROOT::RVec<TypeLorentzVector> &vP4)
                   { return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4)
                                              { return p4.Pt(); }); },
                   {{nameCol}})
           .Define(prefNameCol + "Eta", [](const ROOT::RVec<TypeLorentzVector> &vP4)
                   { return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4)
                                              { return p4.Eta(); }); },
                   {{nameCol}})
           .Define(prefNameCol + "Phi", [](const ROOT::RVec<TypeLorentzVector> &vP4)
                   { return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4)
                                              { return p4.Phi(); }); },
                   {{nameCol}})
           .Define(prefNameCol + "E", [](const ROOT::RVec<TypeLorentzVector> &vP4)
                   { return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4)
                                              { return p4.E(); }); },
                   {{nameCol}})
           .Define(prefNameCol + "Et", [](const ROOT::RVec<TypeLorentzVector> &vP4)
                   { return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4)
                                              { return p4.Et(); }); },
                   {{nameCol}})
           .Define(prefNameCol + "M", [](const ROOT::RVec<TypeLorentzVector> &vP4)
                   { return ROOT::VecOps::Map(vP4, [](const TypeLorentzVector &p4)
                                              { return p4.M(); }); },
                   {{nameCol}});
}

Double_t GetMTTwo(const TypeLorentzVector p4DDA, const TypeLorentzVector p4DDB, const Double_t ptMet, const Double_t phiMet)
{
  const ROOT::Math::XYVector p2DDA(p4DDA.X(), p4DDA.Y());
  const ROOT::Math::XYVector p2DDB(p4DDB.X(), p4DDB.Y());
  const ROOT::Math::XYVector basep2Met(TMath::Cos(phiMet), TMath::Sin(phiMet));
  const ROOT::Math::XYVector basep2MetPerp(-basep2Met.Y(), basep2Met.X());
  const ROOT::Math::XYVector p2DDARotated(p2DDA.Dot(basep2Met), p2DDA.Dot(basep2MetPerp));
  const ROOT::Math::XYVector p2DDBRotated(p2DDB.Dot(basep2Met), p2DDB.Dot(basep2MetPerp));
  if (p4DDA.M2() >= p4DDB.M2() + (p4DDB.Et() - p2DDBRotated.X()) * ptMet * 2)
  {
    return p4DDA.M();
  }
  if (p4DDB.M2() >= p4DDA.M2() + (p4DDA.Et() + p2DDBRotated.X()) * ptMet * 2)
  {
    return p4DDB.M();
  }
  const Double_t ptX1aInit =
      (p4DDA.M2() - p4DDB.M2() + p2DDBRotated.X() * ptMet - p4DDB.Et() * ptMet * 2) / (-(p4DDA.Et() + p4DDB.Et()) * 2 + p2DDARotated.X() + p2DDBRotated.X());
  Double_t pParaX1a = ptX1aInit;
  Double_t pVertX1a = 0.;
  // ROOT::Math::XYVector gradMtSqA = basep2Met * (p4DDA.Et() * 2) - p2DDA;
  // ROOT::Math::XYVector gradMtSqB = basep2Met * (p4DDB.Et() * 2) - p2DDB;
  // ROOT::Math::XYVector direction2D(-gradMtSqA.X()+gradMtSqA.Y(), gradMtSqB.X()-gradMtSqB.Y());
  ROOT::Math::XYVector gradMtSqARotated(ROOT::Math::XYVector(1., 0.) - p2DDARotated);
  ROOT::Math::XYVector gradMtSqBRotated(ROOT::Math::XYVector(1., 0) - p2DDBRotated);
  ROOT::Math::XYVector direction2DRotated(-gradMtSqARotated.X() + gradMtSqARotated.Y(), gradMtSqBRotated.X() - gradMtSqBRotated.Y());
  Double_t slopeIntersection = DetXY(gradMtSqARotated, gradMtSqBRotated);
  slopeIntersection /= direction2DRotated.R();
  direction2DRotated /= direction2DRotated.R();
  const Double_t sgnInit = std::signbit(slopeIntersection) ? -1 : 1;
  constexpr Double_t step = 0.1;
  direction2DRotated *= -sgnInit;
  slopeIntersection *= -sgnInit;
  std::function<Double_t(Double_t pParaX1a)> fDeltaMT2Para = [
                                                                 // mutable
                                                                 &pVertX1a,
                                                                 // constant
                                                                 &p4DDA, &p4DDB, &p2DDARotated, &p2DDBRotated, &ptMet](Double_t pParaX1a) -> Double_t
  {
    return (p4DDA.M2() + p4DDA.Et() * (pParaX1a * pParaX1a + pVertX1a * pVertX1a) * 2 - (p2DDARotated.X() * pParaX1a + p2DDARotated.Y() * pVertX1a) * 2) - (p4DDB.M2() + p4DDB.Et() * ((ptMet - pParaX1a) * (ptMet - pParaX1a) + pVertX1a * pVertX1a) * 2 - (p2DDBRotated.X() * (ptMet - pParaX1a) - p2DDBRotated.Y() * pVertX1a) * 2);
  };
  std::function<Double_t(Double_t pParaX1a)> fGradDeltaMT2Para = [
                                                                     // mutable
                                                                     &pVertX1a,
                                                                     // constant
                                                                     &p4DDA, &p4DDB, &p2DDARotated, &p2DDBRotated, &ptMet](Double_t pParaX1a) -> Double_t
  {
    return (
               (p4DDA.Et() * pParaX1a / (pParaX1a * pParaX1a + pVertX1a * pVertX1a) - p2DDARotated.X()) - (p4DDB.Et() * (ptMet - pParaX1a) / ((ptMet - pParaX1a) * (ptMet - pParaX1a) + pVertX1a * pVertX1a) - p2DDBRotated.X())) *
           2;
  };
  ROOT::Math::XYVector p2X1aRotated(ptX1aInit, 0), p2X1bRotated(ptMet - ptX1aInit, 0);
  ROOT::Math::XYVector p2X1aRotatedPrev(p2X1aRotated);
  Double_t slopeIntersectionPrev = slopeIntersection;
  while (slopeIntersection < 0 && -slopeIntersection <= __FLT_EPSILON__)
  {
    pVertX1a += direction2DRotated.Y() * step;
    Double_t pParaX1aPred = pParaX1a + direction2DRotated.X() * step;
    auto *finder = new ROOT::Math::RootFinder();
    finder->SetMethod(ROOT::Math::RootFinder::kBRENT);
    ROOT::Math::Functor1D f(fDeltaMT2Para);
    ROOT::Math::GradFunctor1D g(fDeltaMT2Para, fGradDeltaMT2Para);
    finder->SetFunction(f, -step * 100, step * 100);
    finder->SetFunction(g, pParaX1aPred);
    finder->Solve();
    pParaX1a = finder->Root();
    delete finder;
    p2X1aRotated.SetXY(pParaX1a, pVertX1a);
    p2X1bRotated.SetXY(ptMet - pParaX1a, pVertX1a);
    gradMtSqARotated = p2X1aRotated * (p4DDA.Et() * 2 / p2X1aRotated.R()) - p2DDARotated * 2;
    gradMtSqBRotated = p2X1bRotated * (p4DDB.Et() * 2 / p2X1bRotated.R()) - p2DDBRotated * 2;
    direction2DRotated.SetXY((-gradMtSqARotated.X() + gradMtSqARotated.Y()) * (-sgnInit), (gradMtSqBRotated.X() - gradMtSqBRotated.Y()) * (-sgnInit));
    slopeIntersection = DetXY(gradMtSqARotated, gradMtSqBRotated) * (-sgnInit);
    slopeIntersection /= direction2DRotated.R();
    direction2DRotated /= direction2DRotated.R();
    p2X1aRotatedPrev = p2X1aRotated;
    slopeIntersectionPrev = slopeIntersection;
  }
  if (TMath::Abs(slopeIntersection) > __FLT_EPSILON__)
  {
    p2X1aRotated =
        p2X1aRotatedPrev * (-slopeIntersection / (slopeIntersection - slopeIntersectionPrev)) + p2X1bRotated * (slopeIntersectionPrev / (slopeIntersection - slopeIntersectionPrev));
  }
  return TMath::Sqrt(p4DDA.M2() + (p4DDA.Et() * p2X1aRotated.R() - p2DDARotated.Dot(p2X1aRotated)) * 2);
}

void xAna_monoZ_inconservation(const std::vector<std::string> vFileIn, const std::string fileOut, const size_t nThread = 1, const int debug = 0)
{
  const std::string typenameLorentzVector = "ROOT::Math::PtEtaPhiMVector";
  ROOT::EnableImplicitMT(nThread);
  constexpr Double_t massZ = 91.1876; //< static mass of Z (constant)
  constexpr Double_t massElectron =
      0.0005109989461; //< static mass of electron (constant)
  constexpr Double_t massDown = 0.0048;
  constexpr Int_t pdgZ = 23;       //< pdgid of Z
  constexpr Int_t pdgZp = 55;      //< pdgid of Z'
  constexpr Int_t pdgX2 = 18;      //< pdgid of x2
  constexpr Int_t pdgX1 = 5000522; //< pdgid of x1
  constexpr Int_t pdgDown = 1;     //< pdgid of down quark
  constexpr Int_t pdgElectron = 11;
  constexpr Int_t pdgMuon = 13;
  constexpr Int_t pdgTau = 15;
  const Bool_t isSignal = TString(fileOut).Contains("ignal");
  Double_t mX2Truth = 100.;
  {
    std::string strMX2;
    if (RefgetParamFilenameNum(strMX2, vFileIn[0], "Mchi2"))
    {
      mX2Truth = std::stod(strMX2);
    }
  }
  std::vector<ROOT::RDF::RResultPtr<TH1D>> vHistViewGenUnion = {};
  ROOT::RDF::RNode df = ROOT::RDataFrame("tree/treeMaker", vFileIn);
  df = df
           .Define("mcWeightSgn", "static_cast<Int_t>((mcWeight > 0) - (0 > mcWeight))");
  df = df
           .Define("genZpIdx", [&pdgZp](const ROOT::RVec<Int_t> genParId, const ROOT::RVec<Int_t> genMomParId) -> Int_t
                   {
      // std::cerr << "genZpIdx: Calculating ..." << std::flush;
      // auto &&result = static_cast<Int_t>(ROOT::VecOps::Nonzero(genParId == pdgZp && genMomParId /*== 10002*/!= pdgZp)[0]);
      // std::cerr << " Done, result: " << result << std::endl;
      // return result;
      const auto vResult = ROOT::VecOps::Nonzero(genParId == pdgZp && genMomParId == 10002);
      return vResult.size() ? static_cast<Int_t>(vResult[0]) : -1; },
                   {"genParId", "genMomParId"})
           .Define("genX2Idx", [](const Int_t nGenPar, const ROOT::RVec<Int_t> genParId, const ROOT::RVec<Int_t> genMomParId)
                   {
      ROOT::RVec<Int_t> result(2, -1);
      for (int i = 0; i < nGenPar; ++i) {
        if (TMath::Abs(genParId[i]) == pdgX2 && TMath::Abs(genMomParId[i]) == pdgZp) {
          result[(genParId[i] < 0)] = i;
        }
      }
      return result; },
                   {"nGenPar", "genParId", "genMomParId"})
           .Define("genX1Idx", [](const Int_t nGenPar, const ROOT::RVec<Int_t> genParId, const ROOT::RVec<Int_t> genMomParId)
                   {
      ROOT::RVec<Int_t> result(2, -1);
      for (int i = 0; i < nGenPar; ++i) {
        if (TMath::Abs(genParId[i]) == pdgX1 && TMath::Abs(genMomParId[i]) == pdgX2) {
          result[(genMomParId[i] < 0)] = i;
        }
      }
      return result; },
                   {"nGenPar", "genParId", "genMomParId"})
           .Define("genDIdx", [](const Int_t nGenPar, const ROOT::RVec<Int_t> genParId, const ROOT::RVec<Int_t> genMomParId)
                   {
      ROOT::RVec<Int_t> result(4, -1);
      for (int i = 0; i < nGenPar; ++i) {
        if (TMath::Abs(genParId[i]) == pdgDown && TMath::Abs(genMomParId[i]) == pdgX2) {
          result[((genMomParId[i] < 0) << 1) + (genParId[i] < 0)] = i;
        }
      }
      return result; },
                   {"nGenPar", "genParId", "genMomParId"});
  for (const char *parName : {"Zp", "X2"})
  {
    for (const char *suff : {"Px", "Py", "Pz", "E"})
    {
      df = df
               .Define(Form("gen%s%s", parName, suff), Form("ROOT::VecOps::Take(genPar%s, gen%sIdx)", suff, parName));
    }
  }
  for (const char *suff : {"Px", "Py", "Pz", "E"})
  {
    df = df
             .Define(Form("genZp%sChildDiff", suff), Form("genZp%s - ROOT::VecOps::Sum(genX2%s)", suff, suff));
    vHistViewGenUnion.emplace_back(df.Histo1D(Form("genZp%sChildDiff", suff), "mcWeightSgn"));
  }
  TFile *tfOut = TFile::Open(fileOut.c_str(), "RECREATE");
  tfOut->mkdir("GenUnion")->cd();
  for (auto &histview: vHistViewGenUnion) {
    histview->Write();
  }
  tfOut->Close();
}

#if true
int main(int argc, char **argv)
{
  size_t nThread = 1;
  int debug = 0;
  bool useMultiple = false;
  struct optparse_long longopts[] = {
      {"help", 'h', OPTPARSE_NONE},
      {"multiple", 'm', OPTPARSE_NONE},
      {"debug", 'v', OPTPARSE_NONE},
      {"threads", 'j', OPTPARSE_REQUIRED},
      {NULL} // NULL termination
  };
  int option;
  struct optparse options;
  optparse_init(&options, argv);
  // argv is NULL terminated, so no need to check argc here
  while ((option = optparse_long(&options, longopts, NULL)) != -1)
  {
    switch (option)
    {
    case 'h':
      std::cout << "Description:\n"
                   "xAna_monoZ_preselect [-v] [-j n] pathFileOut pathFileInGlob1 [pathFileInGlob2 ...]\n"
                   "  -h --help\tDisplay this help message.\n"
                   "  -v --debug\tRun with debug message and examinations.\n"
                   "\tThis can be specified multiple times.\n"
                   "  -j --threads\tSpecify the number of threads to use by EnableImplicitMT.\n"
                   "\t0 means the number of all the logical cores.\n"
                   "\tDefault to 1.\n"
                << std::endl;
      return 0;
      break;
    case 'v':
      ++debug;
      break;
    case 'j':
      nThread = std::atoi(options.optarg);
      break;
    case '?':
      std::fprintf(stderr, "%s: %s\n", argv[0], options.errmsg);
      std::cerr << std::flush;
      return 1;
      break;
    }
  }
  if (debug)
    std::cerr << "Debug level: " << debug << std::endl;
  if (debug)
    std::cerr << "Number of remaining command-line arguments: " << argc - options.optind << std::endl;
  if (argc - options.optind < 2)
  {
    std::cerr << "error: Expect pathFileOut pathFileInGlob1 pathFileInGlob2 ..." << std::endl;
    return 1;
  }
  const std::string pathFileOut(optparse_arg(&options));
  if (debug)
    std::cerr << "pathFileOut: " << pathFileOut << std::endl;
  const Int_t nFileIn = argc - options.optind;
  if (debug)
    std::cerr << "nFileIn: " << nFileIn;
  std::vector<std::string> vPathFileIn(nFileIn, "");
  for (std::string &pathFileIn : vPathFileIn)
  {
    pathFileIn = optparse_arg(&options);
    if (debug)
      std::cerr << "pathFileIn: " << pathFileIn << std::endl;
  }
  xAna_monoZ_inconservation(vPathFileIn, pathFileOut, nThread, debug);
  return 0;
}
#endif
