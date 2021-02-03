#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
// #include <TH1F.h>
#include <TLeaf.h>
#include <TList.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector.h>
#include <TVectorFfwd.h>
// #include <TClonesArray.h>
// #include <TObjArray.h>
#include <TMath.h>
#include <TSystem.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <vector>

#include "HistMerger.h"

#ifndef INTPOW_FUNCTION
#define INTPOW_FUNCTION
template <typename TBase, typename TIndex, typename TRatio>
TRatio intpow(const TBase base, const TIndex index, const TRatio ratio) {
  if (index == 0) {
    return ratio;
  }
  TRatio result = ratio;
  if (index > 0) {
    for (TIndex i = 0; i < index; i++) {
      result *= base;
    }
  } else {
    for (TIndex i = 0; i < -index; i++) {
      result /= base;
    }
  }
  return result;
}
#endif
#ifndef GETNDIGITSM1_FUNCTION
#define GETNDIGITSM1_FUNCTION
Int_t getNDigitsM1(Double_t x) {
  return (Int_t)TMath::Floor(TMath::Log10(TMath::Abs(x)));
}
#endif

// #ifndef TSTRING_UTILS
// #define TSTRING_UTILS
// Bool_t startsWithOfTString(const TString target, const TString pattern,
//                            const Ssiz_t displacement = 0) {
//   return target.SubString(pattern).Start() == displacement;
// }

// Bool_t endsWithOfTString(const TString target, const TString pattern, const
// Ssiz_t displacement=0) {
//   return target.SubString(pattern).Start() == target.Length() -
//   pattern.Length() - displacement;
// }

// Bool_t containsOfTString(const TString target, const TString pattern) {
//   return target.SubString(pattern).Start() != -1;
// }
// #endif

#ifndef HISTMERGER_C
#define HISTMERGER_C

class HistMerger::LeafAnalyzerDefault : public LeafAnalyzerAbstract {
 public:
  static inline TString GetNameLeafModified(TString nameLeaf) {
    // TString nameLeafNew = TString(nameLeaf);
    Ssiz_t indexLeftSquareBracket = nameLeaf.First('[');
    // if (indexLeftSquareBracket != -1) {
    //   nameLeafNew.Resize(indexLeftSquareBracket);
    // }
    // return nameLeafNew;
    return indexLeftSquareBracket == -1 ? nameLeaf
                                        : nameLeaf(0, indexLeftSquareBracket);
  }
  // static inline TString GetNameLeafModified(TLeaf *leaf) {
  //   return GetNameLeafModified((TString)leaf->GetName());
  // }
  static inline TString GetTitleLeaf(TLeaf *leaf) {
    TString titleLeaf = leaf->GetTitle();
    if (titleLeaf == "") {
      titleLeaf = leaf->GetBranch()->GetTitle();
    }
    return titleLeaf;
  }
  static inline TString GetTypeNameLeaf(TLeaf *leaf) {
    return leaf->GetTypeName();
  }

 protected:
  Bool_t debug;
  TString nameTT;
  std::function<Bool_t(Int_t &NBinCorrect, Double_t &lowerCorrect,
                       Double_t &upperCorrect, TString nameTT,
                       TString nameLeafModified, TString typeNameLeaf,
                       TString titleLeaf)>
      assignHistSettingPerLeafTreeExtra;
  std::function<void(Int_t &NBinCorrect, Double_t &lowerCorrect,
                     Double_t &upperCorrect, TString nameTT,
                     TString nameLeafModified, TString typeNameLeaf,
                     TString titleLeaf)>
      adjustHistSettingPerLeafTreeExtra;
  Bool_t dontCheckEmptyness;
  std::vector<TString> vNameLeafFile;
  TString nameLeafModified;
  TString titleLeaf;
  TString typeNameLeaf;
  std::function<TString(TLeaf *leaf)> funTitleLeaf;
  Bool_t isHistSettingAssignedDirectly;
  TString expressionBeforeSetting;
  TCut selection;
  Option_t *optionDraw;
  Long64_t nEntriesToDrawMax;
  Long64_t firstEntryToDraw;
  std::vector<TString> vDeps;
  std::function<Bool_t(TTree *tree)> funHasTargetExtra;
  Bool_t isEverAnalyzed;
  Bool_t allowNeverAnalyzed;
  std::vector<Bool_t> vIsEmpty;
  Bool_t areAllEmpty;
  Int_t nBinsCorrect;
  Double_t lowerCorrect;
  Double_t upperCorrect;
  TString tstrHistSetting;

  enum class HistTypeFlavor { FloatingPoint, Integer, Boolean };
  HistTypeFlavor histTypeFlavor;
  // Bool_t isTypeFloatingPoint;
  Bool_t isFirstNonEmptyLeaf;
  std::vector<Int_t> vNEntriesFile;
  // std::vector<Double_t> vMinimumFile;
  Double_t minimumValueAllFiles;
  // std::vector<Double_t> vMaximumFile;
  Double_t maximumValueAllFiles;
  std::vector<Int_t> vNBinsFile;
  std::vector<Double_t> vLowerFile;
  std::vector<Double_t> vUpperFile;

 public:
  LeafAnalyzerDefault() {
    this->vNameLeafFile.clear();
    this->nameLeafModified = "";
    this->titleLeaf = "";
    this->typeNameLeaf = "";
    this->funTitleLeaf = nullptr;
    this->isHistSettingAssignedDirectly = false;
    this->isEverAnalyzed = false;
    this->allowNeverAnalyzed = true;
    this->vIsEmpty.clear();
    this->areAllEmpty = true;
    this->nBinsCorrect = 0;
    this->lowerCorrect = 0.;
    this->upperCorrect = 0.;
    this->debug = false;
    this->nameTT = "";
    this->assignHistSettingPerLeafTreeExtra = nullptr;
    this->assignHistSettingPerLeafTreeExtra = nullptr;
    this->adjustHistSettingPerLeafTreeExtra = nullptr;
    this->dontCheckEmptyness = false;
    this->expressionBeforeSetting = "";
    this->selection = "";
    this->optionDraw = "";
    this->nEntriesToDrawMax = TTree::kMaxEntries;
    this->firstEntryToDraw = 0;
    this->vDeps.clear();
    this->funHasTargetExtra = nullptr;
    histTypeFlavor = HistTypeFlavor::FloatingPoint;
    // this->isTypeFloatingPoint = false;
    this->isFirstNonEmptyLeaf = true;
    this->vNEntriesFile.clear();
    // this->vMinimumFile.clear();
    // this->vMaximumFile.clear();
    this->minimumValueAllFiles = 0;
    this->maximumValueAllFiles = 0;
    this->vNBinsFile.clear();
    this->vLowerFile.clear();
    this->vUpperFile.clear();
    this->tstrHistSetting = "";
  }
  void SetDebug(Bool_t debug) { this->debug = debug; }
  void SetNameTT(TString nameTT) { this->nameTT = nameTT; }
  TString GetNameTT() const { return this->nameTT; }
  void SetFunAssignHistSettingExtra(
      std::function<Bool_t(Int_t &NBinCorrect, Double_t &lowerCorrect,
                           Double_t &upperCorrect, TString nameTT,
                           TString nameLeafModified, TString typeNameLeaf,
                           TString titleLeaf)>
          assignHistSettingPerLeafTreeExtra) {
    this->assignHistSettingPerLeafTreeExtra = assignHistSettingPerLeafTreeExtra;
  }
  void SetFunAdjustHistSettingExtra(
      std::function<void(Int_t &NBinCorrect, Double_t &lowerCorrect,
                         Double_t &upperCorrect, TString nameTT,
                         TString nameLeafModified, TString typeNameLeaf,
                         TString titleLeaf)>
          adjustHistSettingPerLeafTreeExtra) {
    this->adjustHistSettingPerLeafTreeExtra = adjustHistSettingPerLeafTreeExtra;
  }
  virtual void SetExpressionCustom(TString name, TString typeName,
                                   TString title,
                                   TString expressionBeforeSetting,
                                   TCut selection = "", Option_t *option = "",
                                   Long64_t nentries = TTree::kMaxEntries,
                                   Long64_t firstentry = 0) {
    this->nameLeafModified = name;
    this->typeNameLeaf = typeName;
    this->titleLeaf = title;
    this->expressionBeforeSetting = expressionBeforeSetting;
    this->selection = selection;
    this->optionDraw = option;
    this->nEntriesToDrawMax = nentries;
    this->firstEntryToDraw = firstentry;
  }
  void SetExpressionBeforeSetting(TString expressionBeforeSetting) {
    this->expressionBeforeSetting = expressionBeforeSetting;
  }
  void SetSelection(TCut selection) { this->selection = selection; }
  void SetOptionDraw(Option_t *option) { this->optionDraw = option; }
  void SetNEntriesToDrawMax(Long64_t nentries) {
    this->nEntriesToDrawMax = nentries;
  }
  void SetFirstEntryToDraw(Long64_t firstentry) {
    this->firstEntryToDraw = firstentry;
  }
  TString GetExpressionBeforeSetting() const { return expressionBeforeSetting; }
  virtual void SetHasTarget(
      std::vector<TString> vDeps,
      std::function<Bool_t(TTree *tree)> funHasTargetExtra = nullptr) {
    this->vDeps = vDeps;
    this->funHasTargetExtra = funHasTargetExtra;
  }
  virtual void SetAllowNeverAnalyzed(Bool_t allowNeverAnalyzed) {
    this->allowNeverAnalyzed = allowNeverAnalyzed;
  }
  virtual Bool_t GetAllowNeverAnalyzed() const { return allowNeverAnalyzed; }
  virtual Bool_t GetIsEverAnalyzed() const { return isEverAnalyzed; }
  virtual void AssignHistSetting(Int_t nBinsCorrect, Double_t lowerCorrect,
                                 Double_t upperCorrect,
                                 TString tstrHistSetting = "") {
    this->nBinsCorrect = nBinsCorrect;
    this->lowerCorrect = lowerCorrect;
    this->upperCorrect = upperCorrect;
    if (tstrHistSetting.Length()) {
      this->tstrHistSetting = tstrHistSetting;
    }
    this->isHistSettingAssignedDirectly = true;
  }
  void SetDontCheckEmptyness(Bool_t dontCheckEmptyness) {
    this->dontCheckEmptyness = dontCheckEmptyness;
  }
  void SetFunTitleLeaf(std::function<TString(TLeaf *leaf)> funTitleLeaf) {
    this->funTitleLeaf = funTitleLeaf;
  }
  Bool_t GetDontCheckEmptyness() const {
    return this->isHistSettingAssignedDirectly && this->dontCheckEmptyness;
  }
  void AnalyzeLeafBasic(TLeaf *leaf) {
    TString nameLeaf = leaf->GetName();
    this->nameLeafModified = this->GetNameLeafModified(nameLeaf);
    if (debug) std::cout << "bla" << std::endl;
    this->titleLeaf =
        funTitleLeaf == nullptr ? GetTitleLeaf(leaf) : funTitleLeaf(leaf);
    this->typeNameLeaf = GetTypeNameLeaf(leaf);
    this->AnalyzeLeafBasicExtra();
  }

 protected:
  virtual inline void AnalyzeLeafBasicExtra() {
    // isTypeFloatingPoint = containsOfTString(typeNameLeaf, "loat") ||
    //                containsOfTString(typeNameLeaf, "ouble");
    if (typeNameLeaf.Contains("Bool") || typeNameLeaf.Contains("bool")) {
      histTypeFlavor = HistTypeFlavor::Boolean;
    } else if (typeNameLeaf.Contains("loat") ||
               typeNameLeaf.Contains("ouble")) {
      histTypeFlavor = HistTypeFlavor::FloatingPoint;
    } else {
      histTypeFlavor = HistTypeFlavor::Integer;
    }
    if (histTypeFlavor == HistTypeFlavor::Boolean) {
      AssignHistSetting(2, 0, 2);
      // tstrHistSetting = TString::Format("(%d,%d,%d)", nBinsCorrect,
      // (Int_t)lowerCorrect, (Int_t)upperCorrect);
    }
    if (!isHistSettingAssignedDirectly &&
        assignHistSettingPerLeafTreeExtra != nullptr) {
      isHistSettingAssignedDirectly = assignHistSettingPerLeafTreeExtra(
          nBinsCorrect, lowerCorrect, upperCorrect, nameTT, nameLeafModified,
          typeNameLeaf, titleLeaf);
    }
    if (debug && isHistSettingAssignedDirectly)
      std::cout << "Histogram setting is set directly to"
                << " nBinsCorrect: " << nBinsCorrect
                << " lowerCorrect: " << lowerCorrect
                << " upperCorrect: " << upperCorrect << " " << tstrHistSetting
                << std::endl;
  }

 public:
  TString GetNameLeafModified() const { return this->nameLeafModified; }
  TString GetTitleLeaf() const { return this->titleLeaf; }
  TString GetTypeNameLeaf() const { return this->typeNameLeaf; }
  static inline Int_t GetNBins(TH1 *hist) { return hist->GetNbinsX(); }
  static inline Double_t GetMinimumValue(TH1 *hist) {
    return hist->GetBinCenter(hist->FindFirstBinAbove()) -
           hist->GetBinWidth(1) / 2;
  }
  static inline Double_t GetMaximumValue(TH1 *hist) {
    return hist->GetBinCenter(hist->FindLastBinAbove()) +
           hist->GetBinWidth(1) / 2;
  }
  Bool_t GetHasTarget(TTree *tree) const {
    for (auto nameDeps : vDeps) {
      if (tree->GetLeaf(nameDeps) == nullptr) {
        if (debug)
          std::cout << nameDeps << " not found in tree " << tree->GetName()
                    << "(" << tree << ")"
                    << " in the current file." << std::endl;
        return false;
      }
    }
    if (funHasTargetExtra != nullptr) {
      return funHasTargetExtra(tree);
    }
    return true;
  }
  virtual void AnalyzeLeaf(TLeaf *leaf) {
    TString nameLeaf = leaf->GetName();
    vNameLeafFile.push_back(nameLeaf);
    if (isHistSettingAssignedDirectly) {
      if (dontCheckEmptyness) {
        return;
      } else if (!areAllEmpty) {
        return;
      }
    }
    if (debug)
      std::cout << "Analyzing leaf: " << nameLeaf << " (" << leaf
                << ")"
                   " tree: "
                << nameTT << " (" << leaf->GetBranch()->GetTree() << ")"
                << std::endl;
    if (leaf->GetLen() == 0) {
      if (debug) std::cout << "leaf->GetLen() returns 0" << std::endl;
      vIsEmpty.push_back(true);
      return;
    }
    expressionBeforeSetting = nameLeaf;
    EvaluateAndAnalyze(leaf->GetBranch()->GetTree());
  }
  void EvaluateAndAnalyze(TTree *tree) {
    if (debug)
      std::cout << "Analyzing expression: " << expressionBeforeSetting
                << " tree: " << nameTT << " (" << tree << ")"
                << " selection: " << selection << " option: " << optionDraw
                << " nEntriesToDrawMax: " << nEntriesToDrawMax
                << " firstEntryToDraw: " << firstEntryToDraw << std::endl;
    isEverAnalyzed = true;
    if (tree->GetEntriesFast() == 0) {
      if (debug)
        std::cout << "Tree " << tree->GetName() << " has 0 entry." << std::endl;
      vIsEmpty.push_back(true);
      vNEntriesFile.push_back(0);
      return;
    }
    TString nameHistAutogen = "h" + this->nameLeafModified + "Autogen" + "Temp";
    tree->Draw(expressionBeforeSetting + ">>" + nameHistAutogen, selection,
               optionDraw, nEntriesToDrawMax, firstEntryToDraw);
    TH1 *histAutogen = (TH1 *)gDirectory->Get(nameHistAutogen);
    AnalyzeHist(histAutogen);
  }

 protected:
  static inline TString GetTStrHistSetting(
      Int_t nBinsCorrect, Double_t lowerCorrect, Double_t upperCorrect,
      HistTypeFlavor histTypeFlavor = HistTypeFlavor::FloatingPoint) {
    if (histTypeFlavor == HistTypeFlavor::FloatingPoint) {
      return TString::Format("(%d,%F,%F)", nBinsCorrect, lowerCorrect,
                             upperCorrect);
    } else {
      return TString::Format("(%d,%lld,%lld)", nBinsCorrect,
                             (Long64_t)lowerCorrect, (Long64_t)upperCorrect);
    }
  }
  void CalculateTStrHistSetting() {
    this->tstrHistSetting = GetTStrHistSetting(nBinsCorrect, lowerCorrect,
                                               upperCorrect, histTypeFlavor);
  }
  void AnalyzeHist(TH1 *histAutogen) {
    if (debug) {
      if (histAutogen == nullptr) {
        Fatal("HistMerger::LeafAnalyzerDefault::AnalyzeHist",
              "histAutogen is nullptr!\n"
              "nameTT: %s, "
              "expressionBeforeSetting: %s, "
              "selection: %s, "
              "optionDraw: %s, "
              "nEntriesToDrawMax: %lld, "
              "firstEntryToDraw: %lld",
              nameTT.Data(), expressionBeforeSetting.Data(),
              selection.GetTitle(), optionDraw, nEntriesToDrawMax,
              firstEntryToDraw);
      }
    }
    Int_t nEntries = histAutogen->GetEntries();
    vNEntriesFile.push_back(nEntries);
    if (nEntries == 0) {
      if (debug) std::cout << "Histogram is empty." << std::endl;
      // Empty histogram
      vIsEmpty.push_back(true);
      if (isFirstNonEmptyLeaf && !isHistSettingAssignedDirectly &&
          !tstrHistSetting.Length()) {
        // nBinsCorrect = GetNBins(histAutogen);
        // lowerCorrect = histAutogen->GetBinCenter(1) - 0.5;
        // upperCorrect =
        //     histAutogen->GetBinCenter(histAutogen->GetNbinsX()) + 0.5;
        // CalculateTStrHistSetting();
      }
    } else {
      vIsEmpty.push_back(false);
      areAllEmpty = false;
      // if (histTypeFlavor == HistTypeFlavor::Boolean) {
      if (isHistSettingAssignedDirectly) {
      } else {
        Double_t minimumValue = GetMinimumValue(histAutogen);
        Double_t maximumValue = GetMaximumValue(histAutogen);
        if (isFirstNonEmptyLeaf) {
          minimumValueAllFiles = minimumValue;
          maximumValueAllFiles = maximumValue;
          isFirstNonEmptyLeaf = false;
        } else {
          if (minimumValueAllFiles > minimumValue)
            minimumValueAllFiles = minimumValue;
          if (maximumValueAllFiles < maximumValue)
            maximumValueAllFiles = maximumValue;
        }
        if (debug)
          std::cout << ", minimumValue: " << minimumValue
                    << ", maximumValue: " << maximumValue << std::endl;
        if (histTypeFlavor == HistTypeFlavor::Integer) {
        } else {
          Int_t nBins = histAutogen->GetNbinsX();
          vNBinsFile.push_back(nBins);
          vLowerFile.push_back(histAutogen->GetBinCenter(1) -
                               histAutogen->GetBinWidth(1) / 2);
          vUpperFile.push_back(histAutogen->GetBinCenter(nBins) +
                               histAutogen->GetBinWidth(1) / 2);
          if (debug)
            std::cout << "bin number: " << vNBinsFile.back() << std::endl;
        }
      }
    }
    // histAutogen->Delete();
    delete histAutogen;
  }

 public:
  std::vector<TString> GetVNameLeafFile() const { return this->vNameLeafFile; }
  std::vector<Bool_t> GetVIsEmpty() const { return this->vIsEmpty; }
  void Summarize() {
    if (!isEverAnalyzed && !allowNeverAnalyzed) {
      Fatal("HistMerger::LeafAnalyzerDefault::Summarize",
            "Expression %s has never been analyzed!\n"
            "nameTT: %s, expressionBeforeSetting: %s",
            expressionBeforeSetting.Data(), GetNameTT().Data(),
            expressionBeforeSetting.Data());
    }
    if (debug && !isHistSettingAssignedDirectly)
      std::cout << "minimumValueAllFiles, maximumValueAllFiles: "
                << minimumValueAllFiles << ", " << maximumValueAllFiles
                << std::endl;
    if (isHistSettingAssignedDirectly) {
      if (tstrHistSetting.Length() == 0) {
        CalculateTStrHistSetting();
      }
    } else {
      if (areAllEmpty) {
        if (!tstrHistSetting.Length()) {
          nBinsCorrect = 0;
          lowerCorrect = 0.;
          upperCorrect = 1.;
        }
      } else {
        tstrHistSetting = "";
        if (histTypeFlavor == HistTypeFlavor::Integer) {
          minimumValueAllFiles =
              (TMath::Abs(minimumValueAllFiles) <= __INT32_MAX__)
                  ? TMath::Nint(minimumValueAllFiles)
                  : std::round(minimumValueAllFiles);
          maximumValueAllFiles =
              (TMath::Abs(maximumValueAllFiles) <= __INT32_MAX__)
                  ? TMath::Nint(maximumValueAllFiles)
                  : std::round(maximumValueAllFiles);
          if (minimumValueAllFiles < 0 && maximumValueAllFiles > 0) {
            lowerCorrect = minimumValueAllFiles;
            upperCorrect = maximumValueAllFiles + 1;
          } else if (minimumValueAllFiles > 0 &&
                     (minimumValueAllFiles < 20 ||
                      (minimumValueAllFiles <
                           (((Long64_t)maximumValueAllFiles) >> 1) &&
                       minimumValueAllFiles <
                           intpow(10, getNDigitsM1(maximumValueAllFiles),
                                  (Double_t)1.)))) {
            lowerCorrect = 0;
            upperCorrect = maximumValueAllFiles + 1;
          } else if (maximumValueAllFiles < 0 &&
                     (maximumValueAllFiles > -20 ||
                      (maximumValueAllFiles >
                           -(((Long64_t)minimumValueAllFiles) >> 1) &&
                       maximumValueAllFiles >
                           intpow(10, getNDigitsM1(minimumValueAllFiles),
                                  (Double_t)-1.)))) {
            lowerCorrect = minimumValueAllFiles;
            upperCorrect = 0;
          } else {
            lowerCorrect = minimumValueAllFiles;
            upperCorrect = maximumValueAllFiles + 1;
          }
          if (upperCorrect - lowerCorrect > __INT32_MAX__) {
            ULong64_t l64NBinsCorrect = upperCorrect - lowerCorrect;
            UShort_t nbitsReduced = 0;
            while (l64NBinsCorrect & ~((ULong64_t)__INT32_MAX__)) {
              nbitsReduced++;
              l64NBinsCorrect >>= 1;
            }
            nBinsCorrect = l64NBinsCorrect;
          } else {
            nBinsCorrect = upperCorrect - lowerCorrect;
          }
          // tstrHistSetting =
          //     TString::Format("(%d,%lld,%lld)", nBinsCorrect,
          //                     (Long64_t)lowerCorrect,
          //                     (Long64_t)upperCorrect);
        } else {
          UInt_t nHists = vNBinsFile.size();
          std::sort(vNBinsFile.begin(), vNBinsFile.end());
          nBinsCorrect = vNBinsFile[nHists >> 1];
          if (minimumValueAllFiles <= 0 && maximumValueAllFiles >= 0) {
            if (minimumValueAllFiles != 0) {
              Int_t nDigitsM1Minimum = getNDigitsM1(minimumValueAllFiles);
              lowerCorrect = intpow(10, nDigitsM1Minimum,
                                    TMath::Floor(intpow(10, -nDigitsM1Minimum,
                                                        minimumValueAllFiles)));
            } else {
              lowerCorrect = 0.;
            }
            if (maximumValueAllFiles != 0) {
              Int_t nDigitsM1Maximum = getNDigitsM1(maximumValueAllFiles);
              upperCorrect = intpow(10, nDigitsM1Maximum,
                                    TMath::Ceil(intpow(10, -nDigitsM1Maximum,
                                                       maximumValueAllFiles)));
            } else {
              upperCorrect = 0.;
            }
            if (lowerCorrect == 0 && upperCorrect == 0) {
              lowerCorrect -= static_cast<Double_t>(__FLT_EPSILON__);
              upperCorrect += static_cast<Double_t>(__FLT_EPSILON__);
            }
            // } else if (minimumValueAllFiles == 0) {
            //   lowerCorrect = minimumValueAllFiles;
            //   upperCorrect = maximumValueAllFiles;
            //   if (lowerCorrect == upperCorrect) {
            //     upperCorrect += __FLT_EPSILON__;
            //   }
            // } else if (maximumValueAllFiles == 0) {
            //   lowerCorrect = minimumValueAllFiles;
            //   upperCorrect = maximumValueAllFiles;
            //   if (lowerCorrect == upperCorrect) {
            //     lowerCorrect -= __FLT_EPSILON__;
            //   }
          } else if (minimumValueAllFiles == maximumValueAllFiles) {
            lowerCorrect = minimumValueAllFiles - TMath::Min(TMath::Abs(minimumValueAllFiles),
                                       static_cast<Double_t>(__FLT_EPSILON__));
            upperCorrect += maximumValueAllFiles + TMath::Min(TMath::Abs(maximumValueAllFiles),
                                       static_cast<Double_t>(__FLT_EPSILON__));
          } else {
            Bool_t isPositive = maximumValueAllFiles > 0;
            Double_t valueOuter =
                isPositive ? maximumValueAllFiles : -minimumValueAllFiles;
            Double_t valueInner =
                isPositive ? minimumValueAllFiles : -maximumValueAllFiles;
            Int_t nDigitsM1Outer = getNDigitsM1(valueOuter);
            Int_t orderDiff = TMath::Nint(
                TMath::Log10(maximumValueAllFiles - minimumValueAllFiles));
            Int_t mDigitsToRound = TMath::Min(nDigitsM1Outer, orderDiff);
            // Previous merging error is caused by
            // mistakenly exchanging of TMath::Floor and TMath::Ceilling
            // causing valueOuter to be equal to valueInner when mDigitsToRound
            // == nDigitsM1Outer e.g. When minimumAllFile is 20 and
            // maximumAllFile is 1016 lowerCorrect and upperCorrect were
            // mistakenly assigned to be 1000 and 1000
            valueOuter =
                intpow(10, mDigitsToRound,
                       TMath::Ceil(intpow(10, -mDigitsToRound, valueOuter)));
            valueInner =
                intpow(10, mDigitsToRound,
                       TMath::Floor(intpow(10, -mDigitsToRound, valueInner)));
            if (valueInner == valueOuter) {
              valueOuter += TMath::Min(valueOuter,
                                       static_cast<Double_t>(__FLT_EPSILON__));
              valueInner -= TMath::Min(valueInner,
                                       static_cast<Double_t>(__FLT_EPSILON__));
            }
            if (isPositive) {
              lowerCorrect = valueInner;
              upperCorrect = valueOuter;
            } else {
              lowerCorrect = -valueOuter;
              upperCorrect = -valueInner;
            }
          }
          // tstrHistSetting = TString::Format("(%d,%F,%F)", nBinsCorrect,
          //                                   lowerCorrect, upperCorrect);
        }
        // if (tstrHistSetting.Length() == 0) {
        //   tstrHistSetting = TString::Format("(%d,%F,%F)", nBinsCorrect,
        //                                     lowerCorrect, upperCorrect);
        // }
      }
      if (adjustHistSettingPerLeafTreeExtra != nullptr) {
        adjustHistSettingPerLeafTreeExtra(
            nBinsCorrect, lowerCorrect, upperCorrect, nameTT, nameLeafModified,
            typeNameLeaf, titleLeaf);
      }
      CalculateTStrHistSetting();
    }
    if (debug && tstrHistSetting.Length() == 0) {
      Fatal("HistMerger::LeafAnalyzerDefault::Summarize",
            "tstrHistSetting has not been set when function Summarize reaches "
            "the end");
    }
    if (debug)
      std::cout << TString::Format(
                       "nBinsCorrect: %d, lowerCorrect: %F, upperCorrect: %f, "
                       "tstrHistSetting: %s",
                       nBinsCorrect, lowerCorrect, upperCorrect,
                       tstrHistSetting.Data())
                << std::endl;
  }
  Bool_t GetAreAllEmpty() const { return areAllEmpty; }
  Int_t GetNBinsCorrect() const { return nBinsCorrect; }
  Double_t GetLowerCorrect() const { return lowerCorrect; }
  Double_t GetUpperCorrect() const { return upperCorrect; }
  TString GetHistSetting() const { return tstrHistSetting; }
  TH1 *DrawHistCorrected(TString nameHist, TTree *tree, Int_t iHist,
                         Bool_t isToClone) {
    if (debug && !GetHistSetting().Length()) {
      Fatal("HistMerger::LeafAnalyzerDefault::DrawHistCorrected", "GetHistSetting() returns \"\"");
    }
    if (debug && iHist >= 0 && GetVNameLeafFile().size() < static_cast<size_t>(iHist) + 1) {
      Fatal("HistMerger::LeafAnalyzerDefault::DrawHistCorrected",
          "iHist (%d) exceeds the boundary of vNameLeafFile (size: %ld)",
          iHist, vNameLeafFile.size());
    }
    TString expressionFull =
        (iHist < 0 ? expressionBeforeSetting : GetVNameLeafFile()[iHist]) +
        ">>" + nameHist + GetHistSetting();
    TH1 *histNotCloned = nullptr;
    if (debug)
      std::cout << "On tree " << tree << " " << tree->GetName() << std::endl;
    if ((dontCheckEmptyness || iHist < 0 || vIsEmpty[iHist]) &&
        tree->GetEntriesFast() == 0) {
      if (debug)
        std::cout << "Tree " << tree->GetName()
                  << " has 0 entries.\nReturning a new empty histogram."
                  << std::endl;

      return new TH1F(nameHist, titleLeaf, nBinsCorrect, lowerCorrect,
                      upperCorrect);
    } else {
      if (debug && tree->GetEntriesFast() == 0) {
        Fatal("HistMerger::LeafAnalyzerDefault",
              "Tree has 0 etry but histogram isn't marked empty.");
      }
      if (debug)
        std::cout << "Drawing " << expressionFull << " ..." << std::flush;
      tree->Draw(expressionFull, selection, optionDraw, nEntriesToDrawMax,
                 firstEntryToDraw);
      if (debug)
        std::cout << " done."
                  << "\n";
      if (debug) std::cout << "Getting ...";
      histNotCloned = (TH1 *)gDirectory->Get(nameHist);
    }
    if (histNotCloned == nullptr) {
      Fatal("HistMerger::LeafAnalyzerDefault::DrawHistCorrected",
            "tree->Draw() returns nullptr!\n"
            "tree: %s (%p), "
            "expressionFull: %s, "
            "selection: %s, "
            "optionDraw: %s"
            "nEntriesToDrawMax: %lld"
            "firstEntryToDraw: %lld",
            nameTT.Data(), tree, expressionFull.Data(), selection.GetTitle(),
            optionDraw, nEntriesToDrawMax, firstEntryToDraw);
    }
    if (debug) std::cout << " done." << std::endl;
    if (isToClone) {
      TH1 *histCloned = (TH1 *)histNotCloned->Clone(nameHist);
      histNotCloned->Delete();
      delete histNotCloned;
      return histCloned;
    }
    return histNotCloned;
  }
};

void HistMerger::InitializeHidden() {
  this->funNameLeafModified = nullptr;
  this->supplyLeafAnalyzer = nullptr;
  this->nTT = 0;
  this->vWeightDataset.clear();
  this->vWeightFile.clear();
  this->vvvIsHistFileLeafTree.clear();
  this->vvvIsHistFileLeafTreeCustom.clear();
  this->vvHistResultLeafTree.clear();
  this->vTFOut.clear();
  this->vNLeavesTree.clear();
  this->vNAnalyzersTreeCustom.clear();
  this->vvNameModifiedLeafTree.clear();
  this->vvAnalyzerLeafTree.clear();
}

inline void HistMerger::ProcessAnalyzerNew(LeafAnalyzerAbstract *analyzer,
                                           TString nameTT) {
  analyzer->SetDebug(debug);
  analyzer->SetNameTT(nameTT);
  analyzer->SetFunAssignHistSettingExtra(assignHistSettingPerLeafTreeExtra);
  analyzer->SetFunAdjustHistSettingExtra(adjustHistSettingPerLeafTreeExtra);
  analyzer->SetFunTitleLeaf(funTitleLeaf);
}

void HistMerger::InitializeWhenRun() {
  if (funPathTFIn == nullptr) {
    Fatal("HistMerger::InitializeWhenRun", "funPathTFIn is not specified.");
  }
  if (dirTFTemp == "" && nLeavesToUseCorrectedTempFileMin) {
    Warning("HistMerger::InitializeWhenRun",
            "dirTFTemp is empty. Use \"temp\"");
    dirTFTemp = "temp";
  }
  if (funNameTFTemp == nullptr && nLeavesToUseCorrectedTempFileMin) {
    Warning("HistMerger::InitializeWhenRun",
            "funNameTFTemp is not specified. \n"
            "Assign to default value from mkstemp.");
    char *charsTemp = nullptr;
    strcat(charsTemp, "XXXXXX");
    mkstemp(charsTemp);
    SetNameTFTemp((TString) "%s" + charsTemp + ".root");
  }
  if (dirTFOut == "") {
    Warning("HistMerger::InitializeWhenRun", "dirTFOut is empty. Use \".\"");
    dirTFOut = ".";
  }
  if (funNameTFTemp == nullptr) {
    Warning("HistMerger::InitializeWhenRun",
            "funNameTFOut is not specified.\n"
            "Set with format \"output_%%s.root\"");
    SetNameTFOut("output_%s.root");
  }
  if (supplyLeafAnalyzer == nullptr || funNameLeafModified == nullptr) {
    SetLeafAnalyzer<LeafAnalyzerDefault>();
  }
  // if (supplyLeafAnalyzer == nullptr) {
  //   supplyLeafAnalyzer = []() -> LeafAnalyzerAbstract * {
  //     return (LeafAnalyzerAbstract *)(new LeafAnalyzerDefault());
  //   };
  // }
  // if (funNameLeafModified == nullptr) {
  //   funNameLeafModified = [](TString nameLeaf) -> TString {
  //     return LeafAnalyzerDefault::GetNameLeafModified(nameLeaf);
  //   };
  //   // funNameLeafModified =
  //   static_cast<TString(*)(TString)>(LeafAnalyzerDefault::GetNameLeafModified);
  // }
  nTT = vNameTT.size();
  if (nTT == 0) {
    Warning("HistMerger::InitializeWhenRun",
            "No tree name specified. Nothing will be done.");
  }
  if (vvAnalyzerLeafTreeCustom.size() < nTT) {
    Warning(
        "HistMerger::InitializeWhenRun",
        "Size of vvAnalizerLeafTreeCustom is smaller than the number of trees\n"
        "Fill the remaining space with empty vectors");
    vvAnalyzerLeafTreeCustom.reserve(nTT);
    for (UInt_t i = vvAnalyzerLeafTreeCustom.size(); i < nTT; i++) {
      std::vector<LeafAnalyzerAbstract *> vAnalyzerLeafCustom;
      vAnalyzerLeafCustom.clear();
      vvAnalyzerLeafTreeCustom.push_back(vAnalyzerLeafCustom);
    }
  } else if (vvAnalyzerLeafTreeCustom.size() > nTT) {
    Warning(
        "HistMerger::InitializeWhenRun",
        "Size of vvAnalizerLeafTreeCustom is larger than the number of trees\n"
        "Truncate to fit the correct number");
    vvAnalyzerLeafTreeCustom.resize(nTT);
  }
  if (debug) {
    std::cout << "Initial custom analyzer ";
    for (UInt_t indexNameCustom=0; indexNameCustom<nTT; indexNameCustom++) {
      std::cout << vvAnalyzerLeafTreeCustom[indexNameCustom].size();
      if (indexNameCustom != nTT-1) {
        std::cout << "+";
      }
    }
    std::cout << std::endl;
  }
  for (LeafAnalyzerAbstract *analyzer : vAnalyzerCustomByName) {
    if (analyzer->GetNameTT().Length() == 0) {
      for (std::vector<LeafAnalyzerAbstract *> &vAnalyzerLeafCustom :
           vvAnalyzerLeafTreeCustom) {
        vAnalyzerLeafCustom.push_back(analyzer);
      }
      continue;
    }
    typename std::vector<TString>::iterator iterNameTT = vNameTT.begin();
    while (true) {
      iterNameTT = std::find(iterNameTT, vNameTT.end(), analyzer->GetNameTT());
      if (iterNameTT == vNameTT.end()) {
        break;
      }
      UInt_t iTree = std::distance(vNameTT.begin(), iterNameTT);
      vvAnalyzerLeafTreeCustom[iTree].push_back(analyzer);
      iterNameTT++;
    }
  }
  for (UInt_t iTree = 0; iTree < vNameTT.size(); iTree++) {
    for (LeafAnalyzerAbstract *analyzer : vvAnalyzerLeafTreeCustom[iTree]) {
      ProcessAnalyzerNew(analyzer, vNameTT[iTree]);
    }
  }
  for (UInt_t iTree = 0; iTree < nTT; iTree++) {
    std::vector<std::vector<Bool_t>> vvIsHistFileLeaf;
    vvIsHistFileLeaf.clear();
    vvvIsHistFileLeafTree.push_back(vvIsHistFileLeaf);
    std::vector<std::vector<Bool_t>> vvIsHistFileLeafCustom(
        vvAnalyzerLeafTreeCustom[iTree].size(), std::vector<Bool_t>());
    vvvIsHistFileLeafTreeCustom.push_back(vvIsHistFileLeafCustom);
    // std::vector<std::vector<TH1 *>> vvHistFileLeaf;
    // vvHistFileLeaf.clear();
    std::vector<LeafAnalyzerAbstract *> vAnalyzerLeaf;
    vAnalyzerLeaf.clear();
    vvAnalyzerLeafTree.push_back(vAnalyzerLeaf);
    std::vector<TString> vNameModifiedLeaf;
    vNameModifiedLeaf.clear();
    vvNameModifiedLeafTree.push_back(vNameModifiedLeaf);
    // std::vector<TString> vTitleLeaf;
    // vTitleLeaf.clear();
    // vvTitleLeafTree.push_back(vTitleLeaf);
    // std::vector<TString> vTypeNameLeaf;
    // vTypeNameLeaf.clear();
    // vvTypeNameLeafTree.push_back(vTypeNameLeaf);
    // TList *tltlHistFileLeaf = new TList;
    // (*toatltlHistFileLeafTree)[iTree] = tltlHistFileLeaf;
  }
}

void HistMerger::Run() {
  InitializeWhenRun();

  Bool_t areSomeInfilesClosed = false;
  UInt_t nInfileOpen = 0;
  // std::function<TString(TString)> modifyNameLeaf =
  //     [](TString nameLeaf) -> TString {
  //   TString nameLeafNew = TString(nameLeaf);
  //   Ssiz_t indexLeftSquareBracket = nameLeaf.First('[');
  //   if (indexLeftSquareBracket != -1) {
  //     nameLeafNew.Resize(indexLeftSquareBracket);
  //   }
  //   return nameLeafNew;
  // };

  Long64_t nEntryOriginal;

  TTree *arrTT[nTT];

  auto mountVarsFromTDir = [this, &nEntryOriginal, &arrTT](TDirectory *tdir) {
    nEntryOriginal = (*((TVectorD *)tdir->Get("tvdNEntryOriginal")))[0];
    for (UInt_t i = 0; i < nTT; i++) {
      arrTT[i] = (TTree *)tdir->Get(this->vNameTT[i]);
    }
  };

  UInt_t nDataset = vNumberFile.size();
  if (nDataset != vCrossSection.size()) {
    std::cerr << "The length of fileNumbers (" << nDataset
              << ") does not match that of the crossSection ("
              << vCrossSection.size() << ")" << std::endl;
  }
  nDataset = TMath::Min(nDataset, (UInt_t)vCrossSection.size());

  UInt_t nFileTot = 0;
  for (UInt_t iDataset = 0; iDataset < nDataset; iDataset++) {
    nFileTot += vNumberFile[iDataset];
  }

  if (debug) std::cout << "nFileTot: " << nFileTot << std::endl;
  TCanvas *c1;

  const UInt_t nFileTotOriginal = nFileTot;
  const std::vector<UInt_t> vNumberFileOriginal = vNumberFile;
  TFile *arrFile[nFileTotOriginal];
  Bool_t arrIsOpenableFile[nFileTotOriginal];  //< Array of whether each file is
                                               // openable or not (The per-file
                                               // vectors cares only about
                                               // openable files).

  auto closeTFilesIn = [this, &nInfileOpen, &arrIsOpenableFile, &arrFile](
                           Int_t iFileOriginalToClose,
                           std::function<void()> funDoBeforeClose = nullptr) {
    // tfAutogenHist->Write();
    // areSomeInfilesClosed = true;
    if (funDoBeforeClose != nullptr) funDoBeforeClose();
    while (nInfileOpen > 0 && iFileOriginalToClose >= 0) {
      if (arrIsOpenableFile[iFileOriginalToClose]) {
        arrFile[iFileOriginalToClose]->Close();
        delete arrFile[iFileOriginalToClose];
        nInfileOpen--;
      }
      iFileOriginalToClose--;
    }
    if (this->debug && nInfileOpen > 0) {
      std::cerr << "Closing wrong number of files" << std::endl;
    }
  };

  UInt_t nLeavesAllTrees = 0;
  auto getNameHistOfFileFromAnalyzerAndIFile =
      [](LeafAnalyzerAbstract *analyzer, UInt_t iFile, TString suffixStage) {
        return "h" + analyzer->GetNameTT() + analyzer->GetNameLeafModified() +
               suffixStage + iFile;
      };
  // auto getNameHistOfFile = [this](UInt_t iTree, UInt_t indexName, UInt_t
  // iFile,
  //                                 TString suffixStage) -> TString {
  //   return (TString) "h" + vNameTT[iTree] +
  //          vvNameModifiedLeafTree[iTree][indexName] + suffixStage + iFile;
  // };

  if (dirTFTemp.Length() > 1 && dirTFTemp.EndsWith(seperatorPath)) {
    dirTFTemp.Resize(dirTFTemp.Length() - 1);
  }
  if (!gSystem->AccessPathName(dirTFTemp, kFileExists)) {
    gSystem->mkdir(dirTFTemp);
  }
  TString pathTFAutogenHist =
      (dirTFTemp == "" ? "" : dirTFTemp + seperatorPath) +
      funNameTFTemp("Autogen");
  TString pathTFCorrectedHist =
      (dirTFTemp == "" ? "" : dirTFTemp + seperatorPath) +
      funNameTFTemp("Corrected");
  const std::function<TString(TString nameTT, TString nameLeafModified)>
      getNameHistResultFromNames =
          [](TString nameTT, TString nameLeafModified) {
            return "h" + nameLeafModified;
          };
  // const std::function<TString(UInt_t iTree, UInt_t indexName)>
  //     getNameHistResult = [this, getNameHistResultFromNames](
  //                             UInt_t iTree, UInt_t indexName) -> TString {
  //   return getNameHistResultFromNames(vNameTT[iTree],
  //                                     vvNameModifiedLeafTree[iTree][indexName]);
  // };
  if (dirTFOut.Length() > 1 && dirTFOut.EndsWith(seperatorPath)) {
    dirTFOut.Resize(dirTFOut.Length() - 1);
  }
  if (!gSystem->AccessPathName(dirTFOut, kFileExists)) {
    gSystem->mkdir(dirTFOut);
  }
  const std::function<TString(UInt_t)> funPathTFOutHist =
      [this](UInt_t iTree) -> TString {
    return dirTFOut + seperatorPath + funNameTFOut(vNameTT[iTree]);
  };
  // const std::function<void(Int_t &, Double_t &, Double_t &, UInt_t, UInt_t)>
  //     adjustHistSettingPerLeafTreeExtraIndex = [this](Int_t &nBinCorrect,
  //                                                     Double_t &lowerCorrect,
  //                                                     Double_t &upperCorrect,
  //                                                     UInt_t iTree,
  //                                                     UInt_t indexName) {
  //       TString nameTT = vNameTT[iTree];
  //       TString &nameLeafModified = vvNameModifiedLeafTree[iTree][indexName];
  //       // TString typeName = vvTypeNameLeafTree[iTree][indexName];
  //       // TString title = vvTitleLeafTree[iTree][indexName];
  //       LeafAnalyzerAbstract *&analyzer =
  //       vvAnalyzerLeafTree[iTree][indexName]; TString typeName =
  //       analyzer->GetTypeNameLeaf(); TString title =
  //       analyzer->GetTitleLeaf(); return adjustHistSettingPerLeafTreeExtra(
  //           nBinCorrect, lowerCorrect, upperCorrect, nameTT,
  //           nameLeafModified, typeName, title);
  //     };
  // UInt_t iFile = 0;
  // TFile *tfAutogenHist, *tfCorrectedHist;
  TFile *tfCorrectedHist;
  // tfAutogenHist =
  //     TFile::Open(pathTFAutogenHist, toRecreateOutFile ? "recreate" :
  //     "update");
  // TFile *tfCorrectedHist = TFile::Open(pathTFCorrectedHist, toRecreateOutFile
  // ? "recreate" : "update");
  for (UInt_t iDataset = 0, iFile = 0, nFileMissingTot = 0; iDataset < nDataset;
       iDataset++) {
    Long64_t nEntryOriginalDatasetCurrent = 0;
    if (debug) std::cout << "iDataset: " << iDataset << std::endl;
    Bool_t isFileOpenable = true;
    for (UInt_t jFile = 0, nFileMissingThisDataset = 0;
         jFile < vNumberFile[iDataset];
         (isFileOpenable ? jFile++ : nFileMissingThisDataset++),
                (isFileOpenable ? iFile++ : nFileMissingTot++)) {
      // TString filenameCurrent = filenameGeneralCurrentCluster + "_" + (iFile
      // + nFileMissingTot) + ".root";
      TString pathFileIn = funPathTFIn(iFile + nFileMissingTot);
      if (debug)
        std::cout << "Opening TFile " << pathFileIn << " ..." << std::endl;
      gSystem->ExpandPathName(pathFileIn);
      TFile *tfCurrent = nullptr;
      // gSystem->AccessPathName(pathFileIn,kFileExists) returns FALSE if the
      // file of pathFileIn EXISTS!
      isFileOpenable = !gSystem->AccessPathName(pathFileIn, kFileExists);
      if (isFileOpenable) {
        tfCurrent = TFile::Open(pathFileIn, "READ");
        isFileOpenable = tfCurrent != nullptr && !tfCurrent->IsZombie();
      } else {
        std::cerr << "File not exists!" << std::endl;
      }
      arrIsOpenableFile[iFile + nFileMissingTot] = isFileOpenable;
      arrFile[iFile + nFileMissingTot] = tfCurrent;
      if (!isFileOpenable) {
        std::cerr << "File not openable: " << pathFileIn << " (" << tfCurrent
                  << ")" << std::endl;
        if (allowMissingInputFiles) {
          vNumberFile[iDataset]--;
          nFileTot--;
          continue;
        }
      }
      nInfileOpen++;
      if (debug) std::cout << "Done." << std::endl;
      if (debug)
        std::cout << "Mounting variables from the file ..." << std::endl;
      mountVarsFromTDir(tfCurrent);
      if (debug) std::cout << "Done." << std::endl;
      if (debug) std::cout << "nEntryOriginal: " << nEntryOriginal << std::endl;
      nEntryOriginalDatasetCurrent += nEntryOriginal;
      for (UInt_t iTree = 0; iTree < nTT; iTree++) {
        if (debug) std::cout << "iTree: " << iTree << std::endl;
        if (arrTT[iTree] == nullptr) {
          Error("HistMerger::Run",
                "Tree %s (iTree=%d) does not exist in file iFile=%d "
                "iFileOriginal=%d",
                vNameTT[iTree].Data(), iTree, iFile, iFile + nFileMissingTot);
          for (auto &vvIsHistFile : vvvIsHistFileLeafTree[iTree]) {
            vvIsHistFile.push_back(false);
          }
          for (auto &vvIsHistFileCustom : vvvIsHistFileLeafTreeCustom[iTree]) {
            vvIsHistFileCustom.push_back(false);
          }
          continue;
        }
        if (debug)
          std::cout << "Getting leafnames from the tree ..." << std::endl;
        // TList *tltlHistFileLeaf = (TList *)
        // (*toatltlHistFileLeafTree)[iTree]; std::vector<TString>
        // vNameModifiedLeafOld = vvNameModifiedLeafTree[iTree]; const UInt_t
        // nNameModifiedLeafOld = vNameModifiedLeafOld.size();
        const UInt_t nNameModifiedLeafOld =
            vvNameModifiedLeafTree[iTree].size();
        Bool_t arrIsHistLeafOld[nNameModifiedLeafOld];
        for (UInt_t indexNameOld = 0; indexNameOld < nNameModifiedLeafOld;
             indexNameOld++) {
          arrIsHistLeafOld[indexNameOld] = false;
        }  //< For possible missing file and changing nFileTot
        for (TObject *leafRaw : *(arrTT[iTree]->GetListOfLeaves())) {
          TLeaf *leaf = (TLeaf *)leafRaw;
          TString nameLeaf = leaf->GetName();
          TString nameLeafModified = funNameLeafModified(nameLeaf);
          if (funIsToVetoLeaf != nullptr &&
              funIsToVetoLeaf(vNameTT[iTree], nameLeafModified))
            continue;
          LeafAnalyzerAbstract *analyzer = nullptr;
          Bool_t isIncluded = false;
          UInt_t indexName = 0;
          // for (TString nameIncluded : vvNameModifiedLeafTree[iTree]) {
          //   // if (debug) std::cout << "indexName : " << indexName;
          //   if (nameLeafModified == nameIncluded) {
          //     isIncluded = true;
          //     // if (debug) std::cout << " is included.";
          //     break;
          //   }
          //   indexName++;
          // }
          indexName = std::distance(
              vvNameModifiedLeafTree[iTree].begin(),
              std::find(vvNameModifiedLeafTree[iTree].begin(),
                        vvNameModifiedLeafTree[iTree].end(), nameLeafModified));
          isIncluded = indexName != vvNameModifiedLeafTree[iTree].size();
          // if (debug) std::cout << std::endl;
          if (!isIncluded) {
            if (debug) std::cout << "Found new leafname." << std::endl;
            if (debug)
              std::cout << "nameLeaf: " << nameLeaf
                        << ", nameLeafModified: " << nameLeafModified;
            vvNameModifiedLeafTree[iTree].push_back(nameLeafModified);
            analyzer = supplyLeafAnalyzer();
            if (debug) std::cout << ", analyzer: " << analyzer << std::endl;
            vvAnalyzerLeafTree[iTree].push_back(analyzer);
            ProcessAnalyzerNew(analyzer, vNameTT[iTree]);
            analyzer->SetAllowNeverAnalyzed(false);
            analyzer->AnalyzeLeafBasic(leaf);
            // TString titleLeaf = leaf->GetTitle();
            // if (titleLeaf == "") {
            //   titleLeaf = leaf->GetBranch()->GetTitle();
            // // }
            // vvTitleLeafTree[iTree].push_back(titleLeaf);
            // vvTypeNameLeafTree[iTree].push_back(leaf->GetTypeName());
            std::vector<Bool_t> vIsHistFile(nFileTot);
            vIsHistFile.clear();
            std::vector<TH1 *> vHistFile(nFileTot);
            vHistFile.clear();
            // Commented out for possible missing file and changing nFileTot
            // for (UInt_t i=0; i<nFileTot; i++) {
            //   vIsHistFile.push_back(false);
            // }
            // Added for possible missing file and changing nFileTot
            for (UInt_t i = 0; i < iFile; i++) {
              vIsHistFile.push_back(false);
            }
            vIsHistFile.push_back(true);
            // Now vIsHistFile has iFile+1 elements
            vvvIsHistFileLeafTree[iTree].push_back(vIsHistFile);
            // vvvHistFileLeafTree[iTree].push_back(vHistFile);
            // TList *tlHistFileCurrent = new TList;
            // tltlHistFileLeaf->AddLast(tlHistFileCurrent);
          } else {
            // indexName should be smaller than current
            // vNameModifiedLeafTree[iTree].size() unless the nameLeafModified
            // is duplicated
            if (indexName >= vvNameModifiedLeafTree[iTree].size() ||
                arrIsHistLeafOld[indexName]) {
              Error("HistMerger::Run",
                    "nameLeafModified duplication (%s, %s) (indexName: %d) at "
                    "tree %s (iTree: %d), iFile: %d",
                    nameLeaf.Data(), nameLeafModified.Data(), indexName,
                    vNameTT[iTree].Data(), iTree, iFile);
            } else {
              arrIsHistLeafOld[indexName] = true;
              analyzer = vvAnalyzerLeafTree[iTree][indexName];
            }
          }
          // vvvIsHistFileLeafTree[iTree][indexName][iFile] = true;

          if (analyzer == nullptr) {
            continue;
          }

          analyzer->AnalyzeLeaf(leaf);

          // TString nameHistAutogen =
          //     getNameHistOfFile(iTree, indexName, iFile, "Autogen");
          // if (debug)
          //   std::cout << "Drawing " << nameLeaf << " as " << nameHistAutogen
          //             << " ...";
          // arrTT[iTree]->Draw(nameLeaf + ">>" + nameHistAutogen);
          // TH1 *histAutogen = (TH1 *)gDirectory->Get(nameHistAutogen);
          // histAutogen->SetTitle(nameLeaf);
          // histAutogen->SetName(nameHistAutogen);
          // histAutogen->SetDirectory(tfAutogenHist);
          // // histAutogen->Write(nameHistAutogen);
          // if (debug && histAutogen == nullptr)
          //   std::cerr << "Fatal: " << nameHistAutogen << " is nullptr!"
          //             << std::endl;
          // if (debug) std::cout << " " << histAutogen;
          // if (debug) std::cout << " Done." << std::endl;
        }
        // As for now, vvvIsHistFileLeafTree[iTree][indexName]
        // has iFile elements if indexName < nNameModifiedLeafOld
        // and iFile + 1 elements if indexName >= nNameModifiedLeafOld
        // Now we are going to push_back one element to every
        // vvvIsHistFileLeafTree[iTree][indexName]
        // that indexName < nNameModifiedLeafOld
        // according to arrIsHistLeafOld
        for (UInt_t indexNameOld = 0; indexNameOld < nNameModifiedLeafOld;
             indexNameOld++) {
          vvvIsHistFileLeafTree[iTree][indexNameOld].push_back(
              arrIsHistLeafOld[indexNameOld]);
        }
        if (debug) {
          for (UInt_t indexName = 0;
               indexName < vvvIsHistFileLeafTree[iTree].size(); indexName++) {
            if (vvvIsHistFileLeafTree[iTree][indexName].size() != iFile + 1) {
              std::cerr << "vvvIsHistFileLeafTree[" << iTree << "]["
                        << indexName << "].size() != iFile + 1"
                        << "("
                        << "size: "
                        << vvvIsHistFileLeafTree[iTree][indexName].size()
                        << ", iFile: " << iFile << ")" << std::endl;
            }
          }
        }
        // UInt_t nAnalyzerLeafCustomOld =
        // vvAnalyzerLeafTreeCustom[iTree].size(); for (UInt_t
        // indexNameCustomOld = 0;
        //      indexNameCustomOld < nAnalyzerLeafCustomOld;
        //      indexNameCustomOld++) {
        //   vvvIsHistFileLeafTreeCustom[iTree][indexNameCustomOld].push_back(
        //       vvAnalyzerLeafTreeCustom[iTree][indexNameCustomOld]->GetHasTarget(
        //           arrTT[iTree]));
        //   if (debug) std::cout << "Custom analyzer "
        //   <<
        //   vvAnalyzerLeafTreeCustom[iTree][indexNameCustomOld]->GetNameLeafModified()
        //   << (vvvIsHistFileLeafTreeCustom[iTree][indexNameCustomOld].back() ?
        //   "can" : "cannot")
        //   << " be applied to this file"
        //   << std::endl;
        // }
        if (debug && vvAnalyzerLeafTreeCustom[iTree].size() !=
                         vvvIsHistFileLeafTreeCustom[iTree].size()) {
          Fatal("HistMerger::Run",
                "Size of vvAnalyzerLeafTreeCustom[%d]: %ld and "
                "vvvIsHistFileLeafTreeCustom[%d]: %ld do not match!",
                iTree, vvAnalyzerLeafTreeCustom[iTree].size(), iTree,
                vvvIsHistFileLeafTreeCustom[iTree].size());
        }
        if (pushCustomAnalyzersWhenRun != nullptr) {
          if (debug)
            std::cout << "Pushing custom analyzers for tree " << vNameTT[iTree]
                      << " in file (iFile: " << iFile << " ..." << std::endl;
          UInt_t nAnalyzerLeafCustomOld =
              vvAnalyzerLeafTreeCustom[iTree].size();
          // vvAnalyzerLeafTreeCustom[iTree].push_back(nullptr);

          pushCustomAnalyzersWhenRun(
              arrTT[iTree], vvAnalyzerLeafTreeCustom[iTree],
              vvAnalyzerLeafTree[iTree],
              nNameModifiedLeafOld,
              [this, iTree](LeafAnalyzerAbstract *analyzer) {
                vvAnalyzerLeafTreeCustom[iTree].push_back(analyzer);
              });
          if (vvAnalyzerLeafTreeCustom[iTree].size() > nAnalyzerLeafCustomOld) {
            if (debug) {
              std::cout << "nAnalyzerLeafCustomOld: " << nAnalyzerLeafCustomOld << std::endl;
              for (UInt_t indexNameCustom = nAnalyzerLeafCustomOld;
                   indexNameCustom < vvAnalyzerLeafTreeCustom[iTree].size();
                   indexNameCustom++) {
                std::cout << "Push "
                          << vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]
                                 ->GetNameLeafModified()
                          << "("
                          << vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]
                          << ")"
                          << " indexNameCustom: " << indexNameCustom
                          << std::endl;
              }
            }
            std::vector<Bool_t> vIsHistFileCustomBlank(iFile, false);
            vvvIsHistFileLeafTreeCustom[iTree].resize(
                vvAnalyzerLeafTreeCustom[iTree].size(), vIsHistFileCustomBlank);
          }
          if (debug) std::cout << "Done." << std::endl;
        }
        if (debug) std::cout << "Analyzing custom leaves ..." << std::endl;
        for (UInt_t indexNameCustom = 0;
             indexNameCustom < vvAnalyzerLeafTreeCustom[iTree].size();
             indexNameCustom++) {
          if (debug) std::cout << "On analyzer " << vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]
          << " " << vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]->GetNameLeafModified() << std::endl;
          vvvIsHistFileLeafTreeCustom[iTree][indexNameCustom].push_back(
              vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]->GetHasTarget(
                  arrTT[iTree]));
          if (debug &&
              vvvIsHistFileLeafTreeCustom[iTree][indexNameCustom].size() !=
                  iFile + 1) {
            Fatal("HistMerger::Run",
                  "Wrong size of vvvIsHistFileLeafTreeCustom[%d][%d] (expected "
                  "(iFile+1): %d, got: %ld)",
                  iTree, indexNameCustom, iFile + 1,
                  vvvIsHistFileLeafTreeCustom[iTree][indexNameCustom].size());
          }
          if (debug)
            std::cout
                << vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]
                       ->GetNameLeafModified()
                << " is "
                << (vvvIsHistFileLeafTreeCustom[iTree][indexNameCustom].back()
                        ? ""
                        : " not")
                << " available in this file" << std::endl;
          if (vvvIsHistFileLeafTreeCustom[iTree][indexNameCustom].back()) {
            vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]
                ->EvaluateAndAnalyze(arrTT[iTree]);
          }
        }
      }
      if (debug) std::cout << "Done getting names." << std::endl;
      if (nInfileOpenMax && nInfileOpen >= nInfileOpenMax) {
        closeTFilesIn(iFile + nFileMissingTot, [&areSomeInfilesClosed]() {
          areSomeInfilesClosed = true;
        });
      }
    }
    if (debug)
      std::cout << "nEntryOriginalDatasetCurrent: "
                << nEntryOriginalDatasetCurrent << std::endl;
    if (debug)
      std::cout << "vCrossSection[" << iDataset
                << "]: " << vCrossSection[iDataset] << std::endl;
    Double_t weight =
        (nEntryOriginalDatasetCurrent > 0)
            ? vCrossSection[iDataset] / nEntryOriginalDatasetCurrent
            : 0;
    vWeightDataset.push_back(weight);
    if (debug)
      std::cout << "vWeightDataset[" << iDataset
                << "]: " << vWeightDataset[iDataset] << "," << std::endl;
    if (debug) std::cout << "composed of files index";
    for (UInt_t jFile = 0; jFile < vNumberFile[iDataset]; jFile++) {
      UInt_t indexFile = iFile - vNumberFile[iDataset] + jFile;
      if (debug) std::cout << " " << indexFile;
      vWeightFile.push_back(weight);
    }
    if (debug) std::cout << "." << std::endl;
    if (debug && vWeightDataset.size() != iDataset + 1)
      std::cerr << "Error: vWeightDataset.size() != iDataset + 1 "
                << " (size:" << vWeightDataset.size()
                << ", iDataset:" << iDataset + 1 << ")" << std::endl;
  }
  // tfAutogenHist->Write();

  // Close all the input files if and only if the total number of opened input
  // files are greater than nInfileOpenMax
  if (areSomeInfilesClosed) {
    closeTFilesIn(nFileTotOriginal - 1);
  }

  const std::function<void(TH1 *, LeafAnalyzerAbstract *, TFile *)>
      processHistResult = [getNameHistResultFromNames](
                              TH1 *histResult, LeafAnalyzerAbstract *analyzer,
                              TFile *tfOutHist) {
        histResult->SetName(getNameHistResultFromNames(
            analyzer->GetNameTT(), analyzer->GetNameLeafModified()));
        histResult->SetDirectory(tfOutHist);
        histResult->SetTitle(analyzer->GetNameLeafModified());
        histResult->SetXTitle(analyzer->GetTitleLeaf());
      };

  nLeavesAllTrees = 0;
  for (auto vAnalyzerLeaf : vvAnalyzerLeafTree) {
    UInt_t nLeavesCurrent = vAnalyzerLeaf.size();
    vNLeavesTree.push_back(nLeavesCurrent);
    nLeavesAllTrees += nLeavesCurrent;
  }
  for (auto vAnalyzerLeafCustom : vvAnalyzerLeafTreeCustom) {
    UInt_t nAnalyzersCurrent = vAnalyzerLeafCustom.size();
    vNAnalyzersTreeCustom.push_back(nAnalyzersCurrent);
    nLeavesAllTrees += nAnalyzersCurrent;
  }
  Bool_t isToUseCorrectedTempFile =
      nLeavesToUseCorrectedTempFileMin &&
      nLeavesAllTrees >= nLeavesToUseCorrectedTempFileMin;
  tfCorrectedHist = nullptr;
  if (isToUseCorrectedTempFile) {
    if (debug) std::cout << "Opening tfCorrectedHist ...";
    tfCorrectedHist = TFile::Open(pathTFCorrectedHist,
                                  toRecreateOutFile ? "recreate" : "update");
    if (debug) std::cout << " Done." << std::endl;
  } else {
    vvHistResultLeafTree.clear();
    vvHistResultLeafTree.reserve(nTT);
    for (UInt_t iTree = 0; iTree < nTT; iTree++) {
      std::vector<TH1 *> vHistResultLeaf(
          vNLeavesTree[iTree] + vNAnalyzersTreeCustom[iTree], nullptr);
      vvHistResultLeafTree.push_back(vHistResultLeaf);
    }
    vTFOut.clear();
    vTFOut.reserve(nTT);
    for (UInt_t iTree = 0; iTree < nTT; iTree++) {
      vTFOut.push_back(TFile::Open(funPathTFOutHist(iTree),
                                   toRecreateOutFile ? "recreate" : "update"));
    }
  }
  gStyle->SetOptStat(111111);
  if (debug) {
    std::cout << "Drawing histograms with correct settings ..." << std::endl;
    if (isToUseCorrectedTempFile) {
      std::cout << "Using tfCorrectedHist to store corrected histograms."
                << std::endl;
    } else {
      std::cout << "Not using tempfile tfCorrectedHist, \n"
                << "Hists will be scaled and added on the fly." << std::endl;
    }
  }
  for (auto vAnalyzerLeaf : vvAnalyzerLeafTree) {
    for (LeafAnalyzerAbstract *analyzer : vAnalyzerLeaf) {
      analyzer->Summarize();
    }
  }
  for (auto vAanlyzerLeafCostum : vvAnalyzerLeafTreeCustom) {
    for (LeafAnalyzerAbstract *analyzer : vAanlyzerLeafCostum) {
      analyzer->Summarize();
    }
  }

  std::vector<std::vector<UInt_t>> vvIHistLeafTree(nTT);
  vvIHistLeafTree.clear();
  if (debug) std::cout << "Number of analyzers in each tree: ";
  for (UInt_t iTree = 0; iTree < nTT; iTree++) {
    std::cout << vNLeavesTree[iTree] << "+" << vNAnalyzersTreeCustom[iTree]
              << " ";
    std::vector<UInt_t> vIHistLeaf(
        vNLeavesTree[iTree] + vNAnalyzersTreeCustom[iTree], 0);
    vvIHistLeafTree.push_back(vIHistLeaf);
  }

  if (debug) std::cout << std::endl;
  for (UInt_t iFile = 0, iFileOriginal = 0;
       iFile < nFileTot || iFileOriginal < nFileTotOriginal;
       iFile++, iFileOriginal++) {
    if (debug) std::cout << "iFile: " << iFile;
    TFile *tfInCurrent = nullptr;
    while (iFileOriginal < nFileTotOriginal &&
           !arrIsOpenableFile[iFileOriginal]) {
      iFileOriginal++;
    }
    if (iFileOriginal >= nFileTotOriginal) {
      break;
    }
    if (debug) std::cout << " iFileOriginal: " << iFileOriginal << std::endl;
    if (areSomeInfilesClosed) {
      if (debug) std::cout << "Opening input rootfiles ...";
      arrFile[iFileOriginal] = TFile::Open(funPathTFIn(iFileOriginal));
      nInfileOpen++;
      if (debug) std::cout << " Done." << std::endl;
    }
    tfInCurrent = arrFile[iFileOriginal];
    if (debug)
      std::cout << "tfInCurrent (" << tfInCurrent
                << "), isOpen: " << arrFile[iFileOriginal]->IsOpen()
                << std::endl;
    if (debug) std::cout << "Mounting variables ...";
    mountVarsFromTDir(tfInCurrent);
    if (debug) std::cout << " Done." << std::endl;
    for (UInt_t iTree = 0; iTree < nTT; iTree++) {
      if (arrTT[iTree] == nullptr) {
        continue;
      }
      UInt_t nNameModifiedCurrentTree = vNLeavesTree[iTree];
      UInt_t nAnalyzersCurrentCustom = vNAnalyzersTreeCustom[iTree];
      for (UInt_t indexName = 0, indexNameCustom = 0, indexNameBoth = 0;
           indexName < nNameModifiedCurrentTree ||
           indexNameCustom < nAnalyzersCurrentCustom;
           (indexName < nNameModifiedCurrentTree ? indexName++
                                                 : indexNameCustom++),
                  indexNameBoth = indexName + indexNameCustom) {
        Bool_t isCustom = indexName == nNameModifiedCurrentTree;
        // TList *tlHistCorrected = new TList;
        // std::vector<TString> vNameHistCorrected;
        // vNameHistCorrected.clear();
        // UInt_t iFile=0;
        if (!(isCustom
                  ? vvvIsHistFileLeafTreeCustom[iTree][indexNameCustom][iFile]
                  : vvvIsHistFileLeafTree[iTree][indexName][iFile])) {
          std::cerr << "(Tree index, Leaf index) (" << iTree << ", "
                    << indexNameBoth << ") not found in file index " << iFile
                    << " ("
                    << (isCustom
                            ? vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]
                            : vvAnalyzerLeafTree[iTree][indexName])
                           ->GetNameLeafModified()
                    << ")" << std::endl;
          continue;
        }
        if (debug)
          std::cout << "vvIHistLeafTree[" << iTree << "][" << indexNameBoth
                    << "]: " << vvIHistLeafTree[iTree][indexNameBoth]
                    << std::endl;
        if (debug)
          std::cout << "isCustom: " << isCustom << ", indexName: " << indexName
                    << ", indexNameCustom: " << indexNameCustom << std::endl;
        TH1 *histCorrected = nullptr;
        if (debug) std::cout << "Getting analyzer ...";
        LeafAnalyzerAbstract *&analyzer =
            (isCustom ? vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]
                      : vvAnalyzerLeafTree[iTree][indexName]);
        if (debug) std::cout << " (" << analyzer << ") Done." << std::endl;

        // TString nameLeaf =
        //     analyzer->GetVNameLeafFile()[vvIHistLeafTree[iTree][indexNameBoth]];
        if (debug)
          std::cout << "nameLeafModified: " << analyzer->GetNameLeafModified()
                    << std::endl;
        TString tstrHistSetting = analyzer->GetHistSetting();
        if (!analyzer->GetIsEverAnalyzed()) {
          TString tstrWarningCommon = TString::Format(
              "%s (%d), %s (%d) (%s) has never been analyzed",
              analyzer->GetNameTT().Data(), iTree,
              analyzer->GetNameLeafModified().Data(), indexNameBoth,
              analyzer->GetExpressionBeforeSetting().Data());
          if (analyzer->GetAllowNeverAnalyzed()) {
            std::cerr << tstrWarningCommon
                      << "!\n"
                         "It will not be saved."
                      << std::endl;
          } else {
            Fatal("HistMerger::Run", "%s", (tstrWarningCommon + "!").Data());
          }
        } else if (!analyzer->GetDontCheckEmptyness() &&
                   analyzer
                       ->GetVIsEmpty()[vvIHistLeafTree[iTree][indexNameBoth]]) {
          if (vvIHistLeafTree[iTree][indexNameBoth] == 0 &&
              analyzer->GetAreAllEmpty()) {
            if (debug)
              std::cout << "Got " << analyzer->GetNameLeafModified()
                        << " (all-empty leaf)." << std::endl;
            histCorrected = analyzer->GetHistEmptyPreferred();
            TString nameHistCorrected =
                (isToUseCorrectedTempFile
                     ? getNameHistOfFileFromAnalyzerAndIFile(
                           analyzer, 0, "CorrectedAllEmpty")
                     : getNameHistResultFromNames(
                           analyzer->GetNameTT(),
                           analyzer->GetNameLeafModified()));
            if (histCorrected == nullptr) {
              // arrTT[iTree]->Draw(nameHistCorrected +
              //                     (tstrHistSetting.Length() == 0
              //                         ? ""
              //                         : ">>" + tstrHistSetting));
              // histCorrected = (TH1 *)gDirectory->Get(nameHistCorrected);
              histCorrected = analyzer->DrawHistCorrected(
                  nameHistCorrected, arrTT[iTree],
                  isCustom ? -1 : vvIHistLeafTree[iTree][indexNameBoth]);
              histCorrected->Scale(vWeightFile[iFile]);
            }
            if (isToUseCorrectedTempFile) {
              if (debug)
                std::cout << "Writing " << nameHistCorrected << " ("
                          << histCorrected << ") to tfCorrectedHist ...";
              tfCorrectedHist->cd();
              histCorrected->SetDirectory(tfCorrectedHist);
              histCorrected->SetName(nameHistCorrected);
              histCorrected->Write(nameHistCorrected);
              if (debug) std::cout << " Done." << std::endl;
            } else {
              processHistResult(histCorrected, analyzer, vTFOut[iTree]);
              vvHistResultLeafTree[iTree][indexNameBoth] = histCorrected;
            }
            // break; // Lo the bug.
          }
          if (debug && !analyzer->GetAreAllEmpty())
            std::cerr << "Empty histogram found for (" << vNameTT[iTree] << ", "
                      << analyzer->GetNameLeafModified() << ") (" << iTree
                      << ", " << indexNameBoth << ")" << std::endl;
        } else {
          TString nameHistCorrected = getNameHistOfFileFromAnalyzerAndIFile(
              analyzer, iFile, "Corrected");
          // arrTT[iTree]->Draw(vvvNameFileLeafTree[iTree][indexName][iFile] +
          // ">>" +
          //                    nameHistCorrected +
          //                    vvHistsettingLeafTree[iTree][indexName]);

          // arrTT[iTree]->Draw(nameLeaf + ">>" + nameHistCorrected +
          //                    tstrHistSetting);
          // if (debug) std::cout << " Getting ...";
          // histCorrected = (TH1 *)gDirectory->Get(nameHistCorrected);
          histCorrected = analyzer->DrawHistCorrected(
              nameHistCorrected, arrTT[iTree],
              isCustom ? -1 : vvIHistLeafTree[iTree][indexNameBoth]);
          if (debug) std::cout << " (" << histCorrected << ") ";
          if (debug) std::cout << " Scaling ...";
          histCorrected->Scale(vWeightFile[iFile]);
          if (debug) std::cout << " Done." << std::endl;
          if (isToUseCorrectedTempFile) {
            if (debug)
              std::cout << "Writing to tfCorrectedHist (" << tfCorrectedHist
                        << ") ...";
            histCorrected->SetDirectory(tfCorrectedHist);
            histCorrected->SetName(nameHistCorrected);
            histCorrected->Write(nameHistCorrected);
            if (debug) std::cout << " Done." << std::endl;
          } else {
            if (vvHistResultLeafTree[iTree][indexNameBoth] == nullptr) {
              if (debug)
                std::cout << "Storing to vvHistResultLeafTree[" << iTree << "]["
                          << indexNameBoth << "]";
              ;
              TString nameHistResult = getNameHistResultFromNames(
                  analyzer->GetNameTT(), analyzer->GetNameLeafModified());
              vvHistResultLeafTree[iTree][indexNameBoth] =
                  (TH1 *)histCorrected->Clone(nameHistResult);
              processHistResult(vvHistResultLeafTree[iTree][indexNameBoth],
                                analyzer, vTFOut[iTree]);
              if (debug)
                std::cout << " (" << vvHistResultLeafTree[iTree][indexNameBoth]
                          << ")" << std::endl;
              if (debug) std::cout << " Done." << std::endl;
            } else {
              if (debug)
                std::cout << "Adding to vvHistResultLeafTree[" << iTree << "]["
                          << indexNameBoth << "]";
              histCorrected->SetName(nameHistCorrected);
              vTFOut[iTree]->cd();
              vvHistResultLeafTree[iTree][indexNameBoth]->Add(
                  (TH1 *)histCorrected->Clone(nameHistCorrected + "ToAdd"));
              if (debug) std::cout << " Done." << std::endl;
            }
          }
        }
        vvIHistLeafTree[iTree][indexNameBoth]++;
      }
    }
    if (areSomeInfilesClosed && nInfileOpen >= nInfileOpenMax) {
      closeTFilesIn(iFileOriginal);
    }
  }
  if (nInfileOpenMax) {
    if (debug)
      std::cout << "Closing openable input ROOT files ..." << std::endl;
    closeTFilesIn(nFileTotOriginal - 1);
    if (debug) std::cout << "Done." << std::endl;
  }

  if (isToUseCorrectedTempFile) {
    //   if (debug) std::cout << "Closing tfCorrectedHist ...";
    //   tfCorrectedHist->Close();
    //   delete tfCorrectedHist;
    //   if (debug) std::cout << " Done." << std::endl;
    //   tfCorrectedHist = TFile::Open(pathTFCorrectedHist, "read");
    tfCorrectedHist->ReOpen("read");
  }
  // tfAutogenHist->ReOpen("read");

  // if (debug) std::cout << "Closing openable input ROOT files ..." <<
  // std::endl; for (UInt_t iFileOriginal = 0; iFileOriginal < nFileTotOriginal;
  //      iFileOriginal++) {
  //   if (arrIsOpenableFile[iFileOriginal]) {
  //     if (debug) std::cout << "Closing iFileOriginal: " << iFileOriginal << "
  //     "; arrFile[iFileOriginal]->Close();
  //   }
  // }
  // if (debug) std::cout << "Done." << std::endl;

  if (debug) std::cout << "Saving results ..." << std::endl;
  if (debug && isToUseCorrectedTempFile)
    std::cout << "Using tfCorrectedHist,\n"
              << "will merge before save" << std::endl;
  // tfAutogenHist = TFile::Open(pathTFAutogenHist, "read");
  // tfCorrectedHist = TFile::Open(pathTFCorrectedHist, "read");
  for (UInt_t iTree = 0; iTree < nTT; iTree++) {
    if (debug)
      std::cout << "iTree: " << iTree << " (" << vNameTT[iTree] << ")"
                << std::endl;
    TFile *tfOutHistMerge = nullptr;
    if (isToUseCorrectedTempFile) {
      tfOutHistMerge = TFile::Open(funPathTFOutHist(iTree),
                                   toRecreateOutFile ? "recreate" : "update");
      if (debug)
        std::cout << "tfOutHist: " << tfOutHistMerge
                  << " (isOpen: " << tfOutHistMerge->IsOpen() << ")"
                  << std::endl;
    } else {
      vTFOut[iTree]->cd();
    }
    for (UInt_t indexName = 0, indexNameCustom = 0, indexNameBoth = 0;
         indexName < vNLeavesTree[iTree] ||
         indexNameCustom < vNAnalyzersTreeCustom[iTree];
         (indexName < vNLeavesTree[iTree] ? indexName++ : indexNameCustom++),
                indexNameBoth = indexName + indexNameCustom) {
      Bool_t isCustom = indexName == vNLeavesTree[iTree];
      TH1 *histResult;
      LeafAnalyzerAbstract *analyzer =
          isCustom ? vvAnalyzerLeafTreeCustom[iTree][indexNameCustom]
                   : vvAnalyzerLeafTree[iTree][indexName];
      TString nameHistResult = getNameHistResultFromNames(
          analyzer->GetNameTT(), analyzer->GetNameLeafModified());
      if (debug)
        std::cout << "Generating " << nameHistResult << " ..." << std::endl;
      if (!analyzer->GetIsEverAnalyzed()) {
        Warning("HistMerger::Run", "%s is never analyzed and will not be saved. tree: %s (%d), indexName: %d, indexNameCustom: %d",
        analyzer->GetNameLeafModified().Data(), analyzer->GetNameTT().Data(), iTree, indexName, indexNameCustom);
        continue;
      }
      if (!isToUseCorrectedTempFile) {
        if (debug) std::cout << "Writing histogram " << vvHistResultLeafTree[iTree][indexNameBoth] << " " << vvHistResultLeafTree[iTree][indexNameBoth]->GetName() << " ...";
        vvHistResultLeafTree[iTree][indexNameBoth]->Write();
        if (debug) std::cout << " Done." << std::endl;
        // histResult = vvHistResultLeafTree[iTree][indexName];
      } else if (analyzer->GetAreAllEmpty()) {
        // UInt_t iFileFirst;
        // for (iFileFirst = 0;
        //      !vvvIsHistFileLeafTree[iTree][indexName][iFileFirst];
        //      iFileFirst++)
        //   ;
        // histResult = (TH1 *)tfAutogenHist
        //                  ->Get(getNameHistOfFile(iTree, indexName,
        //                  iFileFirst,
        //                                          "Autogen"))
        //                  ->Clone(nameHistResult);
        // histResult->SetName(nameHistResult);
        // if (debug)
        //   std::cout << "Got " << nameHistResult << " (empty histogram)."
        //             << std::endl;
        TString nameHistCorrected = getNameHistOfFileFromAnalyzerAndIFile(
            analyzer, 0, "CorrectedAllEmpty");
        histResult = (TH1 *)tfCorrectedHist->Get(nameHistCorrected);
        if (debug) std::cout << "histResult: " << histResult;
        processHistResult(histResult, analyzer, tfOutHistMerge);
        if (debug)
          std::cout << " (name: " << histResult->GetName() << ")" << std::endl;
        histResult->Write(nameHistResult);
      } else {
        // Bool_t isErrorUnexpected = false;
        TList *tlHistCorrectedToMerge = new TList;
        UInt_t nHist = 0;
        for (UInt_t iFile = 0; iFile < nFileTot; iFile++) {
          if (!(isCustom
                    ? vvvIsHistFileLeafTreeCustom[iTree][indexNameCustom][iFile]
                    : vvvIsHistFileLeafTree[iTree][indexName][iFile])) {
            continue;
          }
          if (!analyzer->GetDontCheckEmptyness() &&
              analyzer->GetVIsEmpty()[nHist]) {
            nHist++;
            continue;
          }
          if (debug) std::cout << "iFile: " << iFile << " ";
          TString nameHistCorrected = getNameHistOfFileFromAnalyzerAndIFile(
              analyzer, iFile, "Corrected");
          TH1 *histCorrected = (TH1 *)tfCorrectedHist->Get(nameHistCorrected);
          if (histCorrected == nullptr) {
            std::cerr << "Fatal: Histogram not found: " << nameHistCorrected
                      << std::endl;
            // isErrorUnexpected = true;
            break;
          }
          tlHistCorrectedToMerge->Add(histCorrected->Clone(nameHistCorrected));
          nHist++;
        }
        if (!nHist) {
          // std::cerr << "Fatal: "
          //           << "No histograms found for (" << iTree << ", " <<
          //           indexName
          //           << ") (" << vvNameModifiedLeafTree[iTree][indexName] <<
          //           ")"
          //           << std::endl;
          histResult = nullptr;
          // isErrorUnexpected = true;
        } else {
          if (debug) std::cout << "Merging ...";
          histResult =
              (TH1 *)tlHistCorrectedToMerge->First()->Clone(nameHistResult);
          processHistResult(histResult, analyzer, tfOutHistMerge);
          histResult->Clear();
          histResult->Merge(tlHistCorrectedToMerge);
          processHistResult(histResult, analyzer, tfOutHistMerge);
          histResult->Write(nameHistResult);
        }
        // if (isErrorUnexpected) continue;
        if (debug) std::cout << " Done." << std::endl;
        tlHistCorrectedToMerge->Delete();
        delete tlHistCorrectedToMerge;
      }
    }
    if (isToUseCorrectedTempFile) {
      tfOutHistMerge->Close();
    } else {
      vTFOut[iTree]->Close();
    }
  }
  if (isToUseCorrectedTempFile) {
    tfCorrectedHist->Close();
  }
  // tfAutogenHist->Close();
  for (auto analyzerLeaf : vvAnalyzerLeafTree) {
    for (auto analyzer : analyzerLeaf) {
      analyzer->Finalize();
    }
  }
  for (auto analyzerLeaf : vvAnalyzerLeafTreeCustom) {
    for (auto analyzer : analyzerLeaf) {
      analyzer->Finalize();
    }
  }
}

#endif
