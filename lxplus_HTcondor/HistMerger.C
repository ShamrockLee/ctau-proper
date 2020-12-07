#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>

#include "mergeToHists.h"

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

#include <functional>
#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

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

#ifndef TSTRING_UTILS
#define TSTRING_UTILS
Bool_t startsWithOfTString(const TString target, const TString pattern,
                           const Ssiz_t displacement = 0) {
  return target.SubString(pattern).Start() == displacement;
}

// Bool_t endsWithOfTString(const TString target, const TString pattern, const
// Ssiz_t displacement=0) {
//   return target.SubString(pattern).Start() == target.Length() -
//   pattern.Length() - displacement;
// }

Bool_t containsOfTString(const TString target, const TString pattern) {
  return target.SubString(pattern).Start() != -1;
}
#endif

#ifndef HISTMERGER_C
#define HISTMERGER_C

class HistMerger::LeafAnalyzerDefault : LeafAnalyzerAbstract {
 public:
  static inline TString GetNameLeafModified(TString nameLeaf) {
    TString nameLeafNew = TString(nameLeaf);
    Ssiz_t indexLeftSquareBracket = nameLeaf.First('[');
    if (indexLeftSquareBracket != -1) {
      nameLeafNew.Resize(indexLeftSquareBracket);
    }
    return nameLeafNew;
  }
  static inline TString GetNameLeafModified(TLeaf *leaf) {
    return GetNameLeafModified((TString)leaf->GetName());
  }
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
  LeafAnalyzerDefault() {
    this->vNameLeafFile.clear();
    this->nameLeafModified = "";
    this->titleLeaf = "";
    this->typeNameLeaf = "";
    this->vIsEmpty.clear();
    this->areAllEmpty = true;
    this->nBinsCorrect = 0;
    this->lowerCorrect = 0.;
    this->upperCorrect = 0.;
    this->debug = false;
    this->nameTT = "";
    this->adjustHistSettingPerLeafTreeExtra = nullptr;
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
  void SetFunAdjustHistSettingExtra(
      std::function<void(Int_t &NBinCorrect, Double_t &lowerCorrect,
                         Double_t &upperCorrect, TString nameTT,
                         TString nameLeafModified, TString typeNameLeaf,
                         TString titleLeaf)>
          adjustHistSettingPerLeafTreeExtra) {
    this->adjustHistSettingPerLeafTreeExtra = adjustHistSettingPerLeafTreeExtra;
  }
  void AnalyzeLeafBasic(TLeaf *leaf) {
    TString nameLeaf = leaf->GetName();
    this->nameLeafModified = GetNameLeafModified(nameLeaf);
    this->titleLeaf = GetTitleLeaf(leaf);
    this->typeNameLeaf = GetTypeNameLeaf(leaf);
    this->AnalyzeLeafBasicExtra();
  }
  TString GetNameLeafModified() { return this->nameLeafModified; }
  TString GetTitleLeaf() { return this->titleLeaf; }
  TString GetTypeNameLeaf() { return this->typeNameLeaf; }
  static inline Int_t GetNBins(TH1 *hist) { return hist->GetNbinsX(); }
  static inline Double_t GetMinimumValue(TH1 *hist) {
    return hist->GetBinCenter(hist->FindFirstBinAbove()) -
           hist->GetBinWidth(1) / 2;
  }
  static inline Double_t GetMaximumValue(TH1 *hist) {
    return hist->GetBinCenter(hist->FindLastBinAbove()) +
           hist->GetBinWidth(1) / 2;
  }
  void AnalyzeLeaf(TLeaf *leaf) {
    vNameLeafFile.push_back(leaf->GetName());
    TString nameHistAutogen = "h" + this->nameLeafModified + "Autogen" + "Temp";
    leaf->Draw(nameHistAutogen);
    TH1 *histAutogen = (TH1 *)gDirectory->Get(nameHistAutogen);
    Int_t nEntries = histAutogen->GetEntries();
    vNEntriesFile.push_back(nEntries);
    if (nEntries == 0) {
      // Empty histogram
      vIsEmpty.push_back(true);
      if (isFirstNonEmptyLeaf && tstrHistSetting.Length() != 0) {
        if (histTypeFlavor == HistTypeFlavor::Boolean) {
          lowerCorrect = 0;
          upperCorrect = 2;
          nBinsCorrect = 2;
          tstrHistSetting = TString::Format("(%d,%d,%d)", nBinsCorrect,
                                            lowerCorrect, upperCorrect);
        } else {
          nBinsCorrect = GetNBins(histAutogen);
          lowerCorrect = GetMinimumValue(histAutogen);
          upperCorrect = GetMaximumValue(histAutogen);
          tstrHistSetting = TString::Format("(%d,%F,%F)", nBinsCorrect,
                                            lowerCorrect, upperCorrect);
        }
      }
    } else {
      vIsEmpty.push_back(false);
      areAllEmpty = false;
      if (histTypeFlavor == HistTypeFlavor::Boolean) {
      } else {
        Double_t minimumValue = GetMinimumValue(histAutogen);
        Double_t maximumValue = GetMaximumValue(histAutogen);
        if (isFirstNonEmptyLeaf) {
          minimumValueAllFiles == minimumValue;
          maximumValueAllFiles == maximumValue;
          isFirstNonEmptyLeaf = false;
        } else {
          if (minimumValueAllFiles > minimumValue)
            minimumValueAllFiles = minimumValue;
          if (maximumValueAllFiles < maximumValue)
            maximumValueAllFiles = maximumValue;
        }
        if (histTypeFlavor == HistTypeFlavor::Integer) {
        } else {
          Int_t nBins = histAutogen->GetNbinsX();
          vNBinsFile.push_back(nBins);
          vLowerFile.push_back(histAutogen->GetBinCenter(1) -
                               histAutogen->GetBinWidth(1) / 2);
          vUpperFile.push_back(histAutogen->GetBinCenter(nBins) +
                               histAutogen->GetBinWidth(1) / 2);
        }
      }
    }
    histAutogen->Delete();
    delete histAutogen;
  }
  std::vector<TString> &GetVNameLeafFile() { return this->vNameLeafFile; }
  std::vector<Bool_t> &GetVIsEmpty() { return this->vIsEmpty; }
  void Summarize() {
    if (debug && histTypeFlavor != HistTypeFlavor::Boolean)
      std::cout << "minimumValueAllFiles, maximumValueAllFiles: "
                << minimumValueAllFiles << ", " << maximumValueAllFiles
                << std::endl;
    if (histTypeFlavor == HistTypeFlavor::Boolean) {
      nBinsCorrect = 2;
      lowerCorrect = 0;
      upperCorrect = 2;
      tstrHistSetting = TString::Format("(%d,%d,%d)", nBinsCorrect,
                                        lowerCorrect, upperCorrect);
    } else if (histTypeFlavor == HistTypeFlavor::Integer) {
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
        ULong64_t l64NBinsCorrect = upperCorrect - lowerCorrect - 1;
        UShort_t nbitsReduced = 0;
        while (l64NBinsCorrect & ~((ULong64_t)__INT32_MAX__)) {
          nbitsReduced++;
          l64NBinsCorrect >> 1;
        }
        nBinsCorrect = l64NBinsCorrect;
      } else {
        nBinsCorrect = upperCorrect - lowerCorrect;
      }
      tstrHistSetting =
          TString::Format("(%d,%lld,%lld)", nBinsCorrect,
                          (Long64_t)lowerCorrect, (Long64_t)upperCorrect);
    } else {
      UInt_t nHists = vNBinsFile.size();
      nBinsCorrect = vNBinsFile[nHists >> 1];
      if (minimumValueAllFiles < 0 && maximumValueAllFiles > 0) {
        Int_t nDigitsM1Minimum = getNDigitsM1(minimumValueAllFiles);
        Int_t nDigitsM1Maximum = getNDigitsM1(maximumValueAllFiles);
        lowerCorrect = intpow(
            10, nDigitsM1Minimum,
            TMath::Floor(intpow(10, -nDigitsM1Minimum, minimumValueAllFiles)));
        upperCorrect = intpow(
            10, nDigitsM1Maximum,
            TMath::Ceil(intpow(10, -nDigitsM1Maximum, maximumValueAllFiles)));
      } else if (minimumValueAllFiles == 0) {
        lowerCorrect = minimumValueAllFiles;
        upperCorrect = maximumValueAllFiles;
        if (lowerCorrect == upperCorrect) {
          upperCorrect += __FLT_EPSILON__;
        }
      } else if (maximumValueAllFiles == 0) {
        lowerCorrect = minimumValueAllFiles;
        upperCorrect = maximumValueAllFiles;
        if (lowerCorrect == upperCorrect) {
          lowerCorrect -= __FLT_EPSILON__;
        }
      } else {
        Bool_t isPositive = maximumValueAllFiles > 0;
        Double_t valueOuter =
            isPositive ? maximumValueAllFiles : minimumValueAllFiles;
        Double_t valueInner =
            isPositive ? minimumValueAllFiles : maximumValueAllFiles;
        Int_t nDigitsM1Outer = getNDigitsM1(valueOuter);
        Int_t orderDiff = TMath::Nint(
            TMath::Log10(maximumValueAllFiles - minimumValueAllFiles));
        Int_t mDigitsToRound = TMath::Min(nDigitsM1Outer, orderDiff);
        valueOuter =
            intpow(10, mDigitsToRound,
                   TMath::Floor(intpow(10, -mDigitsToRound, valueOuter)));
        valueInner =
            intpow(10, mDigitsToRound,
                   TMath::Ceil(intpow(10, -mDigitsToRound, valueInner)));
        if (valueInner == valueOuter) {
          valueOuter += isPositive ? __FLT_EPSILON__ : -FLT_EPSILON;
        }
        if (isPositive) {
          lowerCorrect = valueInner;
          upperCorrect = valueOuter;
        } else {
          lowerCorrect = valueOuter;
          upperCorrect = valueInner;
        }
      }
    }
    if (adjustHistSettingPerLeafTreeExtra != nullptr) {
      adjustHistSettingPerLeafTreeExtra(nBinsCorrect, lowerCorrect,
                                        upperCorrect, nameTT, nameLeafModified,
                                        typeNameLeaf, titleLeaf);
    }
  }

 protected:
  Bool_t debug;
  TString nameTT;
  std::function<void(Int_t &NBinCorrect, Double_t &lowerCorrect,
                     Double_t &upperCorrect, TString nameTT,
                     TString nameLeafModified, TString typeNameLeaf,
                     TString titleLeaf)>
      adjustHistSettingPerLeafTreeExtra;
  std::vector<TString> vNameLeafFile;
  TString nameLeafModified;
  TString titleLeaf;
  TString typeNameLeaf;
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

  void AnalyzeLeafBasicExtra() {
    // isTypeFloatingPoint = containsOfTString(typeNameLeaf, "loat") ||
    //                containsOfTString(typeNameLeaf, "ouble");
    if (containsOfTString(typeNameLeaf, "Bool") ||
        containsOfTString(typeNameLeaf, "bool")) {
      histTypeFlavor = HistTypeFlavor::Boolean;
    } else if (containsOfTString(typeNameLeaf, "loat") ||
               containsOfTString(typeNameLeaf, "ouble")) {
      histTypeFlavor = HistTypeFlavor::FloatingPoint;
    } else {
      histTypeFlavor = HistTypeFlavor::Integer;
    }
  }
};

void HistMerger::InitializeHidden() {
  this->vWeightDataset.clear();
  this->vWeightFile.clear();
  this->vvvIsHistFileLeafTree.clear();
  vvHistResultLeafTree.clear();
  vNLeavesTree.clear();
  this->vvNameModifiedLeafTree.clear();
  this->vvAnalyzerLeafTree.clear();
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

  const UInt_t nTT = vNameTT.size();

  TTree *arrTT[nTT];

  auto mountVarsFromTDir = [this, &nEntryOriginal, &arrTT,
                            nTT](TDirectory *tdir) {
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
  auto getNameHistOfFile = [this](
                               UInt_t iTree, UInt_t indexName, UInt_t iFile,
                               TString suffixStage) -> TString {
    return (TString) "h" + vNameTT[iTree] +
           vvNameModifiedLeafTree[iTree][indexName] + suffixStage + iFile;
  };

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
  std::function<TString(UInt_t, UInt_t)> getNameHistResult =
      [this](UInt_t iTree, UInt_t indexName) -> TString {
    return (TString) "h" + vvNameModifiedLeafTree[iTree][indexName];
  };
  if (dirTFOut.Length() > 1 && dirTFOut.EndsWith(seperatorPath)) {
    dirTFOut.Resize(dirTFOut.Length() - 1);
  }
  std::function<TString(UInt_t)> funPathTFOutHist =
      [this](UInt_t iTree) -> TString {
    return dirTFOut + seperatorPath + funNameTFOut(vNameTT[iTree]);
  };
  std::function<void(Int_t &, Double_t &, Double_t &, UInt_t, UInt_t)>
      adjustHistSettingPerLeafTreeExtraIndex =
          [this](
              Int_t &nBinCorrect, Double_t &lowerCorrect,
              Double_t &upperCorrect, UInt_t iTree, UInt_t indexName) {
            TString nameTT = vNameTT[iTree];
            TString &nameLeafModified =
                vvNameModifiedLeafTree[iTree][indexName];
            // TString typeName = vvTypeNameLeafTree[iTree][indexName];
            // TString title = vvTitleLeafTree[iTree][indexName];
            LeafAnalyzerAbstract *&analyzer =
                vvAnalyzerLeafTree[iTree][indexName];
            TString typeName = analyzer->GetTypeNameLeaf();
            TString title = analyzer->GetTitleLeaf();
            return adjustHistSettingPerLeafTreeExtra(
                nBinCorrect, lowerCorrect, upperCorrect, nameTT,
                nameLeafModified, typeName, title);
          };
  for (UInt_t iTree = 0; iTree < nTT; iTree++) {
    std::vector<std::vector<Bool_t>> vvIsHistFileLeaf;
    vvIsHistFileLeaf.clear();
    vvvIsHistFileLeafTree.push_back(vvIsHistFileLeaf);
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
        if (allowMissing) {
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
          TString nameLeafModified = getNameLeafModified(nameLeaf);
          if (getIsToVetoLeaf != nullptr && getIsToVetoLeaf(vNameTT[iTree], nameLeafModified)) continue;
          LeafAnalyzerAbstract *analyzer;
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
          if (debug) std::cout << std::endl;
          if (!isIncluded) {
            if (debug) std::cout << "Found new leafname." << std::endl;
            if (debug)
              std::cout << "nameLeaf: " << nameLeaf
                        << ", nameLeafModified: " << nameLeafModified
                        << std::endl;
            vvNameModifiedLeafTree[iTree].push_back(nameLeafModified);
            analyzer = supplyLeafAnalyzer();
            vvAnalyzerLeafTree[iTree].push_back(analyzer);
            analyzer->SetDebug(debug);
            analyzer->SetNameTT(vNameTT[iTree]);
            analyzer->SetFunAdjustHistSettingExtra(
                adjustHistSettingPerLeafTreeExtra);
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
            if (indexName >= vvNameModifiedLeafTree[iTree].size()) {
              std::cerr << "nameLeafModified duplication (" << nameLeaf << ", "
                        << nameLeafModified << ")" << std::endl;
              std::cerr << " indexName of the original leaf: " << indexName
                        << std::endl;
            } else if (arrIsHistLeafOld[indexName]) {
              std::cerr << "nameLeafModified duplication (" << nameLeaf << ", "
                        << nameLeafModified << ")" << std::endl;
              std::cerr << "indexName of the original leaf: " << indexName
                        << std::endl;
            } else {
              arrIsHistLeafOld[indexName] = true;
              analyzer = vvAnalyzerLeafTree[iTree][indexName];
            }
          }
          // vvvIsHistFileLeafTree[iTree][indexName][iFile] = true;

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

  nLeavesAllTrees = 0;
  for (auto vAnalyzerLeaf : vvAnalyzerLeafTree) {
    UInt_t nLeavesCurrent = vAnalyzerLeaf.size();
    vNLeavesTree.push_back(nLeavesCurrent);
    nLeavesAllTrees += nLeavesCurrent;
  }
  Bool_t isToUseCorrectedTempFile =
      nLeavesToUseCorrectedTempFileMin &&
      nLeavesAllTrees >= nLeavesToUseCorrectedTempFileMin;
  if (isToUseCorrectedTempFile) {
    tfCorrectedHist = TFile::Open(pathTFCorrectedHist,
                                  toRecreateOutFile ? "recreate" : "update");
    tfCorrectedHist->Close();
    delete tfCorrectedHist;
  } else {
    for (UInt_t nLeavesCurrent : vNLeavesTree) {
      std::vector<TH1 *> vHistResultLeaf(nLeavesCurrent, nullptr);
      vvHistResultLeafTree.push_back(vHistResultLeaf);
    }
  }
  // tfAutogenHist = TFile::Open(dirCondorPackCurrent + seperatorPath +
  // "output_" + nameDatagroup + "_" + "Autogen" + "_" + nameClusterID +
  // "_hist.root", "read");
  gStyle->SetOptStat(111111);
  // std::vector<std::vector<Bool_t>> vvAreAllEmptyLeafTree(
  //     nTT);  //< Whether all the histograms have zero entry of each leaf in
  //     each
  //            // tree
  // vvAreAllEmptyLeafTree.clear();
  // std::vector<std::vector<TString>> vvHistsettingLeafTree(
  //     nTT);  //< The histogram settings "(binnumber, lower, upper)" of each
  //     leaf
  //            // in each tree;
  // vvHistsettingLeafTree.clear();
  // std::vector<std::vector<std::vector<TString>>> vvvNameFileLeafTree(nTT);
  // vvvNameFileLeafTree.clear();
  // for (UInt_t iTree = 0; iTree < nTT; iTree++) {
  //   {
  //     std::vector<Bool_t> vAreAllEmptyLeaf;
  //     vAreAllEmptyLeaf.clear();
  //     vvAreAllEmptyLeafTree.push_back(vAreAllEmptyLeaf);
  //     std::vector<TString> vHistsettingLeaf;
  //     vHistsettingLeaf.clear();
  //     vvHistsettingLeafTree.push_back(vHistsettingLeaf);
  //     std::vector<std::vector<TString>> vvNameFileLeaf;
  //     vvNameFileLeaf.clear();
  //     vvvNameFileLeafTree.push_back(vvNameFileLeaf);
  //   }
  //   UInt_t nLeaf = vvNameModifiedLeafTree[iTree].size();
  //   for (UInt_t indexName = 0; indexName < nLeaf; indexName++) {
  //     if (debug) std::cout << "indexName: " << indexName << std::endl;

  //     // Index of the first true element of
  //     // vvvIsHistFileLeafTree[iTree][indexName]
  //     UInt_t iFileFirst = std::distance(
  //         vvvIsHistFileLeafTree[iTree][indexName].begin(),
  //         std::find(vvvIsHistFileLeafTree[iTree][indexName].begin(),
  //                   vvvIsHistFileLeafTree[iTree][indexName].end(), true));
  //     tfAutogenHist->ReOpen("read");
  //     TH1 *histFirst = (TH1 *)tfAutogenHist->Get(
  //         getNameHistOfFile(iTree, indexName, iFileFirst, "Autogen"));
  //     UInt_t nHistAllFile =
  //         std::count(vvvIsHistFileLeafTree[iTree][indexName].begin(),
  //                    vvvIsHistFileLeafTree[iTree][indexName].end(), true);
  //     if (debug) std::cout << "nHistAllFile: " << nHistAllFile << std::endl;
  //     Int_t arrNEntryFile[nHistAllFile];
  //     Double_t arrMinimumFile[nHistAllFile];
  //     Double_t arrMaximumFile[nHistAllFile];
  //     Int_t arrNBinFile[nHistAllFile];
  //     Double_t arrLowerFile[nHistAllFile];
  //     Double_t arrUpperFile[nFileTot];
  //     Bool_t areAllEmpty = true;
  //     // UInt_t iHist = 0;
  //     if (debug)
  //       std::cout << "Analyzing autogen histograms for "
  //                 << vvNameModifiedLeafTree[iTree][indexName] << "...";
  //     // for (auto histFileRaw: *tlHistFile) {
  //     std::vector<TString> vLeafnameFile;
  //     vLeafnameFile.clear();
  //     for (UInt_t iFile = 0, iHist = 0; iFile < nFileTot; iFile++) {
  //       if (!vvvIsHistFileLeafTree[iTree][indexName][iFile]) {
  //         vLeafnameFile.push_back("");
  //         if (debug) std::cout << " Skipping iFile: " << iFile;
  //         continue;
  //       }
  //       if (debug) std::cout << " iHist: " << iHist;
  //       // TH1 *histFile = (TH1 *)histFileRaw;
  //       TString nameHistAutogen =
  //           getNameHistOfFile(iTree, indexName, iFile, "Autogen");
  //       // TH1 *histAutogen = (TH1 *) gDirectory->Get(nameHistAutogen);
  //       TH1 *histAutogen = (TH1 *)tfAutogenHist->Get(nameHistAutogen);
  //       // TH1 *histAutogen = vvvHistFileLeafTree[iTree][indexName][iHist];
  //       if (debug && histAutogen == nullptr)
  //         std::cerr << "Fatal: histAutogen from " << nameHistAutogen
  //                   << " is nullptr!" << std::endl;
  //       // histFile->Scale(vWeightFile[iFile]);
  //       Int_t nEntry = histAutogen->GetEntries();
  //       if (nEntry > 0) {
  //         areAllEmpty = false;
  //       }
  //       arrNEntryFile[iHist] = nEntry;
  //       arrMinimumFile[iHist] =
  //           histAutogen->GetBinCenter(histAutogen->FindFirstBinAbove()) -
  //           histAutogen->GetBinWidth(1) / 2;
  //       arrMaximumFile[iHist] =
  //           histAutogen->GetBinCenter(histAutogen->FindLastBinAbove()) +
  //           histAutogen->GetBinWidth(1) / 2;
  //       Int_t nBin = histAutogen->GetNbinsX();
  //       arrNBinFile[iHist] = nBin;
  //       arrLowerFile[iHist] =
  //           histAutogen->GetBinCenter(1) - histAutogen->GetBinWidth(1) / 2;
  //       arrUpperFile[iHist] =
  //           histAutogen->GetBinCenter(nBin) + histAutogen->GetBinWidth(1) /
  //           2;
  //       vLeafnameFile.push_back(histAutogen->GetTitle());
  //       iHist++;
  //     }
  //     vvvNameFileLeafTree[iTree].push_back(vLeafnameFile);
  //     if (debug) std::cout << " Done." << std::endl;
  //     vvAreAllEmptyLeafTree[iTree].push_back(areAllEmpty);
  //     if (areAllEmpty) {
  //       if (debug)
  //         std::cout << "Empty leaf encountered: "
  //                   << vvNameModifiedLeafTree[iTree][indexName] << std::endl;
  //       vvHistsettingLeafTree[iTree].push_back("");
  //     } else {
  //       if (debug)
  //         std::cout << "Calculating histogram settings ..." << std::endl;
  //       Double_t lowerCorrect;
  //       Double_t upperCorrect;
  //       Int_t nBinCorrect;
  //       if (containsOfTString(vvTypeNameLeafTree[iTree][indexName], "Bool")
  //       ||
  //           containsOfTString(vvTypeNameLeafTree[iTree][indexName], "bool"))
  //           {
  //         lowerCorrect = 0;
  //         upperCorrect = 2;
  //         nBinCorrect = 2;
  //       } else if (containsOfTString(vvTypeNameLeafTree[iTree][indexName],
  //                                    "loat") ||
  //                  containsOfTString(vvTypeNameLeafTree[iTree][indexName],
  //                                    "ouble")) {
  //         Double_t binDensityAverage = 0;
  //         for (UInt_t iHist = 0; iHist < nHistAllFile; iHist++) {
  //           binDensityAverage += arrNBinFile[iHist] /
  //                                (arrUpperFile[iHist] - arrLowerFile[iHist])
  //                                / nHistAllFile;
  //         }
  //         Double_t minimumWholeFile = arrMinimumFile[0];
  //         Double_t maximumWholeFile = arrMaximumFile[0];
  //         for (UInt_t iHist = 1; iHist < nHistAllFile; iHist++) {
  //           if (minimumWholeFile > arrMinimumFile[iHist])
  //             minimumWholeFile = arrMinimumFile[iHist];
  //           if (maximumWholeFile < arrMaximumFile[iHist])
  //             maximumWholeFile = arrMaximumFile[iHist];
  //         }
  //         if (debug)
  //           std::cout << "minimumWholeFile, maximumWholeFile: "
  //                     << minimumWholeFile << ", " << maximumWholeFile
  //                     << std::endl;
  //         if (TMath::Abs(minimumWholeFile) > 1) {
  //           Int_t nDigitM1MinimumWholeFile =
  //               (Int_t)TMath::Floor(TMath::Log10(TMath::Abs(minimumWholeFile)));
  //           lowerCorrect =
  //               intpow(10, nDigitM1MinimumWholeFile,
  //                      TMath::Floor(intpow(10, -nDigitM1MinimumWholeFile,
  //                                          minimumWholeFile)));
  //           if (debug)
  //             std::cout << "nDigitM1MinimumWholeFile, lowerCorrect:"
  //                       << nDigitM1MinimumWholeFile << ", " << lowerCorrect
  //                       << std::endl;
  //         } else {
  //           lowerCorrect = TMath::Floor(minimumWholeFile);
  //         }
  //         if (TMath::Abs(maximumWholeFile) > 1) {
  //           Int_t nDigitM1MaximumWholeFile =
  //               (Int_t)TMath::Floor(TMath::Log10(TMath::Abs(maximumWholeFile)));
  //           upperCorrect =
  //               intpow(10, nDigitM1MaximumWholeFile,
  //                      TMath::Ceil(intpow(10, -nDigitM1MaximumWholeFile,
  //                                         maximumWholeFile)));
  //           if (debug)
  //             std::cout << "nDigitM1MaximumWholeFile, upperCorrect:"
  //                       << nDigitM1MaximumWholeFile << ", " << upperCorrect
  //                       << std::endl;
  //         } else {
  //           upperCorrect = TMath::Ceil(maximumWholeFile);
  //         }
  //         if (lowerCorrect == upperCorrect) {
  //           upperCorrect += 1;
  //           nBinCorrect = TMath::Ceil(binDensityAverage);
  //         } else {
  //           nBinCorrect = (upperCorrect - lowerCorrect) * binDensityAverage;
  //         }
  //         if (nBinCorrect > 1) {
  //           Int_t nDigitM1NBinCorrect =
  //               (Int_t)TMath::Floor(TMath::Log10(nBinCorrect));
  //           nBinCorrect =
  //               intpow((Int_t)10, nDigitM1NBinCorrect,
  //                      TMath::Nint(intpow((Int_t)10, -nDigitM1NBinCorrect,
  //                                         nBinCorrect)));
  //         }
  //         if (nBinCorrect <= 0) {
  //           nBinCorrect = 1;
  //         }
  //       } else {
  //         Long64_t minimumWholeFile = arrMinimumFile[0];
  //         Long64_t maximumWholeFile = arrMaximumFile[0];
  //         for (UInt_t iHist = 1; iHist < nHistAllFile; iHist++) {
  //           if (minimumWholeFile > arrMinimumFile[iHist]) {
  //             minimumWholeFile = arrMinimumFile[iHist];
  //           }
  //           if (maximumWholeFile < arrMaximumFile[iHist]) {
  //             maximumWholeFile = arrMaximumFile[iHist];
  //           }
  //         }
  //         lowerCorrect = minimumWholeFile;
  //         upperCorrect = maximumWholeFile + 1;
  //         nBinCorrect = maximumWholeFile - minimumWholeFile;
  //       }
  //       adjustHistSettingPerLeafTreeExtraIndex(nBinCorrect, lowerCorrect,
  //                                              upperCorrect, iTree,
  //                                              indexName);
  //       TString tstrSettingHistLeaf = (TString) "(" + nBinCorrect + "," +
  //                                     lowerCorrect + "," + upperCorrect +
  //                                     ")";
  //       vvHistsettingLeafTree[iTree].push_back(tstrSettingHistLeaf);
  //       if (debug) std::cout << vvTypeNameLeafTree[iTree][indexName];
  //       if (debug)
  //         std::cout << " (nBinCorrect, lowerCorrect, upperCorrect): "
  //                   << tstrSettingHistLeaf << std::endl;
  //       if (debug) std::cout << " Done." << std::endl;
  //       if (debug) std::cout << " Done." << std::endl;
  //     }
  //   }
  // }
  if (debug) {
    std::cout << "Drawing histograms with correct settings ..." << std::endl;
    if (isToUseCorrectedTempFile) {
      std::cout << "Not using tempfile tfCorrectedHist, \n"
                << "Hists will be scaled and added on the fly." << std::endl;
    } else {
      std::cout << "Using tfCorrectedHist to store corrected histograms."
                << std::endl;
    }
  }
  for (auto vAnalyzerLeaf : vvAnalyzerLeafTree) {
    for (LeafAnalyzerAbstract *analyzer : vAnalyzerLeaf) {
      analyzer->Summarize();
    }
  }

  std::vector<std::vector<UInt_t>> vvIHistLeafTree(nTT);
  for (UInt_t iTree = 0; iTree < nTT; iTree++) {
    std::vector<UInt_t> vIHistLeaf(vNLeavesTree[iTree], 0);
    vvIHistLeafTree.push_back(vIHistLeaf);
  }
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
    if (debug) std::cout << "tfInCurrent" << tfInCurrent << std::endl;
    if (debug) std::cout << "Mounting variables ...";
    mountVarsFromTDir(tfInCurrent);
    if (debug) std::cout << " Done." << std::endl;
    if (isToUseCorrectedTempFile) {
      tfCorrectedHist = TFile::Open(pathTFCorrectedHist, "update");
    }
    for (UInt_t iTree = 0; iTree < nTT; iTree++) {
      UInt_t nNameModifiedCurrentTree = vNLeavesTree[iTree];
      for (UInt_t indexName = 0; indexName < nNameModifiedCurrentTree;
           indexName++) {
        // TList *tlHistCorrected = new TList;
        // std::vector<TString> vNameHistCorrected;
        // vNameHistCorrected.clear();
        // UInt_t iFile=0;
        if (!vvvIsHistFileLeafTree[iTree][indexName][iFile]) {
          std::cerr << "(Tree index, Leaf index) (" << iTree << ", "
                    << indexName << ") not found in file index " << iFile
                    << " (" << vvNameModifiedLeafTree[iTree][indexName] << ")"
                    << std::endl;
          continue;
        }
        TH1 *histCorrected = nullptr;
        LeafAnalyzerAbstract *&analyzer = vvAnalyzerLeafTree[iTree][indexName];
        TString nameLeaf =
            analyzer->GetVNameLeafFile()[vvIHistLeafTree[iTree][indexName]];
        TString tstrHistSetting = analyzer->GetHistSetting();
        if (analyzer->GetVIsEmpty()[vvIHistLeafTree[iTree][indexName]]) {
          if (vvIHistLeafTree[iTree][indexName] == 0 &&
              analyzer->GetAreAllEmpty()) {
            histCorrected = analyzer->GetHistEmptyPreferred();
            if (debug)
              std::cout << "Got " << nameLeaf << " (all-empty leaf)."
                        << std::endl;
            TString nameHistCorrected =
                getNameHistOfFile(iTree, indexName, 0, "CorrectedAllEmpty");
            if (histCorrected == nullptr) {
              arrTT[iTree]->Draw(nameHistCorrected +
                                 (tstrHistSetting.Length() == 0
                                      ? ""
                                      : ">>" + tstrHistSetting));
              histCorrected = (TH1 *)gDirectory->Get(nameHistCorrected);
              histCorrected->Scale(vWeightFile[iFile]);
            }
            if (isToUseCorrectedTempFile) {
              histCorrected->SetDirectory(tfCorrectedHist);
              histCorrected->SetName(nameHistCorrected);
              histCorrected->Write(nameHistCorrected);
            } else {
              TString nameHistResult = getNameHistResult(iTree, indexName);
              histCorrected->SetName(nameHistResult);
              vvHistResultLeafTree[iTree][indexName] = histCorrected;
            }
          }
          if (debug && !analyzer->GetAreAllEmpty())
            std::cerr << "Empty histogram found for (" << vNameTT[iTree] << ", "
                      << vvNameModifiedLeafTree[iTree][indexName] << ") ("
                      << iTree << ", " << indexName << ")" << std::endl;
        } else {
          TString nameHistCorrected =
              getNameHistOfFile(iTree, indexName, iFile, "Corrected");
          // arrTT[iTree]->Draw(vvvNameFileLeafTree[iTree][indexName][iFile] +
          // ">>" +
          //                    nameHistCorrected +
          //                    vvHistsettingLeafTree[iTree][indexName]);

          arrTT[iTree]->Draw(nameLeaf + ">>" + nameHistCorrected +
                             tstrHistSetting);
          histCorrected = (TH1 *)gDirectory->Get(nameHistCorrected);
          histCorrected->Scale(vWeightFile[iFile]);
          if (isToUseCorrectedTempFile) {
            histCorrected->SetDirectory(tfCorrectedHist);
            histCorrected->SetName(nameHistCorrected);
            histCorrected->Write(nameHistCorrected);
          } else {
            if (vvHistResultLeafTree[iTree][indexName] == nullptr) {
              histCorrected->SetName(getNameHistResult(iTree, indexName));
              vvHistResultLeafTree[iTree][indexName] = histCorrected;
            } else {
              histCorrected->SetName(nameHistCorrected);
              vvHistResultLeafTree[iTree][indexName]->Add(histCorrected);
            }
          }
        }
        vvIHistLeafTree[iTree][indexName]++;
      }
    }
    if (areSomeInfilesClosed && nInfileOpen >= nInfileOpenMax) {
      closeTFilesIn(iFileOriginal);
    }
    if (isToUseCorrectedTempFile) {
      tfCorrectedHist->Close();
    }
  }
  if (nInfileOpenMax) {
    if (debug)
      std::cout << "Closing openable input ROOT files ..." << std::endl;
    closeTFilesIn(nFileTotOriginal - 1);
    if (debug) std::cout << "Done." << std::endl;
  }
  delete tfCorrectedHist;
  // tfAutogenHist->ReOpen("read");
  tfCorrectedHist = TFile::Open(pathTFCorrectedHist, "read");
  // tfCorrectedHist->ReOpen("read");

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
    TFile *tfOutHist = TFile::Open(funPathTFOutHist(iTree),
                                   toRecreateOutFile ? "recreate" : "update");
    for (UInt_t indexName = 0; indexName < vNLeavesTree[iTree]; indexName++) {
      TH1 *histResult;
      LeafAnalyzerAbstract *analyzer = vvAnalyzerLeafTree[iTree][indexName];
      TString nameHistResult = getNameHistResult(iTree, indexName);
      if (debug)
        std::cout << "Generating " << nameHistResult << " ..." << std::endl;
      if (!isToUseCorrectedTempFile) {
        histResult = vvHistResultLeafTree[iTree][indexName];
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
        TString nameHistCorrected =
            getNameHistOfFile(iTree, indexName, 0, "CorrectedAllEmpty");
        histResult = (TH1 *)tfCorrectedHist->Get(nameHistCorrected);
      } else {
        Bool_t isErrorUnexpected = false;
        TList *tlHistCorrectedToMerge = new TList;
        UInt_t nHist = 0;
        for (UInt_t iFile = 0; iFile < nFileTot; iFile++) {
          if (!vvvIsHistFileLeafTree[iTree][indexName][iFile]) {
            continue;
          }
          if (debug) std::cout << "iFile: " << iFile << " ";
          TString nameHistCorrected =
              getNameHistOfFile(iTree, indexName, iFile, "Corrected");
          TH1 *histCorrected = (TH1 *)tfCorrectedHist->Get(nameHistCorrected);
          if (histCorrected == nullptr) {
            std::cerr << "Fatal: Histogram not found: " << nameHistCorrected
                      << std::endl;
            isErrorUnexpected = true;
            break;
          }
          tlHistCorrectedToMerge->Add(histCorrected->Clone(nameHistCorrected));
          nHist++;
        }
        if (!nHist) {
          std::cerr << "Fatal: "
                    << "No histograms found for (" << iTree << ", " << indexName
                    << ") (" << vvNameModifiedLeafTree[iTree][indexName] << ")"
                    << std::endl;
          histResult = nullptr;
          isErrorUnexpected = true;
        } else {
          if (debug) std::cout << "Merging ...";
          histResult =
              (TH1 *)tlHistCorrectedToMerge->First()->Clone(nameHistResult);
          histResult->SetName(nameHistResult);
          histResult->Clear();
          histResult->Merge(tlHistCorrectedToMerge);
        }
        // if (isErrorUnexpected) continue;
        std::cout << " Done." << std::endl;
      }
      histResult->SetName(nameHistResult);
      histResult->SetDirectory(tfOutHist);
      histResult->SetTitle(vvNameModifiedLeafTree[iTree][indexName]);
      histResult->SetXTitle(analyzer->GetTitleLeaf());
      histResult->Write(nameHistResult);
    }
    // tfOutHist->Write();
    tfOutHist->Close();
  }
  tfCorrectedHist->Close();
  // tfAutogenHist->Close();
}

#endif
