#include <Rtypes.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TString.h>
#include <TTree.h>

#include <cstdlib>
#include <functional>

#ifndef HISTMERGER_H
#define HISTMERGER_H
class HistMerger {
 public:
  class LeafAnalyzerAbstract;
  class LeafAnalyzerDefault;

  void Run();

  Bool_t debug;
  std::vector<TString> vNameTT;
  std::vector<UInt_t> vNumberFile;
  std::vector<Double_t> vCrossSection;
  std::function<TString(UInt_t iFile)> funPathTFIn;
  TString dirTFTemp;
  std::function<TString(TString keywordTFTemp)> funNameTFTemp;
  TString dirTFOut;
  std::function<TString(TString nameTT)> funNameTFOut;
  TString seperatorPath;
  std::function<void(Int_t& NBinCorrect, Double_t& lowerCorrect,
                     Double_t& upperCorrect, TString nameTT,
                     TString nameLeafModified, TString typeNameLeaf,
                     TString titleLeaf)>
      adjustHistSettingPerLeafTreeExtra = nullptr;
  std::function<TString(TString nameLeaf)> getNameLeafModified = nullptr;
  std::function<LeafAnalyzerAbstract*()> supplyLeafAnalyzer = nullptr;
  std::function<Bool_t(TString nameTT, TString nameLeafModified)>
      getIsToVetoLeaf;
  Bool_t toRecreateOutFile;
  Bool_t allowMissing;
  UInt_t nInfileOpenMax;
  UInt_t nLeavesToUseCorrectedTempFileMin;

  void SetNameTFTemp(std::function<TString(TString keyword)> funNameTFTemp) {
    this->funNameTFTemp = funNameTFTemp;
  }
  void SetNameTFTemp(TString fmt) {
    funNameTFTemp = [fmt](TString keyword) -> TString {
      return TString::Format(fmt, keyword.Data());
    };
  }
  void SetNameTFOut(std::function<TString(TString nameTT)> funNameTFOut) {
    this->funNameTFOut = funNameTFOut;
  }
  void SetNameTFOut(TString fmt) {
    funNameTFOut = [fmt](TString nameTT) -> TString {
      return TString::Format(fmt, nameTT.Data());
    };
  }
  template<class LeafAnalyzer>
  void SetLeafAnalyzer() {
    supplyLeafAnalyzer = [](){return (LeafAnalyzerAbstract *)(new LeafAnalyzer());};
    getNameLeafModified = [](TString nameLeaf){return (TString) LeafAnalyzer::GetNameLeafModified(nameLeaf);};
  }

  HistMerger();

 protected:
  /// The weight of each dataset
  std::vector<Double_t> vWeightDataset;
  /// The weight of each openable file
  std::vector<Double_t> vWeightFile;
  /// The existence in each file of each leaf in each tree
  std::vector<std::vector<std::vector<Bool_t>>> vvvIsHistFileLeafTree;
  std::vector<std::vector<TH1*>> vvHistResultLeafTree;
  std::vector<UInt_t> vNLeavesTree;
  /// The modified leafname of each leaf in each tree
  std::vector<std::vector<TString>> vvNameModifiedLeafTree;
  /// The leaf analyzer of each leaf in each tree;
  std::vector<std::vector<LeafAnalyzerAbstract*>> vvAnalyzerLeafTree;

  void InitializeHidden();
  void InitializeWhenRun();
};

HistMerger::HistMerger() {
  this->funPathTFIn = nullptr;
  this->dirTFTemp = "";
  this->funNameTFTemp = nullptr;
  this->dirTFOut = "";
  this->funNameTFOut = nullptr;
  this->seperatorPath = "/";
  this->getNameLeafModified = nullptr;
  this->supplyLeafAnalyzer = nullptr;
  this->getIsToVetoLeaf = nullptr;
  this->toRecreateOutFile = true;
  this->debug = false;
  this->allowMissing = false;
  this->nInfileOpenMax = 30;
  this->nLeavesToUseCorrectedTempFileMin = 100;

  InitializeHidden();
}

class HistMerger::LeafAnalyzerAbstract {
 public:
  void SetDebug(Bool_t debug){};
  void SetNameTT(TString nameTT){};
  void SetFunAdjustHistSettingExtra(
      std::function<void(Int_t& NBinCorrect, Double_t& lowerCorrect,
                         Double_t& upperCorrect, TString nameTT,
                         TString nameLeafModified, TString typeNameLeaf,
                         TString titleLeaf)>
          adjustHistSettingPerLeafTreeExtra){};
  void AnalyzeLeafBasic(TLeaf* leaf);
  // const std::function<TString(TString nameLeaf)> getNameLeafModified;
  TString GetNameLeafModified();
  TString GetTitleLeaf();
  TString GetTypeNameLeaf();
  void AnalyzeLeaf(TLeaf* leaf);
  std::vector<TString>& GetVNameLeafFile();
  std::vector<Bool_t>& GetVIsEmpty();
  void Summarize();
  Bool_t GetAreAllEmpty();
  Int_t GetNBinsCorrect();
  Double_t GetLowerCorrect();
  Double_t GetUpperCorrect();
  TString GetHistSetting();
  TH1* GetHistEmptyPreferred() { return nullptr; };
  void Finalize(){};
};

#endif