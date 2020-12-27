#include <Rtypes.h>
#include <TCut.h>
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

  virtual void Run();

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
      adjustHistSettingPerLeafTreeExtra;
  std::function<TString(TString nameLeaf)> getNameLeafModified;
  std::function<LeafAnalyzerAbstract*()> supplyLeafAnalyzer;
  std::function<Bool_t(TString nameTT, TString nameLeafModified)>
      funIsToVetoLeaf;
  Bool_t toRecreateOutFile;
  Bool_t allowMissing;
  UInt_t nInfileOpenMax;
  UInt_t nLeavesToUseCorrectedTempFileMin;
  std::vector<std::vector<LeafAnalyzerAbstract*>> vvAnalyzerLeafTreeCustom;
  std::vector<LeafAnalyzerAbstract*> vAnalyzerCustomByName;
  std::function<void(
      TTree* tree,
      const std::vector<LeafAnalyzerAbstract*> vAnalyzerLeafCustom,
      const std::function<void(LeafAnalyzerAbstract* analyzerNew)> pushbackNewAnalyzer)>
      pushCustomAnalyzersWhenRun;
  std::function<TString(TLeaf* leaf)> funTitleLeaf;

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
  template <class LeafAnalyzer>
  void SetLeafAnalyzer() {
    supplyLeafAnalyzer = []() {
      return (LeafAnalyzerAbstract*)(new LeafAnalyzer());
    };
    getNameLeafModified = [](TString nameLeaf) {
      return (TString)LeafAnalyzer::GetNameLeafModified(nameLeaf);
    };
  }

  HistMerger();

 protected:
  UInt_t nTT;
  /// The weight of each dataset
  std::vector<Double_t> vWeightDataset;
  /// The weight of each openable file
  std::vector<Double_t> vWeightFile;
  /// The existence in each file of each leaf in each tree
  std::vector<std::vector<std::vector<Bool_t>>> vvvIsHistFileLeafTree;
  std::vector<std::vector<std::vector<Bool_t>>> vvvIsHistFileLeafTreeCustom;
  std::vector<std::vector<TH1*>> vvHistResultLeafTree;
  std::vector<TFile*> vTFOut;
  /// The number of leaves in each tree
  std::vector<UInt_t> vNLeavesTree;
  /// The number of custom analyzers in each tree
  std::vector<UInt_t> vNAnalyzersTreeCustom;
  /// The modified leafname of each leaf in each tree
  std::vector<std::vector<TString>> vvNameModifiedLeafTree;
  /// The leaf analyzer of each leaf in each tree;
  std::vector<std::vector<LeafAnalyzerAbstract*>> vvAnalyzerLeafTree;

  virtual void InitializeHidden();
  virtual void InitializeWhenRun();

  virtual inline void ProcessAnalyzerNew(LeafAnalyzerAbstract *analyzer, TString nameTT);
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
  this->funIsToVetoLeaf = nullptr;
  this->toRecreateOutFile = true;
  this->debug = false;
  this->allowMissing = false;
  this->nInfileOpenMax = 30;
  this->nLeavesToUseCorrectedTempFileMin = 100;
  this->vvAnalyzerLeafTreeCustom.clear();
  this->vAnalyzerCustomByName.clear();
  this->pushCustomAnalyzersWhenRun = nullptr;
  this->funTitleLeaf = nullptr;
  InitializeHidden();
}

class HistMerger::LeafAnalyzerAbstract {
 public:
  virtual void SetDebug(Bool_t debug){};
  virtual void SetNameTT(TString nameTT){};
  virtual TString GetNameTT() const { return ""; };
  virtual void SetExpressionCustom(
      TString name, TString typeName, TString title,
      TString expressionBeforeSetting, TCut selection = "",
      Option_t* option = "", Long64_t nentries = TTree::kMaxEntries,
      Long64_t firstentry = 0) {}
  void SetExpressionBeforeSetting(TString expressionBeforeSetting) {};
  virtual void SetSelection(TCut selection) {};
  virtual void SetOptionDraw(Option_t *option) {};
  virtual void SetNEntriesToDrawMax(Long64_t nentries) {};
  virtual void SetFirstEntryToDraw(Long64_t firstentry) {};
  virtual TString GetExpressionBeforeSetting() const {return "";};
  virtual void SetHasTarget(
      std::vector<TString> vDeps = {},
      std::function<Bool_t(TTree* tree)> funHasTargetExtra = nullptr) {}
  virtual void SetAllowNeverAnalyzed(Bool_t allowNeverAnalyzed) {}
  virtual Bool_t GetAllowNeverAnalyzed() const { return false; }
  virtual Bool_t GetIsEverAnalyzed() const { return true; }
  virtual void AssignHistSetting(Int_t nBinsCorrect, Double_t lowerCorrect,
                                 Double_t upperCorrect,
                                 TString tstrHistSetting = "") {}
  virtual void SetFunAssignHistSettingExtra(
      std::function<Bool_t(Int_t& NBinCorrect, Double_t& lowerCorrect,
                           Double_t& upperCorrect, TString nameTT,
                           TString nameLeafModified, TString typeNameLeaf,
                           TString titleLeaf)>
          assignHistSettingPerLeafTreeExtra) {}
  virtual void SetDontCheckEmptyness(Bool_t dontCheckEmptyness) {}
  virtual Bool_t GetDontCheckEmptyness() const { return false; }
  virtual void SetFunAdjustHistSettingExtra(
      std::function<void(Int_t& NBinCorrect, Double_t& lowerCorrect,
                         Double_t& upperCorrect, TString nameTT,
                         TString nameLeafModified, TString typeNameLeaf,
                         TString titleLeaf)>
          adjustHistSettingPerLeafTreeExtra) {}
  virtual void AnalyzeLeafBasic(TLeaf* leaf) = 0;
  // const std::function<TString(TString nameLeaf)> getNameLeafModified;
  virtual TString GetNameLeafModified() const = 0;
  virtual void SetFunTitleLeaf(std::function<TString(TLeaf* leaf)> funTitleLeaf) {};
  virtual TString GetTitleLeaf() const = 0;
  virtual TString GetTypeNameLeaf() const = 0;
  virtual Bool_t GetHasTarget(TTree* tree) const { return true; }
  virtual void AnalyzeLeaf(TLeaf* leaf) = 0;
  virtual void EvaluateAndAnalyze(TTree *tree) {};
  virtual std::vector<TString>& GetVNameLeafFile() = 0;
  virtual std::vector<Bool_t>& GetVIsEmpty() = 0;
  virtual void Summarize() = 0;
  virtual Bool_t GetAreAllEmpty() const = 0;
  virtual Int_t GetNBinsCorrect() const = 0;
  virtual Double_t GetLowerCorrect() const = 0;
  virtual Double_t GetUpperCorrect() const = 0;
  virtual TString GetHistSetting() const = 0;
  virtual TH1* GetHistEmptyPreferred() const { return nullptr; };
  virtual TH1* DrawHistCorrected(TString nameHist, TTree* tree,
                                 Int_t iHist = -1,
                                 Bool_t isToClone = false) = 0;
  virtual void Finalize(){};
};

#endif