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
  std::vector<TFile *> vTFOut;
  /// The number of leaves in each tree
  std::vector<UInt_t> vNLeavesTree;
  /// The modified leafname of each leaf in each tree
  std::vector<std::vector<TString>> vvNameModifiedLeafTree;
  /// The leaf analyzer of each leaf in each tree;
  std::vector<std::vector<LeafAnalyzerAbstract*>> vvAnalyzerLeafTree;

  virtual void InitializeHidden();
  virtual void InitializeWhenRun();
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
  virtual void SetDebug(Bool_t debug){};
  virtual void SetNameTT(TString nameTT){};
  virtual TString GetNameTT(){return "";};
  virtual void SetExpressionBeforeSetting(TString expressionBeforeSetting,
                                          TString name,
                                          TString typeName="",
                                          TString title="",
                                          std::function<Bool_t (TTree *)> getHasTarget=nullptr,
                                          Int_t nBinsCorrect=-1,
                                          Double_t lowerCorrect=0,
                                          Double_t upperCorrect=0,
                                          Bool_t dontCheckEmptiness=false){};
  virtual void SetFunAssignHistSettingExtra(
      std::function<Bool_t(Int_t& NBinCorrect, Double_t& lowerCorrect,
                         Double_t& upperCorrect, TString nameTT,
                         TString nameLeafModified, TString typeNameLeaf,
                         TString titleLeaf)>
          assignHistSettingPerLeafTreeExtra){}
  virtual void SetDontCheckEmptyness(Bool_t dontCheckEmptyness){}
  virtual Bool_t GetDontCheckEmptyness(){return false;}
  virtual void SetFunAdjustHistSettingExtra(
      std::function<void(Int_t& NBinCorrect, Double_t& lowerCorrect,
                         Double_t& upperCorrect, TString nameTT,
                         TString nameLeafModified, TString typeNameLeaf,
                         TString titleLeaf)>
          adjustHistSettingPerLeafTreeExtra){}
  virtual void AnalyzeLeafBasic(TLeaf* leaf) = 0;
  // const std::function<TString(TString nameLeaf)> getNameLeafModified;
  virtual TString GetNameLeafModified() = 0;
  virtual TString GetTitleLeaf() = 0;
  virtual TString GetTypeNameLeaf() = 0;
  virtual void AnalyzeLeaf(TLeaf* leaf) = 0;
  virtual std::vector<TString>& GetVNameLeafFile() = 0;
  virtual std::vector<Bool_t>& GetVIsEmpty() = 0;
  virtual void Summarize() = 0;
  virtual Bool_t GetAreAllEmpty() = 0;
  virtual Int_t GetNBinsCorrect() = 0;
  virtual Double_t GetLowerCorrect() = 0;
  virtual Double_t GetUpperCorrect() = 0;
  virtual TString GetHistSetting() = 0;
  virtual TH1* GetHistEmptyPreferred() { return nullptr; };
  virtual void Finalize(){};
};

#endif