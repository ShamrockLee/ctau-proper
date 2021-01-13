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

  /// Method that runs the merging
  ///
  /// Execute after the settings are done.
  virtual void Run();

    /// Set funNameTFTemp
  void SetNameTFTemp(std::function<TString(TString keyword)> funNameTFTemp) {
    this->funNameTFTemp = funNameTFTemp;
  }
  /// Set funNameTFTemp using `TString::Format(fmt, keyword.Data())`
  void SetNameTFTemp(TString fmt) {
    funNameTFTemp = [fmt](TString keyword) -> TString {
      return TString::Format(fmt, keyword.Data());
    };
  }
  /// Set funNameTFOut
  void SetNameTFOut(std::function<TString(TString nameTT)> funNameTFOut) {
    this->funNameTFOut = funNameTFOut;
  }
  /// Set funNameTFOut using `TString::Format(fmt, nameTT.Data())`
  void SetNameTFOut(TString fmt) {
    funNameTFOut = [fmt](TString nameTT) -> TString {
      return TString::Format(fmt, nameTT.Data());
    };
  }
  /// Set the LeafAnalyzer class to use
  ///
  /// If not set,
  /// SetLeafAnalyzer<LeafAnalyzerDefault>() will be run.
  template <class LeafAnalyzer>
  void SetLeafAnalyzer() {
    supplyLeafAnalyzer = []() {
      return (LeafAnalyzerAbstract*)(new LeafAnalyzer());
    };
    funNameLeafModified = [](TString nameLeaf) {
      return (TString)LeafAnalyzer::GetNameLeafModified(nameLeaf);
    };
  }
  /// The debug swich for debug messages and extra examinations
  ///
  /// Default to false
  Bool_t debug;
  /// The name of the trees
  ///
  /// Set directly
  ///
  /// REQUIRED
  std::vector<TString> vNameTT;
  /// The number of the files in each dataset
  ///
  /// Set directly
  ///
  /// REQUIRED
  std::vector<UInt_t> vNumberFile;  //<
  /// The cross section of each dataset
  ///
  /// Set directly
  ///
  /// REQUIRED
  std::vector<Double_t> vCrossSection;
  /// `std::function` object that decide the path of the inpum files according
  /// to their index.
  ///
  /// Set directly
  ///
  /// REQUIRED
  std::function<TString(UInt_t iFile)> funPathTFIn;
  /// Directory for the temp files
  ///
  /// Default to `"."`
  TString dirTFTemp;
  /// std::function object that decide the name of the temp files according to
  /// the keyword
  ///
  /// Set directly or using
  /// ~~~ {.C}
  /// void SetNameTFTemp(std::function<TString(TString keyword)> funNameTFTemp)`
  /// void SetNameTFTemp(TString fmt)
  /// ~~~
  ///
  /// Default to nullptr.
  ///
  /// When nullptr is encountered, the following will be run:
  /// ~~~ {.C}
  /// char *charsTemp = nullptr;
  /// strcat(charsTemp, "XXXXXX");
  /// mkstemp(charsTemp);
  /// SetNameTFTemp((TString) "%s" + charsTemp + ".root");
  /// ~~~
  std::function<TString(TString keywordTFTemp)> funNameTFTemp;
  /// Directory for the output files
  ///
  /// Default to `"."`
  TString dirTFOut;
  /// `std::function` object that decide the name of the output file according
  /// to the tree name and the leaf name
  ///
  /// Set directly or using
  /// ~~~ {.C}
  /// void SetNameTFOut(std::function<TString(TString keyword)> funNameTFOut)`
  /// void SetNameTFOut(TString fmt)
  /// ~~~
  ///
  /// Default to `nullptr`.
  ///
  /// When nullptr is encountered, the following will be run:
  /// ~~~ {.C}
  /// funNameTFOut = [fmt](TString nameTT) -> TString {
  ///   return TString::Format(fmt, nameTT.Data());
  /// };
  /// ~~~
  std::function<TString(TString nameTT)> funNameTFOut;
  /// Path seperator
  ///
  /// Default to "/" (the POSIX style)
  TString seperatorPath;
  std::function<Bool_t(Int_t &NBinCorrect, Double_t &lowerCorrect,
                       Double_t &upperCorrect, TString nameTT,
                       TString nameLeafModified, TString typeNameLeaf,
                       TString titleLeaf)>
      assignHistSettingPerLeafTreeExtra;
  /// std::function object that
  /// adjust the automatically decided hist settings
  /// (nBinsCorrect, lowerCorrect, upperCorrect)
  ///
  /// Default to nullptr
  ///
  /// Example:
  /// ~~~ {.C}
  /// adjustHistSettingPerLeafTreeExtra = [](
  ///                Int_t& NBinCorrect, Double_t& lowerCorrect,
  ///                Double_t& upperCorrect, TString nameTT,
  ///                TString nameLeafModified, TString typeNameLeaf,
  ///                TString titleLeaf)->TString{
  ///   if (typeNameLeaf.Contains("loat") || typeNameLeaf.Contains("ouble")) {
  ///     if (nameLeafModified.EndsWith("_phi")) {
  ///       lowerCorrect = -TMath::Pi;
  ///       upperCorrect = TMath::Pi;
  ///     }
  ///   }
  /// }
  /// ~~~
  std::function<void(Int_t& NBinCorrect, Double_t& lowerCorrect,
                     Double_t& upperCorrect, TString nameTT,
                     TString nameLeafModified, TString typeNameLeaf,
                     TString titleLeaf)>
      adjustHistSettingPerLeafTreeExtra;
  /// std::function object that modifies the leaf name.
  ///
  /// Set directly or by SetLeafAnalyzer<LeafAnalyzer>()
  ///
  /// Default to `static_cast<TString(*)(TString)>(LeafAnalyzerDefault::GetNameLeafModified)`,
  /// which is equivalant to
  /// ~~~ {.C}
  /// [](TString nameLeaf)->TString{
  ///   indexLeafBracket = nameLeaf.First('[');
  ///   return indexLeafBracket == -1 ? nameLeaf : nameLeaf(0, indexFirstBracket - 1);
  /// }
  /// ~~~
  std::function<Bool_t(TString nameTT, TString nameLeafModified)>
      funIsToVetoLeaf;
  /// Whether to recreate the output files or to append.
  ///
  /// Default to `true` (the former).
  Bool_t toRecreateOutFile;
  /// Whether to allow missing input files.
  ///
  /// Set directly
  ///
  /// Default to `false`, for errors shouldn't be passed silently.
  /// However, it is useful to switch to `true` because
  /// it is hard to guarentee all the jobs on the GRID are passed,
  /// and missing one or two files does not nessesarily influence
  /// the accuracy of the result.
  Bool_t allowMissingInputFiles;
  /// The maximum number of input files to stay open when run.
  /// When the number of opening input files reaches this number,
  /// The then-opening files will be closed
  /// To prevent too many input files opened simultaniously.
  ///
  /// Set directly
  ///
  /// Default to 30
  ///
  /// Set to 0 to opt out all the input-file-closing operations.
  UInt_t nInfileOpenMax;
  /// The minimum number of leaves to use temporary files to store the
  /// histograms waiting for merging.
  /// This is to prevent too many output files being opened simultaniously.
  ///
  /// Set directly
  ///
  /// Default to 100
  ///
  /// Set to 0 to opt out the use of temp files
  UInt_t nLeavesToUseCorrectedTempFileMin;
  /// Custom analyzers placed inside nested vectors
  /// with vvAnalyzerLeafTreeCustom[iTree] representing
  /// the vector of custom analyzers to apply on tree named vNameTT[iTree]
  ///
  /// A custom analyzers are set to store
  /// the expression, the selection, the option, the maximum entries, and the first entry
  /// to be use as tree->Draw(...) arguments,
  /// allowing custom expressions to be drawn automatically.
  ///
  /// Set directly
  ///
  /// Default to {}
  ///
  /// If the size does not match vNameTT when Run(), it will be resize(nTT, {})
  std::vector<std::vector<LeafAnalyzerAbstract*>> vvAnalyzerLeafTreeCustom;
  /// Custom analyzes placed inside a vector
  /// to be added to vvAnalyzerLeafTreeCustom when Run()
  /// according to the tree name set inside the analyzer
  std::vector<LeafAnalyzerAbstract*> vAnalyzerCustomByName;
  std::function<void(
      TTree* tree, const std::vector<LeafAnalyzerAbstract*> vAnalyzerLeafCustom,
      const std::function<void(LeafAnalyzerAbstract* analyzerNew)>
          pushbackNewAnalyzer)>
      pushCustomAnalyzersWhenRun;
  /// `std::function` object to get a title from a leaf
  ///
  /// Set directly
  ///
  /// Default to
  /// ~~~ {.C}
  /// [](TLeaf *leaf)->TString{ return leaf->GetTitle() == "" ? leaf->GetBranch()->GetTitle() : leaf->GetTitle(); }
  /// ~~~
  std::function<TString(TLeaf* leaf)> funTitleLeaf;
  

  HistMerger();
  virtual ~HistMerger(){};

protected:
  /// std::function that opts out certain leaf according to the tree name and the leaf name.
  /// Not applied to the custom analyzers
  ///
  /// Set through SetLeafAnalyzer<LeafAnalyzer>()
  ///
  /// Default is set through SetLeafAnalyzer<LeafAnalyzerDefault>() to
  /// ~~~ {.C}
  /// [](TString nameLeaf) {
  ///   return (TString)LeafAnalyzer::GetNameLeafModified(nameLeaf);
  /// }
  /// ~~~
  std::function<TString(TString nameLeaf)> funNameLeafModified;
  /// std::function that supply the pointer to a new LeafAnalyzerAbstract instance
  ///
  /// Set through SetLeafAnalyzer<LeafAnalyzer>()
  ///
  /// Default is set through SetLeafAnalyzer<LeafAnalyzerDefault>() to
  /// ~~~ {.C}
  /// []()->LeafAnalyzerAbstract{ return (LeafAnalyzerAbstract *)(new LeafAnalyzerDefault); }
  /// ~~~
  std::function<LeafAnalyzerAbstract*()> supplyLeafAnalyzer;
 protected:
  /// Number of trees
  ///
  /// Will be assign to be vNameTT.size() when Run()
  UInt_t nTT;
  /// The weight of each dataset
  std::vector<Double_t> vWeightDataset;
  /// The weight of each openable file
  std::vector<Double_t> vWeightFile;
  /// The existence of leaves in each file of each leaf in each tree
  std::vector<std::vector<std::vector<Bool_t>>> vvvIsHistFileLeafTree;
  /// The availability of custom analyzers in each file of each leaf in each tree
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

  virtual inline void ProcessAnalyzerNew(LeafAnalyzerAbstract* analyzer,
                                         TString nameTT);
};

HistMerger::HistMerger() {
  this->funPathTFIn = nullptr;
  this->dirTFTemp = "";
  this->funNameTFTemp = nullptr;
  this->dirTFOut = "";
  this->funNameTFOut = nullptr;
  this->seperatorPath = "/";
  this->funIsToVetoLeaf = nullptr;
  this->toRecreateOutFile = true;
  this->debug = false;
  this->allowMissingInputFiles = false;
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
  virtual void SetExpressionCustom(TString name, TString typeName,
                                   TString title,
                                   TString expressionBeforeSetting,
                                   TCut selection = "", Option_t* option = "",
                                   Long64_t nentries = TTree::kMaxEntries,
                                   Long64_t firstentry = 0) {}
  void SetExpressionBeforeSetting(TString expressionBeforeSetting){};
  virtual void SetSelection(TCut selection){};
  virtual void SetOptionDraw(Option_t* option){};
  virtual void SetNEntriesToDrawMax(Long64_t nentries){};
  virtual void SetFirstEntryToDraw(Long64_t firstentry){};
  virtual TString GetExpressionBeforeSetting() const { return ""; };
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
  virtual void SetFunTitleLeaf(
      std::function<TString(TLeaf* leaf)> funTitleLeaf){};
  virtual TString GetTitleLeaf() const = 0;
  virtual TString GetTypeNameLeaf() const = 0;
  virtual Bool_t GetHasTarget(TTree* tree) const { return true; }
  virtual void AnalyzeLeaf(TLeaf* leaf) = 0;
  virtual void EvaluateAndAnalyze(TTree* tree){};
  virtual std::vector<TString> GetVNameLeafFile() const = 0;
  virtual std::vector<Bool_t> GetVIsEmpty() const = 0;
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
  virtual ~LeafAnalyzerAbstract(){};
};

#endif
