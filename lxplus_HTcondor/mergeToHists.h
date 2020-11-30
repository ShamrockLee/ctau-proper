#ifndef MERGETOHIST_H
#define MERGETOHIST_H

#include <Rtypes.h>
#include <TString.h>

#include <functional>
#include <vector>

class LeafAnalyzerAbstract {
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
  static TString GetNameLeafModified(TString nameLeaf);
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

void mergeToHists(const std::vector<TString> vNameTT,
                 std::vector<UInt_t> vNumberFile,
                 std::vector<Double_t> vCrossSection,
                 std::function<TString(UInt_t)> funPathTFIn, TString dirTFTemp,
                 std::function<TString(TString)> funNameTFTemp,
                 TString dirTFOut, std::function<TString(TString)> funNameTFOut,
                 TString seperatorPath = "/",
                 std::function<void(Int_t& NBinCorrect, Double_t& lowerCorrect,
                                    Double_t& upperCorrect, TString nameTT,
                                    TString nameLeafModified,
                                    TString typeNameLeaf, TString titleLeaf)>
                     adjustHistSettingPerLeafTreeExtra = nullptr,
                 std::function<TString (TString)> getNameLeafModified = nullptr,
                 std::function<LeafAnalyzerAbstract *()> supplyLeafAnalyzer = nullptr,
                 const Bool_t toRecreateOutFile = true,
                 const Bool_t debug = false, const Bool_t allowMissing = false,
                 UInt_t nInfileOpenMax = 30, UInt_t nLeafWithoutCorrectTempFileMin = 100);

void mergeToHists(const std::vector<TString> vNameTT,
                 std::vector<UInt_t> vNumberFile,
                 std::vector<Double_t> vCrossSection, TString patternPathTFIn,
                 TString dirTFTemp, TString patternNameTFTemp, TString dirTFOut,
                 TString patternNameTFOut, TString seperatorPath = "/",
                 std::function<void(Int_t& NBinCorrect, Double_t& lowerCorrect,
                                    Double_t& upperCorrect, TString nameTT,
                                    TString nameLeafModified,
                                    TString typeNameLeaf, TString titleLeaf)>
                     adjustHistSettingPerLeafTreeExtra = nullptr,
                 std::function<TString (TString)> getNameLeafModified = nullptr,
                 std::function<LeafAnalyzerAbstract *()> supplyLeafAnalyzer = nullptr,
                 const Bool_t toRecreateOutFile = true,
                 const Bool_t debug = false,
                 const Bool_t allowMissing = false,
                 UInt_t nInfileOpenMax = 30, UInt_t nLeafWithoutCorrectTempFileMin = 100) {
  std::function<TString(UInt_t)> funPathTFIn =
      [patternPathTFIn](UInt_t iFileOriginal) -> TString {
    return Form(patternPathTFIn.Data(), iFileOriginal);
  };
  std::function<TString(TString)> funNameTFTemp =
      [patternNameTFTemp](TString keyword) -> TString {
    return Form(patternNameTFTemp.Data(), keyword.Data());
  };
  std::function<TString(TString)> funNameTFOut =
      [patternNameTFOut](TString nameTT) -> TString {
    return Form(patternNameTFOut.Data(), nameTT.Data());
  };
  mergeToHists(vNameTT, vNumberFile, vCrossSection,
              patternPathTFIn == "" ? nullptr : funPathTFIn, dirTFTemp,
              patternNameTFTemp == "" ? nullptr : funNameTFTemp, dirTFOut,
              patternNameTFOut == "" ? nullptr : funNameTFOut, seperatorPath,
              adjustHistSettingPerLeafTreeExtra, getNameLeafModified, supplyLeafAnalyzer,
              toRecreateOutFile, debug,
              allowMissing,
              nInfileOpenMax, nLeafWithoutCorrectTempFileMin);
}
#endif