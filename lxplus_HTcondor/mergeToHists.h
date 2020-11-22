#ifndef MERGETOHIST_H
#define MERGETOHIST_H

#include <Rtypes.h>
#include <TString.h>

#include <functional>
#include <vector>
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
                 const Bool_t toRecreateOutFile = true,
                 const Bool_t debug = false, const Bool_t allowMissing = false);

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
                 const Bool_t toRecreateOutFile = true,
                 const Bool_t debug = false,
                 const Bool_t allowMissing = false) {
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
              adjustHistSettingPerLeafTreeExtra, toRecreateOutFile, debug,
              allowMissing);
}
#endif