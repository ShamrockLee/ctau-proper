#include <TROOT.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TString.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>

#include <functional>
#include <iostream>
#include <cstdio>

void plotAll(TString pathFileIn, TString dirOut, const Bool_t plotSubdir = true, const Bool_t normalize = false, const Bool_t logy = false, const Option_t *optionDraw = "", const std::function<TH1*(TH1*)> funAdjustHist = nullptr);

void plotAll(TDirectory *tdirIn, TString dirOut, const Bool_t plotSubdir, const Bool_t normalize, const Bool_t logy, const Option_t *optionDraw, const std::function<TH1*(TH1*)> funAdjustHist) {
  TString seperatorPath = "/";
  gSystem->mkdir(dirOut);
  if (dirOut.Length() > 1 && dirOut.EndsWith(seperatorPath)) {
    dirOut.Resize(dirOut.Length()-1);
  }
  if (!gPad) {
    gROOT->MakeDefCanvas();
  }
  gPad->SetLogy(logy);
  for (auto&& keyRaw: *tdirIn->GetListOfKeys()) {
    TKey *key = (TKey *) keyRaw;
    if (key == nullptr) continue;
    TString name = key->GetName();
    TObject *obj = key->ReadObj();
    if (obj == nullptr) continue;
    if (plotSubdir && obj->IsA()->InheritsFrom("TDirectory")) {
      plotAll(static_cast<TDirectory*>(obj), dirOut + seperatorPath + name, true, normalize, logy, optionDraw, funAdjustHist);
    }
    if (obj->IsA()->InheritsFrom("TH1")
        && ! obj->IsA()->InheritsFrom("TH2")
        && ! obj->IsA()->InheritsFrom("TH3")
      ) {
      TH1 *hist = (TH1 *) key->ReadObj();
      if (normalize) {
        hist->Scale(1.0 / hist->Integral());
      }
      if (funAdjustHist != nullptr) hist = funAdjustHist(hist);
      const char *optionDrawNew = optionDraw;
      if (hist->GetNbinsX() <= 13 && (! TString(optionDraw).Contains("Text"))) {
        optionDrawNew = (TString(optionDraw) + "Text").Data();
      }
      hist->Draw(optionDrawNew);
      gStyle->SetOptStat(111111);
      gPad->Print(dirOut + seperatorPath + name + ".svg");
    }
  }
  tdirIn->Close();
}

void plotAll(const TString pathFileIn, const TString dirOut, const Bool_t plotSubdir, const Bool_t normalize, const Bool_t logy, const Option_t *optionDraw, const std::function<TH1*(TH1*)> funAdjustHist) {
  // TString basenameFileIn = gSystem->GetFromPipe(
  //     Form("file=%s; test=${file##*/}; echo \"${test%%.root}\"",
  //     pathFileIn.Data()));
  gSystem->mkdir(dirOut, true);
  TFile *fileIn = TFile::Open(pathFileIn);
  plotAll(fileIn, dirOut, plotSubdir, normalize, logy, optionDraw, funAdjustHist);
}

void plotAll(const char* pathFileIn, const char* dirout, const Bool_t plotSubdir = true, const Bool_t normalize = false, const Bool_t logy = false, const Option_t *optionDraw = "", const std::function<TH1*(TH1*)> funAdjustHist = nullptr) {
  plotAll((TString) pathFileIn, (TString) dirout, plotSubdir, normalize, logy, optionDraw, funAdjustHist);
}

TH1* RebinTo100(TH1 * hist) {
  TH1 *resultHist = hist;
  const Int_t nBins = hist->GetNbinsX();
  if (nBins) {
    const Double_t binWidth = hist->GetBinWidth(1);
    if (TMath::Abs(binWidth - TMath::Nint(binWidth)) < 0.001 && TMath::Nint(binWidth) > 0.5 && nBins >= 200) {
      resultHist = hist->Rebin(nBins / 100);
    }
  }
  return resultHist;
}

#if true

#include "nlohmann_json.hpp"

TH1* adjustWithJSONTitle(TH1 *hist) {
  const char* strSetting = hist->GetTitle();
  std::cerr << "name: " << hist->GetName() << ", ";
  std::cerr << "strSetting: "  << strSetting << std::endl;
  const nlohmann::json jSetting = nlohmann::json::parse(strSetting);
  const double alignment = jSetting["alignment"].get<double>();
  const int binDensityOrder = jSetting["binDensityOrder"].get<int>();
  const bool isLowerAssigned = jSetting["isLowerAssigned"].get<bool>();
  const int lowerLimitBins = jSetting["lowerLimitBins"].get<int>();
  const bool isUpperAssigned = jSetting["isUpperAssigned"].get<bool>();
  const int upperLimitBins = jSetting["upperLimitBins"].get<int>();
  if (!isLowerAssigned || !isUpperAssigned) {
    // Round to 12 to 20 instead of 100
    const int biasFactor = 2.;
    int lowerLimitBinsCorrect = lowerLimitBins;
    int upperLimitBinsCorrect = upperLimitBins;
    const int lowerLimitBinsNonzero = hist->FindFirstBinAbove() + lowerLimitBins - 1;
    const int upperLimitBinsNonzero = hist->FindLastBinAbove() + lowerLimitBins;
    if (upperLimitBinsNonzero == lowerLimitBins - 1) {
      if (isLowerAssigned) {
        upperLimitBinsCorrect = lowerLimitBins + 1;
      } else if (isUpperAssigned) {
        lowerLimitBinsCorrect = upperLimitBins - 1;
      } else {
        lowerLimitBinsCorrect = -(alignment > 0);
        upperLimitBinsCorrect = 1 - (alignment > 0);
      }
    }
    if (isLowerAssigned) {
      if (upperLimitBinsNonzero > 0) {
        const int limitOrder = TMath::Ceil(TMath::Log10(upperLimitBinsNonzero / biasFactor + (upperLimitBinsNonzero % biasFactor != 0)));
        const int limitBase = static_cast<int>(TMath::Power(10, limitOrder));
        upperLimitBinsCorrect = limitBase * (upperLimitBinsNonzero / limitBase + (upperLimitBinsNonzero % limitBase != 0));
      } else {
        upperLimitBinsCorrect = (alignment < 0);
      }
    } else if (isUpperAssigned) {
      if (lowerLimitBinsNonzero < 0) {
        const int limitOrder = TMath::Ceil(TMath::Log10(-(lowerLimitBinsNonzero / biasFactor + (lowerLimitBinsNonzero % biasFactor != 0))));
        const int limitBase = static_cast<int>(TMath::Power(10, limitOrder));
        lowerLimitBinsCorrect = limitBase * (lowerLimitBinsNonzero / limitBase + (lowerLimitBinsNonzero % limitBase != 0));
      } else {
        lowerLimitBinsCorrect = -(alignment > 0);
      }
    } else if (lowerLimitBinsNonzero <= 0  && upperLimitBinsNonzero >= 0) {
      int upperLimitOrder = -1, lowerLimitOrder = -1;
      if (upperLimitBinsNonzero > 0) {
        upperLimitOrder = TMath::Max(0, static_cast<int>(TMath::Ceil(TMath::Log10(upperLimitBinsNonzero / biasFactor + (upperLimitBinsNonzero % biasFactor != 0)))));
      }
      if (lowerLimitBinsNonzero < 0) {
        lowerLimitOrder = TMath::Max(0, static_cast<int>(TMath::Ceil(TMath::Log10(-(lowerLimitBinsNonzero / biasFactor + (lowerLimitBinsNonzero % biasFactor != 0))))));
      }
      const int limitOrder = TMath::Max(upperLimitOrder, lowerLimitOrder);
      const int limitBase = static_cast<int>(TMath::Power(10, limitOrder));
      if (limitOrder < 0) {
        std::fprintf(stderr, "lowerLimitBinsNonzero: %d, upperLimitBinsNonzero: %d", lowerLimitBinsNonzero, upperLimitBinsNonzero);
        Fatal("adjustWithJSONTitle", "limitOrder (%d) == -1", limitOrder);
      }
      if (upperLimitBinsNonzero > 0) {
        upperLimitBinsCorrect = limitBase * (upperLimitBinsNonzero / limitBase + (upperLimitBinsNonzero % limitBase != 0));
      } else {
        upperLimitBinsCorrect = (alignment < 0);
      }
      if (lowerLimitBinsNonzero < 0) {
        lowerLimitBinsCorrect = limitBase * (lowerLimitBinsNonzero / limitBase + (lowerLimitBinsNonzero % limitBase != 0));
      } else {
        lowerLimitBinsCorrect = -(alignment > 0);
      }
    } else {
      const int nBinsNonzero = upperLimitBinsNonzero - lowerLimitBinsNonzero;
      const int limitOrder = TMath::Ceil(TMath::Log10(nBinsNonzero / biasFactor + (nBinsNonzero % biasFactor != 0)));
      const int limitBase = static_cast<int>(TMath::Power(10, limitOrder));
      if (upperLimitBinsNonzero > 0) {
        upperLimitBinsCorrect = limitBase * (upperLimitBinsNonzero / limitBase + (upperLimitBinsNonzero % limitBase != 0));
      } else {
        upperLimitBinsCorrect = (alignment < 0);
      }
      if (lowerLimitBinsNonzero < 0) {
        lowerLimitBinsCorrect = limitBase * (lowerLimitBinsNonzero / limitBase + (lowerLimitBinsNonzero % limitBase != 0));
      } else {
        lowerLimitBinsCorrect = -(alignment > 0);
      }
    }
    const double binWidth = TMath::Power(10, -binDensityOrder);
    const double lowerLimit = alignment + binWidth * lowerLimitBinsCorrect;
    const double upperLimit = alignment + binWidth * upperLimitBinsCorrect;
    const double nBinsNew = upperLimitBinsCorrect - lowerLimitBinsCorrect;
    TH1 *histNew = new TH1D(TString(hist->GetName()) + "Adjusted", hist->GetName(), nBinsNew, lowerLimit, upperLimit);
    for (int iBin = TMath::Max(lowerLimitBins, lowerLimitBinsCorrect) + 1; iBin <= TMath::Min(upperLimitBins, upperLimitBinsCorrect); ++iBin) {
      histNew->SetBinContent(iBin - lowerLimitBinsCorrect, hist->GetBinContent(iBin - lowerLimitBins));
      histNew->SetBinError(iBin - lowerLimitBinsCorrect, hist->GetBinError(iBin - lowerLimitBins));
    }
    histNew->SetBinContent(0, hist->GetBinContent(0));
    histNew->SetBinContent(nBinsNew + 1, hist->GetBinContent(upperLimitBins - lowerLimitBins + 1));
    histNew->SetEntries(hist->GetEntries());
    return histNew;
  }
  hist->SetTitle(hist->GetName());
  return hist;
}
#endif
