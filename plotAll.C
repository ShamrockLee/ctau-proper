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

#include <functional>
#include <iostream>

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
    if (obj->IsA()->InheritsFrom("TH1")) {
      TH1 *hist = (TH1 *) key->ReadObj();
      if (normalize) {
        hist->Scale(1.0 / hist->Integral());
      }
      if (funAdjustHist != nullptr) hist = funAdjustHist(hist);
      hist->Draw(optionDraw);
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
