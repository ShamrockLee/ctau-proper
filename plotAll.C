#include<TROOT.h>
#include<TSystem.h>
#include<TDirectory.h>
#include<TFile.h>
#include<TLeaf.h>
#include<TPad.h>
#include<TCanvas.h>
#include<TString.h>
#include<TH1.h>
#include<TLegend.h>

void plotAll(TString pathFileIn, TString dirOut, Bool_t normalize, Bool_t logy, Option_t *optionDraw);

void plotAll(const char* pathFileIn, const char* dirout, Bool_t normalize = false, Bool_t logy = false, Option_t *optionDraw = "") {
  plotAll((TString) pathFileIn, (TString) dirout, normalize, logy, optionDraw);
}

void plotAll(TString pathFileIn, TString dirOut, Bool_t normalize = false, Bool_t logy = false, Option_t *optionDraw = "") {
  TString seperatorPath = "/";
  TFile *fileIn = TFile::Open(pathFileIn);
  TString basenameFileIn = gSystem->GetFromPipe(
      Form("file=%s; test=${file##*/}; echo \"${test%%.root}\"",
      pathFileIn.Data()));
  gSystem->mkdir(dirOut);
  if (dirOut.Length() > 1 && dirOut.EndsWith(seperatorPath)) {
    dirOut.Resize(dirOut.Length()-1);
  }
  if (!gPad) {
    gROOT->MakeDefCanvas();
  }
  gPad->SetLogy(logy);
  for (auto&& keyRaw: *fileIn->GetListOfKeys()) {
    TKey *key = (TKey *) keyRaw;
    if (key == nullptr) continue;
    TString name = key->GetName();
    TObject *obj = key->ReadObj();
    if (obj == nullptr) continue;
    if (!obj->IsA()->InheritsFrom("TH1")) continue;
    TH1 *hist = (TH1 *) key->ReadObj();
    if (normalize) {
      hist->Scale(1.0 / hist->Integral());
    }
    hist->Draw(optionDraw);
    gPad->Print(dirOut + seperatorPath + basenameFileIn + "_" + name + ".svg");
  }
  fileIn->Close();
}

