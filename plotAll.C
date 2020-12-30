#include<TSystem.h>
#include<TDirectory.h>
#include<TFile.h>
#include<TLeaf.h>
#include<TPad.h>
#include<TCanvas.h>
#include<TString.h>
#include<TH1.h>
#include<TLegend.h>

void plotAll(TString pathFileIn, TString dirOut);

void plotAll(const char* pathFileIn, const char* dirout) {
  plotAll((TString) pathFileIn, (TString) dirout);
}

void plotAll(TString pathFileIn, TString dirOut) {
  TString seperatorPath = "/";
  TFile *fileIn = TFile::Open(pathFileIn);
  TString basenameFileIn = gSystem->GetFromPipe(
      Form("file=%s; test=${file##*/}; echo \"${test%%.root}\"",
      pathFileIn.Data()));
  if (dirOut.Length() > 1 && dirOut.EndsWith(seperatorPath)) {
    dirOut.Resize(dirOut.Length()-1);
  }
  for (auto&& keyRaw: *fileIn->GetListOfKeys()) {
    TKey *key = (TKey *) keyRaw;
    if (key == nullptr) continue;
    TString name = key->GetName();
    TObject *obj = key->ReadObj();
    if (obj == nullptr) continue;
    if (!obj->IsA()->InheritsFrom("TH1")) continue;
    TH1 *hist = (TH1 *) key->ReadObj();
    hist->Draw();
    gPad->Print(dirOut + seperatorPath + basenameFileIn + "_" + name + ".svg");
  }
  fileIn->Close();
}

