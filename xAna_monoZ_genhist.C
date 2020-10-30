#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector.h>
#include <TVectorFfwd.h>

#include <iostream>

void xAna_monoZ_genhist(TString outputFileHead, TString outputFileVar,
                     TString outputFileTail, Bool_t toRecreateOutFile = true,
                     Bool_t debug = false) {
  const TString namesLepton[] = {"Electron", "Muon",
                                 "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon",
                                      "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};

  /// Tree for GEN-correct  electron/muon/tau
  TTree *arrTTGen[3];

  /// Tree for events with correct numbers of electron/muon
  TTree *arrTTNumCorrect[2];

  /// Tree for events whose electron/muon pairs pass Z mass cut
  TTree *arrTTZMassCutted[2];

  Long64_t nEntryOriginal;

  auto getTrees =
      [&](TFile *tf) {
        for (Byte_t i = 0; i < 3; i++) {
          arrTTGen[i] = (TTree *)tf->Get((TString) "Gen" + namesLepton[i]);
        }
        for (Byte_t i = 0; i < 2; i++) {
          arrTTNumCorrect[i] =
              (TTree *)tf->Get((TString) "NumCorrect" + namesLepton[i]);
        }
        for (Byte_t i = 0; i < 2; i++) {
          arrTTZMassCutted[i] =
              (TTree *)tf->Get((TString) "ZMassCutted" + namesLepton[i]);
        }
        nEntryOriginal = (*((TVectorD *) tf->Get("tvdNEntryOriginal")))[0];
      };
  // TTree* ttGenElectron = new TTree("GenElectron", "GEN-level electron
  // TString outputFileHist = (TString) "output_" + outputFileHead + "_" +
  //                          outputFileVar + "_" + outputFileTail + "_hist" +
  //                          ".root";

  // TFile *outFileHist = new TFile(outputFileHist.Data(), toRecreateOutFile ?
  // "recreate" : "update");

  TString outImageDir =
      (TString) "../out_images/output_" + outputFileHead + "_" + outputFileTail;

  // TCanvas* c1 = new TCanvas;
  gStyle->SetOptStat(111111);
  auto lambdaPrintHistograms = [/*&c1,*/ /*&outFileHist,*/ toRecreateOutFile,
                                outImageDir, outputFileHead, outputFileVar,
                                outputFileTail](TTree *tt) -> void {
    std::cout << "Printing " << tt->GetName() << std::endl;
    TString nameTree = tt->GetName();
    TString outImageNameHead =
        (TString)outputFileHead + "_" + outputFileVar + "_" + nameTree;
    TString outImageCommonPath = outImageDir + "/" + outImageNameHead + "_";
    TString outputFileHist = (TString) "output_" + outputFileHead + "_" +
                             outputFileVar + "_" + outputFileTail + "_hist_" +
                             nameTree + ".root";
    TFile *outFileHist = new TFile(outputFileHist.Data(),
                                   toRecreateOutFile ? "recreate" : "update");
    for (TObject *leafObject : *(tt->GetListOfLeaves())) {
      TLeaf *leaf = (TLeaf *)leafObject;
      TString nameLeaf = leaf->GetName();
      if (nameLeaf.First("[") >= 0) {
        nameLeaf.Resize(nameLeaf.First("["));
      }
      TString nameHist = "h" + nameLeaf;
      tt->Draw(nameLeaf + ">>" + nameHist);
      TH1 *hist = (TH1 *)gDirectory->Get(nameHist);
      hist->SetTitle((TString)leaf->GetBranch()->GetTitle() + " (" + nameTree +
                     ")");
      TString outImagePath = outImageCommonPath + nameLeaf + ".svg";
      // c1->Clear();
      // hist->Draw();
      // c1->Print(outImagePath);
      outFileHist->WriteObject(hist, nameHist);
    }
    outFileHist->Close();
  };
  if (true) {
    for (Byte_t i = 0; i < 3; i++) {
      lambdaPrintHistograms(arrTTGen[i]);
    }
  }
  // for (Byte_t i=0; i<2; i++) {
  //   lambdaPrintHistograms(arrTTNumCorrect[i]);
  // }
  for (Byte_t i = 0; i < 2; i++) {
    lambdaPrintHistograms(arrTTZMassCutted[i]);
  }
  // c1->Close();
}