#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TVector.h>
#include <TVectorFfwd.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TObjArray.h>

#include <iostream>
#include <fstream>
#include <vector>

TString modifyNameLeaf(TString nameLeaf) {
  TString nameLeafNew = TString(nameLeaf);
  Ssiz_t indexLeftSquareBracket = nameLeaf.First('[');
  if (indexLeftSquareBracket != -1) {
    nameLeafNew.Resize(indexLeftSquareBracket);
  }
  return nameLeafNew;
}

void xAna_monoZ_genhist(TString nameCondorPack, TString nameDatagroup, TString nameClusterID,
                     TString dirCondorPacks = ".", Bool_t toRecreateOutFile = true,
                     Bool_t debug = false) {
  const TString namesLepton[] = {"Electron", "Muon",
                                 "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon",
                                      "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};

  // /// Tree for GEN-correct  electron/muon/tau
  // TTree *arrTTGen[3];

  // /// Tree for events with correct numbers of electron/muon
  // TTree *arrTTNumCorrect[2];

  // /// Tree for events whose electron/muon pairs pass Z mass cut
  // TTree *arrTTZMassCutted[2];

  Long64_t nEntryOriginal;

  // auto mountTrees =
  //     [&arrTTGen, &arrTTNumCorrect, &arrTTZMassCutted, &nEntryOriginal, namesLepton](TFile *tf) {
  //       for (Byte_t i = 0; i < 3; i++) {
  //         arrTTGen[i] = (TTree *)tf->Get((TString) "Gen" + namesLepton[i]);
  //       }
  //       for (Byte_t i = 0; i < 2; i++) {
  //         arrTTNumCorrect[i] =
  //             (TTree *)tf->Get((TString) "NumCorrect" + namesLepton[i]);
  //       }
  //       for (Byte_t i = 0; i < 2; i++) {
  //         arrTTZMassCutted[i] =
  //             (TTree *)tf->Get((TString) "ZMassCutted" + namesLepton[i]);
  //       }
  //       nEntryOriginal = (*((TVectorD *) tf->Get("tvdNEntryOriginal")))[0];
  //     };

  // // TTree* ttGenElectron = new TTree("GenElectron", "GEN-level electron
  // // TString outputFileHist = (TString) "output_" + outputFileHead + "_" +
  // //                          outputFileVar + "_" + outputFileTail + "_hist" +
  // //                          ".root";

  // // TFile *outFileHist = new TFile(outputFileHist.Data(), toRecreateOutFile ?
  // // "recreate" : "update");

  // TString outImageDir =
  //     (TString) "../out_images/output_" + outputFileHead + "_" + outputFileTail;

  std::vector<TString> vNameTT;
  vNameTT.clear();
  for (TString nameLeptonCurrent: namesLepton) {
    vNameTT.push_back((TString) "ZMassCutted" + nameLeptonCurrent);
  }

  const UInt_t nTT = vNameTT.size();

  TTree *arrTT[nTT];

  auto mountVariablesFromFile = [&nEntryOriginal, &arrTT, nTT, vNameTT] (TFile *tf) {
    nEntryOriginal = (*((TVectorD *) tf->Get("tvdNEntryOriginal")))[0];
    for (UInt_t i=0; i<nTT; i++) {
      arrTT[i] = (TTree *) tf->Get(vNameTT[i]);
    }
  };

  // // TCanvas* c1 = new TCanvas;
  // gStyle->SetOptStat(111111);
  // auto lambdaPrintHistograms = [/*&c1,*/ /*&outFileHist,*/ toRecreateOutFile,
  //                               outImageDir, outputFileHead, outputFileVar,
  //                               outputFileTail](TTree *tt) -> void {
  //   std::cout << "Printing " << tt->GetName() << std::endl;
  //   TString nameTree = tt->GetName();
  //   TString outImageNameHead =
  //       (TString)outputFileHead + "_" + outputFileVar + "_" + nameTree;
  //   TString outImageCommonPath = outImageDir + "/" + outImageNameHead + "_";
  //   TString outputFileHist = (TString) "output_" + outputFileHead + "_" +
  //                            outputFileVar + "_" + outputFileTail + "_hist_" +
  //                            nameTree + ".root";
  //   TFile *outFileHist = new TFile(outputFileHist.Data(),
  //                                  toRecreateOutFile ? "recreate" : "update");
  //   for (TObject *leafObject : *(tt->GetListOfLeaves())) {
  //     TLeaf *leaf = (TLeaf *)leafObject;
  //     TString nameLeaf = leaf->GetName();
  //     // if (nameLeaf.First('[') >= 0) {
  //     //   nameLeaf.Resize(nameLeaf.First("["));
  //     // }
  //     TString nameHist = "h" + nameLeaf;
  //     tt->Draw(nameLeaf + ">>" + nameHist);
  //     TH1 *hist = (TH1 *)gDirectory->Get(nameHist);
  //     hist->SetTitle((TString)leaf->GetBranch()->GetTitle() + " (" + nameTree +
  //                    ")");
  //     TString outImagePath = outImageCommonPath + nameLeaf + ".svg";
  //     // c1->Clear();
  //     // hist->Draw();
  //     // c1->Print(outImagePath);
  //     outFileHist->WriteObject(hist, nameHist);
  //   }
  //   outFileHist->Close();
  // };
  // if (true) {
  //   for (Byte_t i = 0; i < 3; i++) {
  //     lambdaPrintHistograms(arrTTGen[i]);
  //   }
  // }
  // // for (Byte_t i=0; i<2; i++) {
  // //   lambdaPrintHistograms(arrTTNumCorrect[i]);
  // // }
  // for (Byte_t i = 0; i < 2; i++) {
  //   lambdaPrintHistograms(arrTTZMassCutted[i]);
  // }
  // // c1->Close();

  std::vector<UInt_t> vecNFile;
  TString dirCondorPackCurrent = dirCondorPacks + "/" + nameCondorPack;
  TString pathTextFileNumbers = dirCondorPackCurrent + "/" + "fileNumbers_" + nameDatagroup + ".txt";
  TString pathTextCrossSections = dirCondorPackCurrent + "/" + "crossSections_" + nameDatagroup + ".txt";
  std::vector<UInt_t> vNumberFile;
  std::vector<Double_t> vCrossSection;
  UInt_t nDataset;
  UInt_t nFileTot;
  {
    std::ifstream infileText;
    infileText.open(pathTextFileNumbers);
    UInt_t valNumberFileCurrent;
    if (!infileText) {
      std::cerr << "Unable to open file " << pathTextFileNumbers << std::endl;
    }
    while (!infileText) {
      if (infileText >> valNumberFileCurrent) {
        vNumberFile.push_back(valNumberFileCurrent);
      } else {
        std::cerr << "One line failed to parse." << std::endl;
      }
    }
    infileText.close();
    infileText.open(pathTextCrossSections);
    Double_t valCrossSectionCurrentt;
    if (!infileText) {
      std::cerr << "Unable to open file " << pathTextCrossSections << std::endl;
    }
    while (!infileText) {
      if (infileText >> valCrossSectionCurrentt) {
        vCrossSection.push_back(valNumberFileCurrent);
      } else {
        std::cerr << "One line failed to parse." << std::endl;
      }
    }
    infileText.close();
    nDataset = vNumberFile.size();
    if (nDataset != vCrossSection.size()) {
      std::cerr << "The length of fileNumbers (" << nDataset << ") does not match that of the crossSection (" << vCrossSection.size() << ")" << std::endl;
    }
    nDataset = TMath::Min(nDataset, (UInt_t) vCrossSection.size());
    nFileTot = 0;
    for (UInt_t iDataset=0; iDataset<nDataset; iDataset++) {
      nFileTot += vNumberFile[iDataset];
    }
  }

  TCanvas *c1;
  TString filenameGeneralCurrentCluster = "output_" + nameDatagroup + "_" + nameClusterID;
  TObjArray* toatltlHistFileLeafTree = new TObjArray(nTT);
  TVectorD tvdWeightDataset(nDataset);
  TVectorD tvdWeightFile(nFileTot);
  std::vector<std::vector<std::vector<Bool_t>>> vvvIsHistFileLeafTree(nTT);
  std::vector<std::vector<TString>> vvNameLeafTree(nTT);
  for (UInt_t iTree=0; iTree<nTT; iTree++) {
    std::vector<std::vector<Bool_t>> vvIsHistFileLeaf;
    vvvIsHistFileLeafTree.push_back(vvIsHistFileLeaf);
    std::vector<TString> vNameLeaf;
    vvNameLeafTree.push_back(vNameLeaf);
    TList *tltlHistFileLeaf = new TList;
    (*toatltlHistFileLeafTree)[iTree] = tltlHistFileLeaf;
  }
  UInt_t iFile = 0;
  for (UInt_t iDataset=0; iDataset<nDataset; iDataset++) {
    Long64_t nEntryOriginalDatasetCurrent = 0;
    for (UInt_t jFile=0; jFile<vNumberFile[iDataset]; jFile++, iFile++) {
      TString filenameCurrent = filenameGeneralCurrentCluster + "_" + iFile + ".root";
      TFile *tfCurrent = TFile::Open(dirCondorPackCurrent + "/" + filenameCurrent, "READ");
      mountVariablesFromFile(tfCurrent);
      nEntryOriginalDatasetCurrent += nEntryOriginal;
      for (UInt_t iTree=0; iTree=nTT; iTree++) {
        TList *tltlHistFileLeaf = (TList *) (*toatltlHistFileLeafTree)[iTree];
        for (TObject *leafRaw: *(arrTT[iTree]->GetListOfLeaves())) {
          TLeaf *leaf = (TLeaf *) leafRaw;
          TString nameLeaf = leaf->GetName();
          TString nameLeafModified = modifyNameLeaf(nameLeaf);
          Bool_t isIncluded = false;
          UInt_t indexName = 0;
          // TList *tlHistFileCurrent;
          for (TString nameIncluded: vvNameLeafTree[iTree]) {
            if (nameLeafModified == nameIncluded) {
              isIncluded = true;
              break;
            }
            indexName++;
          }
          if (!isIncluded) {
            vvNameLeafTree[iTree].push_back(nameLeafModified);
            std::vector<Bool_t> vIsHistFile(nFileTot);
            for (UInt_t i=0; i<nFileTot; i++) {
              vIsHistFile.push_back(false);
            }
            vvvIsHistFileLeafTree[iTree].push_back(vIsHistFile);
            TList *tlHistFileCurrent = new TList;
            tltlHistFileLeaf->AddLast(tlHistFileCurrent);
          }
          vvvIsHistFileLeafTree[iTree][indexName][iFile] = true;
          TString nameHist = (TString)"h" + vNameTT[iTree] + nameLeafModified + iFile;
          arrTT[iTree]->Draw(nameLeaf + ">>" + nameHist);
        }
        Int_t indexName = 0;
        for (auto tlHistFileRaw: *tltlHistFileLeaf) {
          TList *tlHistFile = (TList *) tlHistFileRaw;
          TString nameLeafModified = vvNameLeafTree[iTree][indexName];
          TString nameHist = (TString)"h" + vNameTT[iTree] + nameLeafModified + iFile;
          tlHistFile->AddLast(gDirectory->Get(nameHist));
        }
      }
    }
    tvdWeightDataset[iDataset] = vCrossSection[iDataset]/nEntryOriginalDatasetCurrent;
    for (UInt_t jFile=0; jFile<vNumberFile[iDataset]; jFile++) {
      tvdWeightFile[iFile - vNumberFile[iDataset] + jFile] = tvdWeightDataset[iDataset];
    }
  }

  TObjArray *toatlResult = new TObjArray(nTT);
  for (UInt_t iTree=0; iTree<nTT; iTree++) {
    (*toatlResult)[iTree] = new TList;
    UInt_t indexName = 0;
    for (auto tlHistFileLeafRaw: toatltlHistFileLeafTree[iTree]) {
      TList *tlHistFileLeaf = (TList *) tlHistFileLeafRaw;
      TH1 *histFirst = (TH1 *) tlHistFileLeaf->First();
      TH1 *histResult = (TH1 *) histFirst->Clone("h" + vvNameLeafTree[iTree][indexName]);
      UInt_t iDataset = 0;
      UInt_t iFile = 0;
      for (auto histFileRaw: *tlHistFileLeaf) {
        TH1 *histFile = (TH1 *)histFileRaw;
        histFile->Scale(tvdWeightFile[iFile]);
      }
      histResult->Merge(tlHistFileLeaf);
      ((TList *)(*toatlResult)[iTree])->AddLast(histResult);
    }
  }
  TFile *tfOut = TFile::Open(dirCondorPackCurrent + "/" + "output_" + nameDatagroup + "_" + nameClusterID + "_hist.root", toRecreateOutFile ? "recreate" : "update");
  toatlResult->Write("toatlResult");
  tfOut->Close();

}