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
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

template<typename TBase, typename TIndex, typename TRatio>
TRatio intpow(const TBase base, const TIndex index, const TRatio ratio) {
  if (index == 0) {
    return ratio;
  }
  TRatio result = ratio;
  if (index > 0) {
    for (TIndex i=0; i<index; i++) {
      result *= base;
    }
  } else {
    for (TIndex i=0; i<-index; i++) {
      result /= base;
    }
  }
  return result;
}

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
  for (UShort_t i=0; i<2; i++) {
    TString nameLeptonCurrent = namesLepton[i];
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
      infileText >> valNumberFileCurrent;
      // if (infileText >> valNumberFileCurrent) {
      vNumberFile.push_back(valNumberFileCurrent);
      // } else {
      //   std::cerr << "One line failed to parse." << std::endl;
      // }
    }
    infileText.close();
    vNumberFile = {3, 3}; //TODO
    infileText.open(pathTextCrossSections);
    Double_t valCrossSectionCurrent;
    if (!infileText) {
      std::cerr << "Unable to open file " << pathTextCrossSections << std::endl;
    }
    while (!infileText) {
      infileText >> valCrossSectionCurrent;
      // if (infileText >> valCrossSectionCurrentt) {
      vCrossSection.push_back(valNumberFileCurrent);
      // } else {
      //   std::cerr << "One line failed to parse." << std::endl;
      // }

    }
    infileText.close();
    vCrossSection = {730.6, 730.6}; //TODO
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

  if (debug) std::cout << "nFileTot: " << nFileTot << std::endl;
  TCanvas *c1;
  // TString filenameGeneralCurrentCluster = "output_" + nameDatagroup + "_" + nameClusterID;
  auto getFilenameCurrent = [nameDatagroup, nameClusterID](UInt_t iFile) -> TString {
    return "output_" + nameDatagroup + "_" + nameClusterID + "_" + iFile + ".root";
  };
  TFile *arrFile[nFileTot];
  // TObjArray* toatltlHistFileLeafTree = new TObjArray(nTT); //< The histogram in each file of each leaf in each tree
  // std::vector<std::vector<std::vector<TH1 *>>> vvvHistFileLeafTree(nTT);
  TVectorD tvdWeightDataset(nDataset); //< The weight of each dataset
  TVectorD tvdWeightFile(nFileTot); //< The weight of each file
  std::vector<std::vector<std::vector<Bool_t>>> vvvIsHistFileLeafTree(nTT); //< The existence in each file of each leaf in each tree
  // std::vector<std::vector<std::vector<TString>>> vvvNameEffectiveFileLeafTree(nTT); //< The leafname in each file of each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Int_t>>> vvvNEntryFileLeafTree(nTT); //< The number of entries of the histogram in each file of each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Double_t>>> vvvHistMinFileLeafTree(nTT); //< The minimum value of each histogram in each file of each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Double_t>>> vvvHistMaxFileLeafTree(nTT); //< The maximum value of each histogram in each file of each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Int_t>>> vvvNBinFileLeafTree(nTT); //< The bin-number of the histogram in each file of each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Double_t>>> vvvLowerFileLeafTree(nTT); //< The lower limit of the histogram in each file in each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Double_t>>> vvvUpperFileLeafTree(nTT); //< The upper limit of the histogram in each file in each (existing) leaf in each tree
  std::vector<std::vector<TString>> vvNameModifiedLeafTree(nTT); //< The modified leafname of each leaf in each tree
  std::vector<std::vector<TString>> vvTitleLeafTree(nTT); // The title of the branch of each leaf in each tree
  std::vector<std::vector<TString>> vvTypeNameLeafTree(nTT); // The type name of each each leaf in each tree
  auto getNameHistOfFile = [&vNameTT, &vvNameModifiedLeafTree](UInt_t iTree, UInt_t indexName, UInt_t iFile, TString suffixStage) -> TString {
    return (TString)"h" + vNameTT[iTree] + vvNameModifiedLeafTree[iTree][indexName] + suffixStage + iFile;
  };
  for (UInt_t iTree=0; iTree<nTT; iTree++) {
    std::vector<std::vector<Bool_t>> vvIsHistFileLeaf;
    vvvIsHistFileLeafTree.push_back(vvIsHistFileLeaf);
    std::vector<std::vector<TH1 *>> vvHistFileLeaf;
    // vvvHistFileLeafTree.push_back(vvHistFileLeaf);
    // std::vector<std::vector<TString>> vvNameEffectiveFileLeaf;
    // vvvNameEffectiveFileLeafTree.push_back(vvNameEffectiveFileLeaf);
    // std::vector<std::vector<Int_t>> vvNEntryFileLeaf;
    // vvvNBinFileLeafTree.push_back(vvNEntryFileLeaf);
    // std::vector<std::vector<Double_t>> vvHistMinFileLeaf;
    // vvvHistMinFileLeafTree.push_back(vvHistMinFileLeaf);
    // std::vector<std::vector<Double_t>> vvHistMaxFileLeaf;
    // vvvHistMaxFileLeafTree.push_back(vvHistMaxFileLeaf);
    // std::vector<std::vector<Int_t>> vvNBinFileLeaf;
    // vvvNBinFileLeafTree.push_back(vvNBinFileLeaf);
    // std::vector<std::vector<Double_t>> vvLowerFileLeaf;
    // vvvLowerFileLeafTree.push_back(vvLowerFileLeaf);
    // std::vector<std::vector<Double_t>> vvUpperFileLeaf;
    // vvvUpperFileLeafTree.push_back(vvUpperFileLeaf);
    std::vector<TString> vNameModifiedLeaf;
    vvNameModifiedLeafTree.push_back(vNameModifiedLeaf);
    std::vector<TString> vTitleLeaf;
    vvTitleLeafTree.push_back(vTitleLeaf);
    std::vector<TString> vTypeNameLeaf;
    vvTypeNameLeafTree.push_back(vTypeNameLeaf);
    // TList *tltlHistFileLeaf = new TList;
    // (*toatltlHistFileLeafTree)[iTree] = tltlHistFileLeaf;
  }
  // UInt_t iFile = 0;
  TFile *tfAutogenHist = TFile::Open(dirCondorPackCurrent + "/" + "output_" + nameDatagroup + "_" + "Autogen" + "_" + nameClusterID + "_hist.root", toRecreateOutFile ? "recreate" : "update");
  for (UInt_t iDataset=0, iFile=0; iDataset<nDataset; iDataset++) {
    Long64_t nEntryOriginalDatasetCurrent = 0;
    if (debug) std::cout << "iDataset: " << iDataset << std::endl;
    for (UInt_t jFile=0; jFile<vNumberFile[iDataset]; jFile++, iFile++) {
      // TString filenameCurrent = filenameGeneralCurrentCluster + "_" + iFile + ".root";
      TString filenameCurrent = getFilenameCurrent(iFile);
      if (debug) std::cout << "Opening TFile " << filenameCurrent << " ..." << std::endl;
      TFile *tfCurrent = TFile::Open(dirCondorPackCurrent + "/" + "dataTest_" + nameDatagroup + "/" + filenameCurrent, "READ");
      arrFile[iFile] = tfCurrent;
      if (debug) std::cout << "Done." << std::endl;
      if (debug) std::cout << "Mounting variables from the file ..." << std::endl;
      mountVariablesFromFile(tfCurrent);
      if (debug) std::cout << "Done." << std::endl;
      if (debug) std::cout << "nEntryOriginal: " << nEntryOriginal << std::endl;
      nEntryOriginalDatasetCurrent += nEntryOriginal;
      for (UInt_t iTree=0; iTree<nTT; iTree++) {
        if (debug) std::cout << "iTree: " << iTree << std::endl;
        if (debug) std::cout << "Getting leafnames from the tree ..." << std::endl;
        // TList *tltlHistFileLeaf = (TList *) (*toatltlHistFileLeafTree)[iTree];
        for (TObject *leafRaw: *(arrTT[iTree]->GetListOfLeaves())) {
          TLeaf *leaf = (TLeaf *) leafRaw;
          TString nameLeaf = leaf->GetName();
          TString nameLeafModified = modifyNameLeaf(nameLeaf);
          Bool_t isIncluded = false;
          UInt_t indexName = 0;
          // TList *tlHistFileCurrent;
          for (TString nameIncluded: vvNameModifiedLeafTree[iTree]) {
            // if (debug) std::cout << "indexName : " << indexName;
            if (nameLeafModified == nameIncluded) {
              isIncluded = true;
              // if (debug) std::cout << " is included.";
              break;
            }
            indexName++;
          }
          if (debug) std::cout << std::endl;
          if (!isIncluded) {
            if (debug) std::cout << "Found new leafname." << std::endl;
            if (debug) std::cout << "nameLeaf: " << nameLeaf << ", nameLeafModified:" << nameLeafModified << std::endl;
            vvNameModifiedLeafTree[iTree].push_back(nameLeafModified);
            vvTitleLeafTree[iTree].push_back(leaf->GetBranch()->GetTitle());
            vvTypeNameLeafTree[iTree].push_back(leaf->GetTypeName());
            std::vector<Bool_t> vIsHistFile(nFileTot);
            std::vector<TH1 *> vHistFile(nFileTot);
            for (UInt_t i=0; i<nFileTot; i++) {
              vIsHistFile.push_back(false);
            }
            vvvIsHistFileLeafTree[iTree].push_back(vIsHistFile);
            // vvvHistFileLeafTree[iTree].push_back(vHistFile);
            // TList *tlHistFileCurrent = new TList;
            // tltlHistFileLeaf->AddLast(tlHistFileCurrent);
            
          }
          vvvIsHistFileLeafTree[iTree][indexName][iFile] = true;
          TString nameHistAutogen = getNameHistOfFile(iTree, indexName, iFile, "Autogen");
          if (debug) std::cout << "Drawing " << nameLeaf << " as " << nameHistAutogen << " ...";
          arrTT[iTree]->Draw(nameLeaf + ">>" + nameHistAutogen);
          TH1 *histAutogen = (TH1 *)gDirectory->Get(nameHistAutogen);
          histAutogen->SetTitle(nameLeaf);
          histAutogen->SetName(nameHistAutogen);
          histAutogen->SetDirectory(tfAutogenHist);
          // vvvHistFileLeafTree[iTree][indexName].push_back(histAutogen);
          if (debug && histAutogen == nullptr) std::cerr << "Fatal: " << nameHistAutogen << " is nullptr!" << std::endl;
          if (debug) std::cout << " " << histAutogen;
          if (debug) std::cout << " Done." << std::endl;
        }
        // Int_t indexName = 0;
        // for (const auto&& tlHistFileRaw: *tltlHistFileLeaf) {
        //   // if (debug) std::cout << "Looping." << std::endl;
        //   if (vvvIsHistFileLeafTree[iTree][indexName][iFile]) {
        //     TList *tlHistFile = (TList *) tlHistFileRaw;
        //     TString nameLeafModified = vvNameModifiedLeafTree[iTree][indexName];
        //     TString nameHist = getNameHistOfFile(iTree, indexName, iFile, "Autogen");
        //     if (debug) std::cout << "Picking up " << nameHist << " ...";
        //     TH1* hist = (TH1 *) gDirectory->Get(nameHist);
        //     if (debug) std::cout << " (" << hist->GetNbinsX()
        //     << "," << hist->GetBinContent(1) 
        //     << "," << hist->GetBinContent(hist->GetNbinsX()) << ")";
        //     tlHistFile->AddLast(hist);
        //     if (debug) std::cout << " Done." << std::endl;
        //   }
        //   indexName++;
        // }
      }
      if (debug) std::cout << "Done getting names." << std::endl;
    }
    if (debug) std::cout << "nEntryOriginalDatasetCurrent: " << nEntryOriginalDatasetCurrent << std::endl;
    tvdWeightDataset[iDataset] = vCrossSection[iDataset]/nEntryOriginalDatasetCurrent;
    if (debug) std::cout << "tvdWeightDataset[" << iDataset << "]: " << tvdWeightDataset[iDataset] << "," << std::endl;
    if (debug) std::cout << "composed of files index";
    for (UInt_t jFile=0; jFile<vNumberFile[iDataset]; jFile++) {
      UInt_t indexFile = iFile - vNumberFile[iDataset] + jFile;
      if (debug) std::cout << " " << indexFile;
      tvdWeightFile[indexFile] = tvdWeightDataset[iDataset];
    }
    if (debug) std::cout << "." << std::endl;
  }
  tfAutogenHist->Write();
  tfAutogenHist->Close();
  delete tfAutogenHist;

  tfAutogenHist = TFile::Open(dirCondorPackCurrent + "/" + "output_" + nameDatagroup + "_" + "Autogen" + "_" + nameClusterID + "_hist.root", "read");
  gStyle->SetOptStat(111111);
  // TObjArray *toatlResult = new TObjArray(nTT);
  std::vector<std::vector<Bool_t>> vvIsAllEmptyLeafTree(nTT); //< Whether all the histograms have zero entry of each leaf in each tree
  std::vector<std::vector<TString>> vvHistsettingLeafTree(nTT); //< The histogram settings "(binnumber, lower, upper)" of each leaf in each tree;
  for (UInt_t iTree=0; iTree<nTT; iTree++) {
    TFile *tfOutHist = TFile::Open(dirCondorPackCurrent + "/" + "output_" + nameDatagroup + "_" + vNameTT[iTree] + "_" + nameClusterID + "_hist.root", toRecreateOutFile ? "recreate" : "update");
    TFile *tfCorrectedHist = TFile::Open(dirCondorPackCurrent + "/" + "output_" + nameDatagroup + "_" + "Corrected" + "_" + nameClusterID + "_hist.root", toRecreateOutFile ? "recreate" : "update");
    {
      std::vector<Bool_t> vIsAllEmptyLeaf;
      vvIsAllEmptyLeafTree.push_back(vIsAllEmptyLeaf);
      std::vector<TString> vHistsettingLeaf;
      vvHistsettingLeafTree.push_back(vHistsettingLeaf);
    }
    // (*toatlResult)[iTree] = new TList;
    // UInt_t indexName = 0;
    // for (const auto&& tlHistFileRaw: *((TList *) (*toatltlHistFileLeafTree)[iTree])) {
    UInt_t nLeaf = vvNameModifiedLeafTree[iTree].size();
    for(UInt_t indexName=0; indexName < nLeaf; indexName++) {
      if (debug) std::cout << "indexName: " << indexName << std::endl;
      // TList *tlHistFile = (TList *) tlHistFileRaw;
      // TH1 *histFirst = (TH1 *) tlHistFile->First();

      // Index of the first true element of vvvIsHistFileLeafTree[iTree][indexName]
      UInt_t iFileFirst = std::distance(vvvIsHistFileLeafTree[iTree][indexName].begin(), std::find(vvvIsHistFileLeafTree[iTree][indexName].begin(), vvvIsHistFileLeafTree[iTree][indexName].end(), true));
      //TH1 *histFirst = (TH1 *) gDirectory->Get(getNameHistOfFile(iTree, indexName, iFileFirst, "Autogen"));
      TH1 *histFirst = (TH1 *) tfAutogenHist->Get(getNameHistOfFile(iTree, indexName, iFileFirst, "Autogen"));
      // TH1 *histFirst = (TH1 *) vvvHistFileLeafTree[iTree][indexName][0];
      // TH1 *histResult = (TH1 *) histFirst->Clone("h" + vvNameModifiedLeafTree[iTree][indexName]);
      // histResult->Clear();
      TH1 *histResult;
      TString nameHistResult = "h" + vvNameModifiedLeafTree[iTree][indexName];
      // UInt_t iDataset = 0;
      // UInt_t iFile = 0;
      // UInt_t nHistAllFile = tlHistFile->GetEntries();
      UInt_t nHistAllFile = std::count(vvvIsHistFileLeafTree[iTree][indexName].begin(), vvvIsHistFileLeafTree[iTree][indexName].end(), true);
      if (debug) std::cout << "nHistAllFile: " << nHistAllFile << std::endl;
      Int_t arrNEntryFile[nHistAllFile];
      Double_t arrMinimumFile[nHistAllFile];
      Double_t arrMaximumFile[nHistAllFile];
      Int_t arrNBinFile[nHistAllFile];
      Double_t arrLowerFile[nHistAllFile];
      Double_t arrUpperFile[nFileTot];
      TString arrLeafnameFile[nFileTot];
      Bool_t isAllEmpty = true;
      // UInt_t iHist = 0;
      if (debug) std::cout << "Analyzing autogen histograms for " << vvNameModifiedLeafTree[iTree][indexName] << "...";
      // for (auto histFileRaw: *tlHistFile) {
      for (UInt_t iFile=0, iHist=0; iFile<nFileTot; iFile++) {
        if (!vvvIsHistFileLeafTree[iTree][indexName][iFile]) {
          if (debug) std::cout << " Skipping iFile: " << iFile;
          continue;
        }
        if (debug) std::cout << " iHist: " << iHist;
        // TH1 *histFile = (TH1 *)histFileRaw;
        TString nameHistAutogen = getNameHistOfFile(iTree, indexName, iFile, "Autogen");
        // TH1 *histAutogen = (TH1 *) gDirectory->Get(nameHistAutogen);
        TH1 *histAutogen = (TH1 *) tfAutogenHist->Get(nameHistAutogen);
        // TH1 *histAutogen = vvvHistFileLeafTree[iTree][indexName][iHist];
        if (debug && histAutogen == nullptr) std::cerr << "Fatal: histAutogen from " << nameHistAutogen << " is nullptr!" << std::endl;
        // histFile->Scale(tvdWeightFile[iFile]);-
        arrNEntryFile[iHist] = histAutogen->GetEntries();
        arrMinimumFile[iHist] = histAutogen->GetBinCenter(histAutogen->FindFirstBinAbove()) - histAutogen->GetBinWidth(1) / 2;
        arrMaximumFile[iHist] = histAutogen->GetBinCenter(histAutogen->FindLastBinAbove()) + histAutogen->GetBinWidth(1) / 2;
        Int_t nBin = histAutogen->GetNbinsX();
        arrNBinFile[iHist] = nBin;
        if (nBin > 0) {
          isAllEmpty = false;
        }
        arrLowerFile[iHist] = histAutogen->GetBinCenter(1) - histAutogen->GetBinWidth(1) / 2;
        arrUpperFile[iHist] = histAutogen->GetBinCenter(nBin) + histAutogen->GetBinWidth(1) / 2;
        arrLeafnameFile[iHist] = histAutogen->GetTitle();
        iHist++;
      }
      if (debug) std::cout << " Done." << std::endl;
      if (isAllEmpty) {
        if (debug) std::cout << "Empty leaf encountered: " << vvNameModifiedLeafTree[iTree][indexName] << ">>" << nameHistResult << std::endl;
        // histResult = (TH1 *)((TH1 *)tlHistFile->First())->Clone(nameHistResult);
        histResult = (TH1 *) histFirst->Clone(nameHistResult);
        histResult->SetDirectory(tfOutHist);
      } else {
        if (debug) std::cout << "Calculating histogram settings ..." << std::endl;
        Long64_t lowerCorrect;
        Long64_t upperCorrect;
        Int_t nBinCorrect;
        if (vvTypeNameLeafTree[iTree][indexName].SubString("Bool").Start() != -1
            || vvTypeNameLeafTree[iTree][indexName].SubString("bool").Start() != -1) {
          lowerCorrect = 0;
          upperCorrect = 2;
          nBinCorrect = 2;
        } else if (vvTypeNameLeafTree[iTree][indexName].SubString("loat").Start() != -1
                   || vvTypeNameLeafTree[iTree][indexName].SubString("ouble").Start() != -1) {
          Double_t binDensityAverage = 0;
          for (UInt_t iHist = 0; iHist < nHistAllFile; iHist++) {
            binDensityAverage += arrNBinFile[iHist] /
                                 (arrUpperFile[iHist] - arrLowerFile[iHist]) /
                                 nHistAllFile;
          }
          Double_t minimumWholeFile = arrMinimumFile[0];
          Double_t maximumWholeFile = arrMaximumFile[0];
          for (UInt_t iHist = 1; iHist < nHistAllFile; iHist++) {
            if (minimumWholeFile > arrMinimumFile[iHist])
              minimumWholeFile = arrMinimumFile[iHist];
            if (maximumWholeFile < arrMaximumFile[iHist])
              maximumWholeFile = arrMaximumFile[iHist];
          }
          if (debug) std::cout << "minimumWholeFile, maximumWholeFile: " << minimumWholeFile << ", " << maximumWholeFile << std::endl;
          if (TMath::Abs(minimumWholeFile) > 1) {
            Int_t nDigitM1MinimumWholeFile =
                (Int_t)TMath::Floor(TMath::Log10(TMath::Abs(minimumWholeFile)));
            lowerCorrect =
                intpow(10, nDigitM1MinimumWholeFile,
                       TMath::Floor(intpow(10, -nDigitM1MinimumWholeFile,
                                           minimumWholeFile)));
            if (debug) std::cout << "nDigitM1MinimumWholeFile, lowerCorrect:"
                       << nDigitM1MinimumWholeFile << ", " << lowerCorrect << std::endl;
          } else {
            lowerCorrect = TMath::Floor(minimumWholeFile);
          }
          if (TMath::Abs(maximumWholeFile) > 1) {
            Int_t nDigitM1MaximumWholeFile =
                (Int_t)TMath::Floor(TMath::Log10(TMath::Abs(maximumWholeFile)));
            upperCorrect =
                intpow(10, nDigitM1MaximumWholeFile,
                       TMath::Ceil(intpow(10, -nDigitM1MaximumWholeFile,
                                          maximumWholeFile)));
                if (debug) std::cout << "nDigitM1MaximumWholeFile, upperCorrect:"
                       << nDigitM1MaximumWholeFile << ", " << upperCorrect << std::endl;
          } else {
            upperCorrect = TMath::Ceil(maximumWholeFile);
          }
          if (lowerCorrect == upperCorrect) {
            upperCorrect += 1;
            nBinCorrect = TMath::Ceil(binDensityAverage);
          } else {
            nBinCorrect = (upperCorrect - lowerCorrect) * binDensityAverage;
          }
          if (nBinCorrect > 1) {
            Int_t nDigitM1NBinCorrect =
                (Int_t)TMath::Floor(TMath::Log10(nBinCorrect));
            nBinCorrect = 
                intpow((Int_t)10, nDigitM1NBinCorrect,
                       TMath::Nint(intpow((Int_t)10, -nDigitM1NBinCorrect,
                       nBinCorrect)));
          }
          if (nBinCorrect <= 0) {
            nBinCorrect = 1;
          }
        } else {
          Long64_t minimumWholeFile = arrMinimumFile[0];
          Long64_t maximumWholeFile = arrMaximumFile[0];
          for (UInt_t iHist=1; iHist<nHistAllFile; iHist++) {
            if (minimumWholeFile > arrMinimumFile[iHist]) {
              minimumWholeFile = arrMinimumFile[iHist];
            }
            if (maximumWholeFile < arrMaximumFile[iHist]) {
              maximumWholeFile = arrMaximumFile[iHist];
            }
          }
          lowerCorrect = minimumWholeFile;
          upperCorrect = maximumWholeFile + 1;
          nBinCorrect = maximumWholeFile - minimumWholeFile;
        }
        TString tstrSettingHistLeaf = (TString) "(" + nBinCorrect + "," + lowerCorrect + "," + upperCorrect + ")";
        if (debug) std::cout << vvTypeNameLeafTree[iTree][indexName];
        if (debug) std::cout << " (nBinCorrect, lowerCorrect, upperCorrect): " << tstrSettingHistLeaf << std::endl;
        if (debug) std::cout << "Drawing histograms with correct settings ...";
        // TList *tlHistCorrected = new TList;
        std::vector<TString> vNameHistCorrected;
        // UInt_t iFile=0;
        for (UInt_t iHist=0, iFile=0; iHist<nHistAllFile; iHist++, iFile++) {
          while (!vvvIsHistFileLeafTree[iTree][indexName][iFile] && iFile < nFileTot) {
            std::cerr << "(Tree index, Leaf index) (" << iTree << ", " << indexName
            << ") not found in file index " << iFile
            << " (" << vvNameModifiedLeafTree[iTree][indexName] << ")" << std::endl; 
            iFile++;
          }
          TString nameHistCorrected = getNameHistOfFile(iTree, indexName, iFile, "Corrected");
          vNameHistCorrected.push_back(nameHistCorrected);
          arrTT[iTree]->Draw(arrLeafnameFile[iHist] + ">>" + nameHistCorrected + tstrSettingHistLeaf);
          TH1 *histCorrected = (TH1 *) gDirectory->Get(nameHistCorrected);
          histCorrected->SetDirectory(tfCorrectedHist);
          histCorrected->SetName(nameHistCorrected);
          histCorrected->Scale(tvdWeightFile[iFile]);
          // tlHistCorrected->AddLast(histCorrected);
          if (debug) std::cout << "iFile: " << iFile;
        }
        // tfCorrectedHist->Write();
        // tfCorrectedHist->Close();
        // tfCorrectedHist->Delete();
        // delete tfCorrectedHist;
        // tfCorrectedHist = TFile::Open(dirCondorPackCurrent + "/" + "output_" + nameDatagroup + "_" + "Corrected" + "_" + nameClusterID + "_hist.root", "read");
        if (debug) std::cout << " Done." << std::endl;
        if (debug) std::cout << "Merging to histResult (" << nameHistResult << ")...";
        // TObjArray *toaHistCorrected = new TObjArray(vNameHistCorrected.size());
        // for (UInt_t iHist=0; iHist<vNameHistCorrected.size(); iHist++) {
        //   TString nameHistCorrected = vNameHistCorrected[iHist];
        //   (*toaHistCorrected)[iHist] = (TH1 *) tfCorrectedHist->Get(nameHistCorrected);
        // }
        if (vNameHistCorrected.size() == 0) {
          std::cerr << "Fatal: " << "No histograms found for (" << iTree << ", " << indexName 
          << ") (" << vvNameModifiedLeafTree[iTree][indexName] << ")";
          histResult = nullptr;
        } else {
          // histResult = (TH1 *) (*toaHistCorrected)[0]->Clone(nameHistResult);
          if (debug) std::cout << "Picking up " << vNameHistCorrected[0] << " to histResult ...";
          TH1 *histCorrectedFirst = (TH1 *) tfCorrectedHist->Get(vNameHistCorrected[0]);
          if (debug) std::cout << " histCorrectedFirst: " << histCorrectedFirst;
          histResult = (TH1 *)histCorrectedFirst->Clone(nameHistResult);
          if (debug) std::cout << " histResult: " << histResult << " Done. " << std::endl;
          histResult->SetDirectory(tfOutHist);
          histResult->Clear();
          histResult->SetName(nameHistResult);
          // histResult->Merge(toaHistCorrected);
          for (UInt_t iHist=1; iHist<vNameHistCorrected.size(); iHist++) {
            if (debug) std::cout << "Getting " << vNameHistCorrected[iHist] << " ...";
            TH1 *histCorrected = (TH1 *) tfCorrectedHist->Get(vNameHistCorrected[iHist]);
            if (debug) std::cout << " Adding " << histCorrected;
            histResult->Add(histCorrected);
            if (debug) std::cout << " Done." << std::endl;
          }
        }
        if (debug) std::cout << " Done." << std::endl;
      }
      histResult->SetTitle(vvTitleLeafTree[iTree][indexName]);
      // if (debug) std::cout << "Writing histResult ...";
      // histResult->Write();
      // if (debug) std::cout << " Done." << std::endl;
      // histResult->Merge(tlHistFileLeaf);
      // ((TList *)(*toatlResult)[iTree])->AddLast(histResult);
      indexName++;
    }
    if (debug) std::cout << "Writing into tfOutHist ...";
    tfOutHist->Write();
    if (debug) std::cout << " Done." << std::endl;
    tfOutHist->Close();
    tfCorrectedHist->Close();
  }
  tfAutogenHist->Close();
  // TFile *tfOut = TFile::Open(dirCondorPackCurrent + "/" + "output_" + nameDatagroup + "_" + nameClusterID + "_hist.root", toRecreateOutFile ? "recreate" : "update");
  // toatlResult->Write("toatlResult");
  // tfOut->Close();

}
