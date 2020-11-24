#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>

#include "mergeToHists.h"

// #include <TH1F.h>
#include <TLeaf.h>
#include <TList.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVector.h>
#include <TVectorFfwd.h>
// #include <TClonesArray.h>
// #include <TObjArray.h>
#include <TMath.h>
#include <TSystem.h>

#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#ifndef INTPOW_FUNCTION
#define INTPOW_FUNCTION
template <typename TBase, typename TIndex, typename TRatio>
TRatio intpow(const TBase base, const TIndex index, const TRatio ratio) {
  if (index == 0) {
    return ratio;
  }
  TRatio result = ratio;
  if (index > 0) {
    for (TIndex i = 0; i < index; i++) {
      result *= base;
    }
  } else {
    for (TIndex i = 0; i < -index; i++) {
      result /= base;
    }
  }
  return result;
}
#endif

#ifndef TSTRING_UTILS
#define TSTRING_UTILS
Bool_t startsWithOfTString(const TString target, const TString pattern,
                           const Ssiz_t displacement = 0) {
  return target.SubString(pattern).Start() == displacement;
}

// Bool_t endsWithOfTString(const TString target, const TString pattern, const
// Ssiz_t displacement=0) {
//   return target.SubString(pattern).Start() == target.Length() -
//   pattern.Length() - displacement;
// }

Bool_t containsOfTString(const TString target, const TString pattern) {
  return target.SubString(pattern).Start() != -1;
}
#endif

void mergeToHists(const std::vector<TString> vNameTT,
                  std::vector<UInt_t> vNumberFile,
                  std::vector<Double_t> vCrossSection,
                  std::function<TString(UInt_t)> funPathTFIn, TString dirTFTemp,
                  std::function<TString(TString)> funNameTFTemp,
                  TString dirTFOut,
                  std::function<TString(TString)> funNameTFOut,
                  TString seperatorPath,
                  std::function<void(Int_t &NBinCorrect, Double_t &lowerCorrect,
                                     Double_t &upperCorrect, TString nameTT,
                                     TString nameLeafModified,
                                     TString typeNameLeaf, TString titleLeaf)>
                      adjustHistSettingPerLeafTreeExtra,
                  const Bool_t toRecreateOutFile, const Bool_t debug,
                  const Bool_t allowMissing, UInt_t nInfileOpenMax) {
  Bool_t areSomeInfilesClosed = false;
  UInt_t nInfileOpen = 0;
  std::function<TString(TString)> modifyNameLeaf =
      [](TString nameLeaf) -> TString {
    TString nameLeafNew = TString(nameLeaf);
    Ssiz_t indexLeftSquareBracket = nameLeaf.First('[');
    if (indexLeftSquareBracket != -1) {
      nameLeafNew.Resize(indexLeftSquareBracket);
    }
    return nameLeafNew;
  };

  Long64_t nEntryOriginal;

  const UInt_t nTT = vNameTT.size();

  TTree *arrTT[nTT];

  auto mountVarsFromTDir = [&nEntryOriginal, &arrTT, nTT,
                            vNameTT](TDirectory *tdir) {
    nEntryOriginal = (*((TVectorD *)tdir->Get("tvdNEntryOriginal")))[0];
    for (UInt_t i = 0; i < nTT; i++) {
      arrTT[i] = (TTree *)tdir->Get(vNameTT[i]);
    }
  };

  UInt_t nDataset = vNumberFile.size();
  if (nDataset != vCrossSection.size()) {
    std::cerr << "The length of fileNumbers (" << nDataset
              << ") does not match that of the crossSection ("
              << vCrossSection.size() << ")" << std::endl;
  }
  nDataset = TMath::Min(nDataset, (UInt_t)vCrossSection.size());

  UInt_t nFileTot = 0;
  for (UInt_t iDataset = 0; iDataset < nDataset; iDataset++) {
    nFileTot += vNumberFile[iDataset];
  }

  if (debug) std::cout << "nFileTot: " << nFileTot << std::endl;
  TCanvas *c1;

  const UInt_t nFileTotOriginal = nFileTot;
  const std::vector<UInt_t> vNumberFileOriginal = vNumberFile;
  TFile *arrFile[nFileTotOriginal];
  Bool_t arrIsOpenableFile[nFileTotOriginal];  //< Array of whether each file is
                                               // openable or not (The per-file
                                               // vectors cares only about
                                               // openable files).

  auto closeTFilesIn = [&nInfileOpen, &nInfileOpenMax, &arrIsOpenableFile,
                        &arrFile, debug](
                           Int_t iFileOriginalToClose,
                           std::function<void()> funDoBeforeClose = nullptr) {
    // tfAutogenHist->Write();
    // areSomeInfilesClosed = true;
    if (funDoBeforeClose != nullptr) funDoBeforeClose();
    while (nInfileOpen > 0 && iFileOriginalToClose >= 0) {
      if (arrIsOpenableFile[iFileOriginalToClose]) {
        arrFile[iFileOriginalToClose]->Close();
        delete arrFile[iFileOriginalToClose];
        nInfileOpen--;
      }
      iFileOriginalToClose--;
    }
    if (debug && nInfileOpen > 0) {
      std::cerr << "Closing wrong number of files" << std::endl;
    }
  };

  std::vector<Double_t> vWeightDataset(
      nDataset);  //< The weight of each dataset
  vWeightDataset.clear();
  std::vector<Double_t> vWeightFile(
      nFileTot);  //< The weight of each openable file
  vWeightFile.clear();
  std::vector<std::vector<std::vector<Bool_t>>> vvvIsHistFileLeafTree(
      nTT);  //< The existence in each file of each leaf in each tree
  vvvIsHistFileLeafTree.clear();
  // std::vector<std::vector<std::vector<TString>>>
  // vvvNameEffectiveFileLeafTree(nTT); //< The leafname in each file of each
  // (existing) leaf in each tree std::vector<std::vector<std::vector<Int_t>>>
  // vvvNEntryFileLeafTree(nTT); //< The number of entries of the histogram in
  // each file of each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Double_t>>>
  // vvvHistMinFileLeafTree(nTT); //< The minimum value of each histogram in
  // each file of each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Double_t>>>
  // vvvHistMaxFileLeafTree(nTT); //< The maximum value of each histogram in
  // each file of each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Int_t>>> vvvNBinFileLeafTree(nTT); //<
  // The bin-number of the histogram in each file of each (existing) leaf in
  // each tree std::vector<std::vector<std::vector<Double_t>>>
  // vvvLowerFileLeafTree(nTT); //< The lower limit of the histogram in each
  // file in each (existing) leaf in each tree
  // std::vector<std::vector<std::vector<Double_t>>> vvvUpperFileLeafTree(nTT);
  // //< The upper limit of the histogram in each file in each (existing) leaf
  // in each tree
  std::vector<std::vector<TString>> vvNameModifiedLeafTree(
      nTT);  //< The modified leafname of each leaf in each tree
  vvNameModifiedLeafTree.clear();
  std::vector<std::vector<TString>> vvTitleLeafTree(
      nTT);  // The title of the branch of each leaf in each tree
  vvTitleLeafTree.clear();
  std::vector<std::vector<TString>> vvTypeNameLeafTree(
      nTT);  // The type name of each each leaf in each tree
  vvTypeNameLeafTree.clear();
  auto getNameHistOfFile = [&vNameTT, &vvNameModifiedLeafTree](
                               UInt_t iTree, UInt_t indexName, UInt_t iFile,
                               TString suffixStage) -> TString {
    return (TString) "h" + vNameTT[iTree] +
           vvNameModifiedLeafTree[iTree][indexName] + suffixStage + iFile;
  };

  if (dirTFTemp.Length() > 1 && dirTFTemp.EndsWith(seperatorPath)) {
    dirTFTemp.Resize(dirTFTemp.Length() - 1);
  }
  if (!gSystem->AccessPathName(dirTFTemp, kFileExists)) {
    gSystem->mkdir(dirTFTemp);
  }
  TString pathTFAutogenHist =
      (dirTFTemp == "" ? "" : dirTFTemp + seperatorPath) +
      funNameTFTemp("Autogen");
  TString pathTFCorrectedHist =
      (dirTFTemp == "" ? "" : dirTFTemp + seperatorPath) +
      funNameTFTemp("Corrected");
  std::function<TString(UInt_t, UInt_t)> getNameHistResult =
      [&vvNameModifiedLeafTree](UInt_t iTree, UInt_t indexName) -> TString {
    return (TString) "h" + vvNameModifiedLeafTree[iTree][indexName];
  };
  if (dirTFOut.Length() > 1 && dirTFOut.EndsWith(seperatorPath)) {
    dirTFOut.Resize(dirTFOut.Length() - 1);
  }
  std::function<TString(UInt_t)> funPathTFOutHist =
      [dirTFOut, funNameTFOut, seperatorPath,
       &vNameTT](UInt_t iTree) -> TString {
    return dirTFOut + seperatorPath + funNameTFOut(vNameTT[iTree]);
  };
  std::function<void(Int_t &, Double_t &, Double_t &, UInt_t, UInt_t)>
      adjustHistSettingPerLeafTreeExtraIndex =
          [&vNameTT, &vvNameModifiedLeafTree, &vvTypeNameLeafTree,
           &vvTitleLeafTree, &adjustHistSettingPerLeafTreeExtra](
              Int_t &nBinCorrect, Double_t &lowerCorrect,
              Double_t &upperCorrect, UInt_t iTree, UInt_t indexName) {
            TString nameTT = vNameTT[iTree];
            TString nameLeafModified = vvNameModifiedLeafTree[iTree][indexName];
            TString typeName = vvTypeNameLeafTree[iTree][indexName];
            TString title = vvTitleLeafTree[iTree][indexName];
            return adjustHistSettingPerLeafTreeExtra(
                nBinCorrect, lowerCorrect, upperCorrect, nameTT,
                nameLeafModified, typeName, title);
          };
  for (UInt_t iTree = 0; iTree < nTT; iTree++) {
    std::vector<std::vector<Bool_t>> vvIsHistFileLeaf;
    vvIsHistFileLeaf.clear();
    vvvIsHistFileLeafTree.push_back(vvIsHistFileLeaf);
    std::vector<std::vector<TH1 *>> vvHistFileLeaf;
    vvHistFileLeaf.clear();
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
    vNameModifiedLeaf.clear();
    vvNameModifiedLeafTree.push_back(vNameModifiedLeaf);
    std::vector<TString> vTitleLeaf;
    vTitleLeaf.clear();
    vvTitleLeafTree.push_back(vTitleLeaf);
    std::vector<TString> vTypeNameLeaf;
    vTypeNameLeaf.clear();
    vvTypeNameLeafTree.push_back(vTypeNameLeaf);
    // TList *tltlHistFileLeaf = new TList;
    // (*toatltlHistFileLeafTree)[iTree] = tltlHistFileLeaf;
  }
  // UInt_t iFile = 0;
  TFile *tfAutogenHist, *tfCorrectedHist;
  tfCorrectedHist = TFile::Open(pathTFCorrectedHist,
                                toRecreateOutFile ? "recreate" : "update");
  tfCorrectedHist->Close();
  delete tfCorrectedHist;
  tfAutogenHist =
      TFile::Open(pathTFAutogenHist, toRecreateOutFile ? "recreate" : "update");
  // TFile *tfCorrectedHist = TFile::Open(pathTFCorrectedHist, toRecreateOutFile
  // ? "recreate" : "update");
  for (UInt_t iDataset = 0, iFile = 0, nFileMissingTot = 0; iDataset < nDataset;
       iDataset++) {
    Long64_t nEntryOriginalDatasetCurrent = 0;
    if (debug) std::cout << "iDataset: " << iDataset << std::endl;
    Bool_t isFileOpenable = true;
    for (UInt_t jFile = 0, nFileMissingThisDataset = 0;
         jFile < vNumberFile[iDataset];
         (isFileOpenable ? jFile++ : nFileMissingThisDataset++),
                (isFileOpenable ? iFile++ : nFileMissingTot++)) {
      // TString filenameCurrent = filenameGeneralCurrentCluster + "_" + (iFile
      // + nFileMissingTot) + ".root";
      TString pathFileIn = funPathTFIn(iFile + nFileMissingTot);
      if (debug)
        std::cout << "Opening TFile " << pathFileIn << " ..." << std::endl;
      gSystem->ExpandPathName(pathFileIn);
      TFile *tfCurrent = nullptr;
      // gSystem->AccessPathName(pathFileIn,kFileExists) returns FALSE if the
      // file of pathFileIn EXISTS!
      isFileOpenable = !gSystem->AccessPathName(pathFileIn, kFileExists);
      if (isFileOpenable) {
        tfCurrent = TFile::Open(pathFileIn, "READ");
        isFileOpenable = tfCurrent != nullptr && !tfCurrent->IsZombie();
      } else {
        std::cerr << "File not exists!" << std::endl;
      }
      arrIsOpenableFile[iFile + nFileMissingTot] = isFileOpenable;
      arrFile[iFile + nFileMissingTot] = tfCurrent;
      if (!isFileOpenable) {
        std::cerr << "File not openable: " << pathFileIn << " (" << tfCurrent
                  << ")" << std::endl;
        if (allowMissing) {
          vNumberFile[iDataset]--;
          nFileTot--;
          continue;
        }
      }
      nInfileOpen++;
      if (debug) std::cout << "Done." << std::endl;
      if (debug)
        std::cout << "Mounting variables from the file ..." << std::endl;
      mountVarsFromTDir(tfCurrent);
      if (debug) std::cout << "Done." << std::endl;
      if (debug) std::cout << "nEntryOriginal: " << nEntryOriginal << std::endl;
      nEntryOriginalDatasetCurrent += nEntryOriginal;
      for (UInt_t iTree = 0; iTree < nTT; iTree++) {
        if (debug) std::cout << "iTree: " << iTree << std::endl;
        if (debug)
          std::cout << "Getting leafnames from the tree ..." << std::endl;
        // TList *tltlHistFileLeaf = (TList *)
        // (*toatltlHistFileLeafTree)[iTree]; std::vector<TString>
        // vNameModifiedLeafOld = vvNameModifiedLeafTree[iTree]; const UInt_t
        // nNameModifiedLeafOld = vNameModifiedLeafOld.size();
        const UInt_t nNameModifiedLeafOld =
            vvNameModifiedLeafTree[iTree].size();
        Bool_t arrIsHistLeafOld[nNameModifiedLeafOld];
        for (UInt_t indexNameOld = 0; indexNameOld < nNameModifiedLeafOld;
             indexNameOld++) {
          arrIsHistLeafOld[indexNameOld] = false;
        }  //< For possible missing file and changing nFileTot
        for (TObject *leafRaw : *(arrTT[iTree]->GetListOfLeaves())) {
          TLeaf *leaf = (TLeaf *)leafRaw;
          TString nameLeaf = leaf->GetName();
          TString nameLeafModified = modifyNameLeaf(nameLeaf);
          Bool_t isIncluded = false;
          UInt_t indexName = 0;
          // TList *tlHistFileCurrent;
          for (TString nameIncluded : vvNameModifiedLeafTree[iTree]) {
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
            if (debug)
              std::cout << "nameLeaf: " << nameLeaf
                        << ", nameLeafModified: " << nameLeafModified
                        << std::endl;
            vvNameModifiedLeafTree[iTree].push_back(nameLeafModified);
            vvTitleLeafTree[iTree].push_back(leaf->GetBranch()->GetTitle());
            vvTypeNameLeafTree[iTree].push_back(leaf->GetTypeName());
            std::vector<Bool_t> vIsHistFile(nFileTot);
            vIsHistFile.clear();
            std::vector<TH1 *> vHistFile(nFileTot);
            vHistFile.clear();
            // Commented out for possible missing file and changing nFileTot
            // for (UInt_t i=0; i<nFileTot; i++) {
            //   vIsHistFile.push_back(false);
            // }
            // Added for possible missing file and changing nFileTot
            for (UInt_t i = 0; i < iFile; i++) {
              vIsHistFile.push_back(false);
            }
            vIsHistFile.push_back(true);
            // Now vIsHistFile has iFile+1 elements
            vvvIsHistFileLeafTree[iTree].push_back(vIsHistFile);
            // vvvHistFileLeafTree[iTree].push_back(vHistFile);
            // TList *tlHistFileCurrent = new TList;
            // tltlHistFileLeaf->AddLast(tlHistFileCurrent);
          } else {
            // indexName should be smaller than current
            // vNameModifiedLeafTree[iTree].size() unless the nameLeafModified
            // is duplicated
            if (indexName >= vvNameModifiedLeafTree[iTree].size()) {
              std::cerr << "nameLeafModified duplication (" << nameLeaf << ", "
                        << nameLeafModified << ")" << std::endl;
              std::cerr << " indexName of the original leaf: " << indexName
                        << std::endl;
            } else if (arrIsHistLeafOld[indexName]) {
              std::cerr << "nameLeafModified duplication (" << nameLeaf << ", "
                        << nameLeafModified << ")" << std::endl;
              std::cerr << "indexName of the original leaf: " << indexName
                        << std::endl;
            } else {
              arrIsHistLeafOld[indexName] = true;
            }
          }
          // vvvIsHistFileLeafTree[iTree][indexName][iFile] = true;
          TString nameHistAutogen =
              getNameHistOfFile(iTree, indexName, iFile, "Autogen");
          if (debug)
            std::cout << "Drawing " << nameLeaf << " as " << nameHistAutogen
                      << " ...";
          arrTT[iTree]->Draw(nameLeaf + ">>" + nameHistAutogen);
          TH1 *histAutogen = (TH1 *)gDirectory->Get(nameHistAutogen);
          histAutogen->SetTitle(nameLeaf);
          histAutogen->SetName(nameHistAutogen);
          histAutogen->SetDirectory(tfAutogenHist);
          // histAutogen->Write(nameHistAutogen);
          // vvvHistFileLeafTree[iTree][indexName].push_back(histAutogen);
          if (debug && histAutogen == nullptr)
            std::cerr << "Fatal: " << nameHistAutogen << " is nullptr!"
                      << std::endl;
          if (debug) std::cout << " " << histAutogen;
          if (debug) std::cout << " Done." << std::endl;
        }
        // As for now, vvvIsHistFileLeafTree[iTree][indexName]
        // has iFile elements if indexName < nNameModifiedLeafOld
        // and iFile + 1 elements if indexName >= nNameModifiedLeafOld
        // Now we are going to push_back one element to every
        // vvvIsHistFileLeafTree[iTree][indexName]
        // that indexName < nNameModifiedLeafOld
        // according to arrIsHistLeafOld
        for (UInt_t indexNameOld = 0; indexNameOld < nNameModifiedLeafOld;
             indexNameOld++) {
          vvvIsHistFileLeafTree[iTree][indexNameOld].push_back(
              arrIsHistLeafOld[indexNameOld]);
        }
        if (debug) {
          for (UInt_t indexName = 0;
               indexName < vvvIsHistFileLeafTree[iTree].size(); indexName++) {
            if (vvvIsHistFileLeafTree[iTree][indexName].size() != iFile + 1) {
              std::cerr << "vvvIsHistFileLeafTree[" << iTree << "]["
                        << indexName << "].size() != iFile + 1"
                        << "("
                        << "size: "
                        << vvvIsHistFileLeafTree[iTree][indexName].size()
                        << ", iFile: " << iFile << ")" << std::endl;
            }
          }
        }
      }
      if (debug) std::cout << "Done getting names." << std::endl;
      if (nInfileOpen >= nInfileOpenMax) {
        closeTFilesIn(iFile + nFileMissingTot,
                      [&areSomeInfilesClosed, &tfAutogenHist]() {
                        areSomeInfilesClosed = true;
                        // tfAutogenHist->Write();
                      });
      }
    }
    if (debug)
      std::cout << "nEntryOriginalDatasetCurrent: "
                << nEntryOriginalDatasetCurrent << std::endl;
    if (debug)
      std::cout << "vCrossSection[" << iDataset
                << "]: " << vCrossSection[iDataset] << std::endl;
    Double_t weight =
        (nEntryOriginalDatasetCurrent > 0)
            ? vCrossSection[iDataset] / nEntryOriginalDatasetCurrent
            : 0;
    vWeightDataset.push_back(weight);
    if (debug)
      std::cout << "vWeightDataset[" << iDataset
                << "]: " << vWeightDataset[iDataset] << "," << std::endl;
    if (debug) std::cout << "composed of files index";
    for (UInt_t jFile = 0; jFile < vNumberFile[iDataset]; jFile++) {
      UInt_t indexFile = iFile - vNumberFile[iDataset] + jFile;
      if (debug) std::cout << " " << indexFile;
      vWeightFile.push_back(weight);
    }
    if (debug) std::cout << "." << std::endl;
    if (debug && vWeightDataset.size() != iDataset + 1)
      std::cerr << "Error: vWeightDataset.size() != iDataset + 1 "
                << " (size:" << vWeightDataset.size()
                << ", iDataset:" << iDataset + 1 << ")" << std::endl;
  }
  tfAutogenHist->Write();
  // tfAutogenHist->Close();
  // delete tfAutogenHist;

  // Close all the input files if and only if the total number of opened input files are greater than nInfileOpenMax
  if (areSomeInfilesClosed) {
    closeTFilesIn(nFileTotOriginal - 1);
  }

  // tfAutogenHist = TFile::Open(dirCondorPackCurrent + seperatorPath +
  // "output_" + nameDatagroup + "_" + "Autogen" + "_" + nameClusterID +
  // "_hist.root", "read");
  gStyle->SetOptStat(111111);
  // TObjArray *toatlResult = new TObjArray(nTT);
  std::vector<std::vector<Bool_t>> vvIsAllEmptyLeafTree(
      nTT);  //< Whether all the histograms have zero entry of each leaf in each
             // tree
  vvIsAllEmptyLeafTree.clear();
  std::vector<std::vector<TString>> vvHistsettingLeafTree(
      nTT);  //< The histogram settings "(binnumber, lower, upper)" of each leaf
             // in each tree;
  vvHistsettingLeafTree.clear();
  std::vector<std::vector<std::vector<TString>>> vvvNameFileLeafTree(nTT);
  vvvNameFileLeafTree.clear();
  for (UInt_t iTree = 0; iTree < nTT; iTree++) {
    {
      std::vector<Bool_t> vIsAllEmptyLeaf;
      vIsAllEmptyLeaf.clear();
      vvIsAllEmptyLeafTree.push_back(vIsAllEmptyLeaf);
      std::vector<TString> vHistsettingLeaf;
      vHistsettingLeaf.clear();
      vvHistsettingLeafTree.push_back(vHistsettingLeaf);
      std::vector<std::vector<TString>> vvNameFileLeaf;
      vvNameFileLeaf.clear();
      vvvNameFileLeafTree.push_back(vvNameFileLeaf);
    }
    // (*toatlResult)[iTree] = new TList;
    // UInt_t indexName = 0;
    // for (const auto&& tlHistFileRaw: *((TList *)
    // (*toatltlHistFileLeafTree)[iTree])) {
    UInt_t nLeaf = vvNameModifiedLeafTree[iTree].size();
    for (UInt_t indexName = 0; indexName < nLeaf; indexName++) {
      if (debug) std::cout << "indexName: " << indexName << std::endl;
      // TList *tlHistFile = (TList *) tlHistFileRaw;
      // TH1 *histFirst = (TH1 *) tlHistFile->First();

      // Index of the first true element of
      // vvvIsHistFileLeafTree[iTree][indexName]
      UInt_t iFileFirst = std::distance(
          vvvIsHistFileLeafTree[iTree][indexName].begin(),
          std::find(vvvIsHistFileLeafTree[iTree][indexName].begin(),
                    vvvIsHistFileLeafTree[iTree][indexName].end(), true));
      // tfAutogenHist = TFile::Open(pathTFAutogenHist, "read");
      tfAutogenHist->ReOpen("read");
      // tfCorrectedHist = TFile::Open(pathTFCorrectedHist, toRecreateOutFile ?
      // "recreate" : "update");
      // TH1 *histFirst = (TH1 *) gDirectory->Get(getNameHistOfFile(iTree,
      // indexName, iFileFirst, "Autogen"));
      TH1 *histFirst = (TH1 *)tfAutogenHist->Get(
          getNameHistOfFile(iTree, indexName, iFileFirst, "Autogen"));
      // TH1 *histFirst = (TH1 *) vvvHistFileLeafTree[iTree][indexName][0];
      // TH1 *histResult = (TH1 *) histFirst->Clone("h" +
      // vvNameModifiedLeafTree[iTree][indexName]); histResult->Clear();
      TH1 *histResult;
      TString nameHistResult = "h" + vvNameModifiedLeafTree[iTree][indexName];
      // UInt_t iDataset = 0;
      // UInt_t iFile = 0;
      // UInt_t nHistAllFile = tlHistFile->GetEntries();
      UInt_t nHistAllFile =
          std::count(vvvIsHistFileLeafTree[iTree][indexName].begin(),
                     vvvIsHistFileLeafTree[iTree][indexName].end(), true);
      if (debug) std::cout << "nHistAllFile: " << nHistAllFile << std::endl;
      Int_t arrNEntryFile[nHistAllFile];
      Double_t arrMinimumFile[nHistAllFile];
      Double_t arrMaximumFile[nHistAllFile];
      Int_t arrNBinFile[nHistAllFile];
      Double_t arrLowerFile[nHistAllFile];
      Double_t arrUpperFile[nFileTot];
      Bool_t isAllEmpty = true;
      // UInt_t iHist = 0;
      if (debug)
        std::cout << "Analyzing autogen histograms for "
                  << vvNameModifiedLeafTree[iTree][indexName] << "...";
      // for (auto histFileRaw: *tlHistFile) {
      std::vector<TString> vLeafnameFile;
      vLeafnameFile.clear();
      for (UInt_t iFile = 0, iHist = 0; iFile < nFileTot; iFile++) {
        if (!vvvIsHistFileLeafTree[iTree][indexName][iFile]) {
          vLeafnameFile.push_back("");
          if (debug) std::cout << " Skipping iFile: " << iFile;
          continue;
        }
        if (debug) std::cout << " iHist: " << iHist;
        // TH1 *histFile = (TH1 *)histFileRaw;
        TString nameHistAutogen =
            getNameHistOfFile(iTree, indexName, iFile, "Autogen");
        // TH1 *histAutogen = (TH1 *) gDirectory->Get(nameHistAutogen);
        TH1 *histAutogen = (TH1 *)tfAutogenHist->Get(nameHistAutogen);
        // TH1 *histAutogen = vvvHistFileLeafTree[iTree][indexName][iHist];
        if (debug && histAutogen == nullptr)
          std::cerr << "Fatal: histAutogen from " << nameHistAutogen
                    << " is nullptr!" << std::endl;
        // histFile->Scale(vWeightFile[iFile]);
        Int_t nEntry = histAutogen->GetEntries();
        if (nEntry > 0) {
          isAllEmpty = false;
        }
        arrNEntryFile[iHist] = nEntry;
        arrMinimumFile[iHist] =
            histAutogen->GetBinCenter(histAutogen->FindFirstBinAbove()) -
            histAutogen->GetBinWidth(1) / 2;
        arrMaximumFile[iHist] =
            histAutogen->GetBinCenter(histAutogen->FindLastBinAbove()) +
            histAutogen->GetBinWidth(1) / 2;
        Int_t nBin = histAutogen->GetNbinsX();
        arrNBinFile[iHist] = nBin;
        arrLowerFile[iHist] =
            histAutogen->GetBinCenter(1) - histAutogen->GetBinWidth(1) / 2;
        arrUpperFile[iHist] =
            histAutogen->GetBinCenter(nBin) + histAutogen->GetBinWidth(1) / 2;
        vLeafnameFile.push_back(histAutogen->GetTitle());
        iHist++;
      }
      vvvNameFileLeafTree[iTree].push_back(vLeafnameFile);
      if (debug) std::cout << " Done." << std::endl;
      vvIsAllEmptyLeafTree[iTree].push_back(isAllEmpty);
      if (isAllEmpty) {
        if (debug)
          std::cout << "Empty leaf encountered: "
                    << vvNameModifiedLeafTree[iTree][indexName] << ">>"
                    << nameHistResult << std::endl;
        // histResult = (TH1 *)((TH1
        // *)tlHistFile->First())->Clone(nameHistResult); histResult = (TH1 *)
        // histFirst->Clone(nameHistResult); tfOutHist =
        // TFile::Open(pathTFOutHist, "update");
        // histResult->SetDirectory(tfOutHist);
        // tfAutogenHist->Close();
        // delete tfAutogenHist;
        vvHistsettingLeafTree[iTree].push_back("");
      } else {
        if (debug)
          std::cout << "Calculating histogram settings ..." << std::endl;
        Double_t lowerCorrect;
        Double_t upperCorrect;
        Int_t nBinCorrect;
        if (containsOfTString(vvTypeNameLeafTree[iTree][indexName], "Bool") ||
            containsOfTString(vvTypeNameLeafTree[iTree][indexName], "bool")) {
          lowerCorrect = 0;
          upperCorrect = 2;
          nBinCorrect = 2;
        } else if (containsOfTString(vvTypeNameLeafTree[iTree][indexName],
                                     "loat") ||
                   containsOfTString(vvTypeNameLeafTree[iTree][indexName],
                                     "ouble")) {
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
          if (debug)
            std::cout << "minimumWholeFile, maximumWholeFile: "
                      << minimumWholeFile << ", " << maximumWholeFile
                      << std::endl;
          if (TMath::Abs(minimumWholeFile) > 1) {
            Int_t nDigitM1MinimumWholeFile =
                (Int_t)TMath::Floor(TMath::Log10(TMath::Abs(minimumWholeFile)));
            lowerCorrect =
                intpow(10, nDigitM1MinimumWholeFile,
                       TMath::Floor(intpow(10, -nDigitM1MinimumWholeFile,
                                           minimumWholeFile)));
            if (debug)
              std::cout << "nDigitM1MinimumWholeFile, lowerCorrect:"
                        << nDigitM1MinimumWholeFile << ", " << lowerCorrect
                        << std::endl;
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
            if (debug)
              std::cout << "nDigitM1MaximumWholeFile, upperCorrect:"
                        << nDigitM1MaximumWholeFile << ", " << upperCorrect
                        << std::endl;
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
          for (UInt_t iHist = 1; iHist < nHistAllFile; iHist++) {
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
        adjustHistSettingPerLeafTreeExtraIndex(nBinCorrect, lowerCorrect,
                                               upperCorrect, iTree, indexName);
        TString tstrSettingHistLeaf = (TString) "(" + nBinCorrect + "," +
                                      lowerCorrect + "," + upperCorrect + ")";
        vvHistsettingLeafTree[iTree].push_back(tstrSettingHistLeaf);
        if (debug) std::cout << vvTypeNameLeafTree[iTree][indexName];
        if (debug)
          std::cout << " (nBinCorrect, lowerCorrect, upperCorrect): "
                    << tstrSettingHistLeaf << std::endl;
        // tfCorrectedHist = TFile::Open(pathTFCorrectedHist, "read");
        if (debug) std::cout << " Done." << std::endl;
        // if (debug) std::cout << "Merging to histResult (" << nameHistResult
        // << ")..."; TObjArray *toaHistCorrected = new
        // TObjArray(vNameHistCorrected.size()); for (UInt_t iHist=0;
        // iHist<vNameHistCorrected.size(); iHist++) {
        //   TString nameHistCorrected = vNameHistCorrected[iHist];
        //   (*toaHistCorrected)[iHist] = (TH1 *)
        //   tfCorrectedHist->Get(nameHistCorrected);
        // tfCorrectedHist->Close();
        if (debug) std::cout << " Done." << std::endl;
      }
      // indexName++; // LO THE BUG that cause half of the histogram to
      // disappear
    }
  }
  if (debug) std::cout << "Drawing histograms with correct settings ..." << std::endl;
  for (UInt_t iFile = 0, iFileOriginal = 0; iFile < nFileTot || iFileOriginal < nFileTotOriginal;
       iFile++, iFileOriginal++) {
    if (debug) std::cout << "iFile: " << iFile;
    TFile *tfInCurrent = nullptr;
    while (iFileOriginal < nFileTotOriginal && !arrIsOpenableFile[iFileOriginal]) {
      iFileOriginal++;
    }
    if (iFileOriginal >= nFileTotOriginal) {
      break;
    }
    if (debug) std::cout << " iFileOriginal: " << iFileOriginal << std::endl;
    if (areSomeInfilesClosed) {
      if (debug) std::cout << "Opening input rootfiles ...";
      arrFile[iFileOriginal] = TFile::Open(funPathTFIn(iFileOriginal));
      nInfileOpen++;
      if (debug) std::cout << " Done." << std::endl;
    }
    tfInCurrent = arrFile[iFileOriginal];
    if (debug) std::cout << "tfInCurrent" << tfInCurrent << std::endl;
    if (debug) std::cout << "Mounting variables ...";
    mountVarsFromTDir(tfInCurrent);
    if (debug) std::cout << " Done." << std::endl;
    tfCorrectedHist = TFile::Open(pathTFCorrectedHist, "update");
    for (UInt_t iTree = 0; iTree < nTT; iTree++) {
      UInt_t nNameModifiedCurrentTree = vvNameModifiedLeafTree[iTree].size();
      for (UInt_t indexName = 0; indexName < nNameModifiedCurrentTree;
           indexName++) {
        // TList *tlHistCorrected = new TList;
        // std::vector<TString> vNameHistCorrected;
        // vNameHistCorrected.clear();
        // UInt_t iFile=0;
        if (!vvvIsHistFileLeafTree[iTree][indexName][iFile]) {
          std::cerr << "(Tree index, Leaf index) (" << iTree << ", "
                    << indexName << ") not found in file index " << iFile
                    << " (" << vvNameModifiedLeafTree[iTree][indexName] << ")"
                    << std::endl;
          continue;
        }
        TString nameHistCorrected =
            getNameHistOfFile(iTree, indexName, iFile, "Corrected");
        arrTT[iTree]->Draw(vvvNameFileLeafTree[iTree][indexName][iFile] + ">>" +
                           nameHistCorrected +
                           vvHistsettingLeafTree[iTree][indexName]);
        TH1 *histCorrected = (TH1 *)gDirectory->Get(nameHistCorrected);
        histCorrected->SetDirectory(tfCorrectedHist);
        histCorrected->SetName(nameHistCorrected);
        histCorrected->Scale(vWeightFile[iFile]);
        histCorrected->Write(nameHistCorrected);
      }
    }
    if (areSomeInfilesClosed && nInfileOpen >= nInfileOpenMax) {
      closeTFilesIn(iFileOriginal);
    }
    tfCorrectedHist->Close();
  }
  if (debug) std::cout << "Closing openable input ROOT files ..." << std::endl;
  closeTFilesIn(nFileTotOriginal - 1);
  if (debug) std::cout << "Done." << std::endl;
  delete tfCorrectedHist;
  // tfAutogenHist->Close();
  // delete tfAutogenHist;
  tfAutogenHist->ReOpen("read");
  tfCorrectedHist = TFile::Open(pathTFCorrectedHist, "read");
  // tfCorrectedHist->ReOpen("read");

  // if (debug) std::cout << "Closing openable input ROOT files ..." <<
  // std::endl; for (UInt_t iFileOriginal = 0; iFileOriginal < nFileTotOriginal;
  //      iFileOriginal++) {
  //   if (arrIsOpenableFile[iFileOriginal]) {
  //     if (debug) std::cout << "Closing iFileOriginal: " << iFileOriginal << "
  //     "; arrFile[iFileOriginal]->Close();
  //   }
  // }
  // if (debug) std::cout << "Done." << std::endl;

  if (debug) std::cout << "Merging histograms ..." << std::endl;
  // tfAutogenHist = TFile::Open(pathTFAutogenHist, "read");
  // tfCorrectedHist = TFile::Open(pathTFCorrectedHist, "read");
  for (UInt_t iTree = 0; iTree < nTT; iTree++) {
    if (debug)
      std::cout << "iTree: " << iTree << " (" << vNameTT[iTree] << ")"
                << std::endl;
    TFile *tfOutHist = TFile::Open(funPathTFOutHist(iTree),
                                   toRecreateOutFile ? "recreate" : "update");
    for (UInt_t indexName = 0; indexName < vvNameModifiedLeafTree[iTree].size();
         indexName++) {
      TH1 *histResult;
      TString nameHistResult = getNameHistResult(iTree, indexName);
      if (debug)
        std::cout << "Generating " << nameHistResult << " ..." << std::endl;
      if (vvIsAllEmptyLeafTree[iTree][indexName]) {
        UInt_t iFileFirst;
        for (iFileFirst = 0;
             !vvvIsHistFileLeafTree[iTree][indexName][iFileFirst]; iFileFirst++)
          ;
        histResult = (TH1 *)tfAutogenHist
                         ->Get(getNameHistOfFile(iTree, indexName, iFileFirst,
                                                 "Autogen"))
                         ->Clone(nameHistResult);
        histResult->SetName(nameHistResult);
        if (debug)
          std::cout << "Got " << nameHistResult << " (empty histogram)."
                    << std::endl;
      } else {
        Bool_t isErrorUnexpected = false;
        TList *tlHistCorrectedToMerge = new TList;
        UInt_t nHist = 0;
        for (UInt_t iFile = 0; iFile < nFileTot; iFile++) {
          if (!vvvIsHistFileLeafTree[iTree][indexName][iFile]) {
            continue;
          }
          if (debug) std::cout << "iFile: " << iFile << " ";
          TString nameHistCorrected =
              getNameHistOfFile(iTree, indexName, iFile, "Corrected");
          TH1 *histCorrected = (TH1 *)tfCorrectedHist->Get(nameHistCorrected);
          if (histCorrected == nullptr) {
            std::cerr << "Fatal: Histogram not found: " << nameHistCorrected
                      << std::endl;
            isErrorUnexpected = true;
            break;
          }
          tlHistCorrectedToMerge->Add(histCorrected->Clone(nameHistCorrected));
          nHist++;
        }
        if (!nHist) {
          std::cerr << "Fatal: "
                    << "No histograms found for (" << iTree << ", " << indexName
                    << ") (" << vvNameModifiedLeafTree[iTree][indexName] << ")"
                    << std::endl;
          histResult = nullptr;
          isErrorUnexpected = true;
        } else {
          if (debug) std::cout << "Merging ...";
          histResult =
              (TH1 *)tlHistCorrectedToMerge->First()->Clone(nameHistResult);
          histResult->SetName(nameHistResult);
          histResult->Clear();
          histResult->Merge(tlHistCorrectedToMerge);
        }
        // if (isErrorUnexpected) continue;
        histResult->SetTitle(vvTitleLeafTree[iTree][indexName]);
        histResult->SetDirectory(tfOutHist);
        histResult->Write(nameHistResult);
        std::cout << " Done." << std::endl;
      }
    }
    // tfOutHist->Write();
    tfOutHist->Close();
  }
  tfCorrectedHist->Close();
  tfAutogenHist->Close();
}
