#include <Rtypes.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>
#include <TCut.h>

#include <ROOT/RDataFrame.hxx>

#include <vector>
#include <algorithm>
#include <functional>

class AllignedDrawer {
 public:
  std::vector<TString> vNameTT;
  std::vector<UInt_t> vNumberFile;
  std::vector<UInt_t> vCrossSection;
  std::function<Int_t (Long64_t &nEntriesOriginal, TDirectory *tdirIn, TString nameTT)> funRefgetNEntriesOriginal;
  std::function<TString (UInt_t iFile)> funPathTFIn;
  std::function<TString (UInt_t iFile)> funPathTFOut;
  std::function<void(Int_t& nBinsCorrect, Double_t& lowerCorrect,
                     Double_t& upperCorrect, TString nameTT,
                     TString nameLeafModified, TString typeNameLeaf,
                     TString titleLeaf)>
      adjustHistSettingPerLeafTreeExtra;
  Bool_t debug;

  enum class HistTypeFlavor { FloatingPoint, Integer, Boolean };

  struct HistDescription {
    TString expression;
    TCut selection;
    Option_t *optionDraw;
    Bool_t allowMissing;
    Long64_t nEntriesToDrawMax;
    Long64_t firstEntryToDraw;
    Bool_t isPureLeaf;
  };

 protected:
  UInt_t nDataset;
  UInt_t nFilesTotOriginal;

 public:
  AllignedDrawer() {
    vNameTT = {};
    vNumberFile = {};
    vCrossSection = {};
    funRefgetNEntriesOriginal = nullptr;
    funPathTFIn = nullptr;
    funPathTFOut = nullptr;
    debug = false;
    nFilesTotOriginal = std::accumulate(vNumberFile.begin(), vNumberFile.end(), 0u);
    nDataset = vNumberFile.size();
  }

  void Run () {
    for (const auto nameTT: vNameTT) {
      std::vector<Bool_t> vIsFileAvailable(nFilesTotOriginal, true);
      std::vector<Long64_t> vNEntriesOriginal(nFilesTotOriginal, 0);
      std::vector<Long64_t> vNEntriesOriginalDataset(nDataset, 0);

      for (UInt_t iDataset = 0, iFile = 0; iDataset < nDataset; iDataset++) {
        for (UInt_t jFile = 0; jFile < vNumberFile[iDataset]; ++jFile, ++iFile) {
          auto fileIn = TFile::Open(funPathTFIn(iFile));
          auto keyTree = fileIn->GetKey(nameTT);
          if (keyTree == nullptr || keyTree->IsZombie()) {
            vIsFileAvailable[iFile] = false;
          }
          if (funRefgetNEntriesOriginal(vNEntriesOriginal[iFile], fileIn, nameTT)) {
            vIsFileAvailable[iFile] = false;
            vNEntriesOriginal[iFile] = 0;
          }
          fileIn->Close();
        }
      }

      //
    }
  }
};
