#include <Rtypes.h>
#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TString.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>

#include "HistMerger.C"

#ifndef TSTRING_UTILS
#define TSTRING_UTILS
Bool_t startsWithOfTString(const TString target, const TString pattern, const Ssiz_t displacement=0) {
  return target.SubString(pattern).Start() == displacement;
}

// Bool_t endsWithOfTString(const TString target, const TString pattern, const Ssiz_t displacement=0) {
//   return target.SubString(pattern).Start() == target.Length() - pattern.Length() - displacement;
// }

Bool_t containsOfTString(const TString target, const TString pattern) {
  return target.SubString(pattern).Start() != -1;
}
#endif

void xAna_monoZ_genhist(const TString nameCondorPack, const TString nameDatagroup, const TString nameClusterID,
                     TString dirCondorPacks = ".", const Bool_t toRecreateOutFile = true,
                     const Bool_t debug = false, const Bool_t allowMissing = false) {

  const TString seperatorPath = "/";
  const TString namesLepton[] = {"Electron", "Muon",
                                 "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon",
                                      "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};

  std::vector<TString> vNameTT;
  vNameTT.clear();
  for (UShort_t i=0; i<2; i++) {
    TString nameLeptonCurrent = namesLepton[i];
    vNameTT.push_back((TString) "ZMassCutted" + nameLeptonCurrent);
  }

  std::vector<UInt_t> vNumberFile;
  vNumberFile.clear();
  std::vector<Double_t> vCrossSection;
  vCrossSection.clear();

  TString dirCondorPackCurrent = dirCondorPacks + seperatorPath + nameCondorPack;
  TString pathTextFileNumbers = dirCondorPackCurrent + seperatorPath + "fileNumbers_" + nameDatagroup + ".txt";
  TString pathTextCrossSections = dirCondorPackCurrent + seperatorPath + "crossSections_" + nameDatagroup + ".txt";

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
    //TODO
    if (nameCondorPack == "preselect") {
      if (nameDatagroup == "TT") {
        vNumberFile = {
#include "preselect/fileNumbers_TT_withComma.txt"
      };
        vCrossSection = {730.6, 730.6};
      } else if (nameDatagroup == "DYJets") {
        vNumberFile = {
#include "preselect/fileNumbers_DYJets_withComma.txt"
        };
        vCrossSection = {
#include "preselect/crossSections_DYJets_withComma.txt"
        };
      } else if (nameDatagroup == "signal_Mx2-150_Mv-500_Mx1-1_ctau-1") {
        vNumberFile = {4};
        vCrossSection = {1};
      }
    }
  }

  std::function<TString (TString)> funNameTFTemp = [nameDatagroup, nameClusterID](TString keyword)->TString {
      return "output_" + nameDatagroup + "_" + keyword + "_" + nameClusterID + "_hist.root";
    };

  TString dirTFInAndSlash = dirCondorPackCurrent + seperatorPath + "dataTest_" + nameDatagroup + seperatorPath;
  
  std::function<TString (UInt_t)> funPathTFIn =  [dirTFInAndSlash, nameDatagroup, nameClusterID](UInt_t iFile) {
      return dirTFInAndSlash + "output_" + nameDatagroup + "_" + nameClusterID + "_" + iFile + ".root";
    };
  std::function<TString (TString)> funNameTFOut = [nameDatagroup, nameClusterID] (TString nameTT) -> TString {
      return "output_" + nameDatagroup + "_" + nameTT + "_" + nameClusterID + "_hist.root";
    };

  std::function<void (Int_t&, Double_t&, Double_t&, TString, TString, TString, TString)> adjustHistSetting = [debug](Int_t& nBinCorrect, Double_t& lowerCorrect, Double_t& upperCorrect, TString nameTT, TString nameLeafModified, TString typeNameLeaf, TString titleLeaf) {

    if (containsOfTString(typeNameLeaf, (TString)"loat") || containsOfTString(typeNameLeaf, (TString)"ouble")) {
      if (nameLeafModified.EndsWith("_phi")) {
        if (debug) std::cout << "Found _phi" << std::endl;
        lowerCorrect = -TMath::Pi();
        upperCorrect = TMath::Pi();
        return;
      }
      if (startsWithOfTString(nameLeafModified, (TString)"Jet_btag")) {
        if (startsWithOfTString(nameLeafModified, (TString)"Jet_btagCMVA")) {
          if (debug) std::cout << "Found Jet_btagCMVA" << std::endl;
          lowerCorrect = -1;
          upperCorrect = 1;
        } else {
          if (debug) std::cout << "Found Jet_btag" << std::endl;
          lowerCorrect = 0;
          upperCorrect = 1;
        }
        return;
      }
      if (startsWithOfTString(nameLeafModified, (TString)"Jet_ch")) {
        if (debug) std::cout << "Found Jet_ch" << std::endl;
        lowerCorrect = 0;
        upperCorrect = 1;
        return;
      }
      if (startsWithOfTString(nameLeafModified, (TString)"Jet_qgl")) {
        if (debug) std::cout << "Found Jet_qgl" << std::endl;
        lowerCorrect = 0;
        upperCorrect = 1;
        return;
      }
      if (titleLeaf.EndsWith("Energy Fraction") || titleLeaf.EndsWith("energy fraction")) {
        if (debug) std::cout << "Found Energy Fraction" << std::endl;
        lowerCorrect = 0;
        upperCorrect = 0;
        return;
      }
    }
    if (debug) std::cout << "No additional settings for " << nameLeafModified << " to apply." << std::endl;
  };
  // mergeToHists(vNameTT, vNumberFile, vCrossSection,
  // funPathTFIn,
  // dirCondorPackCurrent, funNameTFTemp,
  // dirCondorPackCurrent, funNameTFOut,
  // seperatorPath,
  // adjustHistSetting, nullptr, nullptr,
  // toRecreateOutFile, debug, allowMissing);
  HistMerger *merger = new HistMerger;
  merger->vNameTT = vNameTT;
  merger->vNumberFile = vNumberFile;
  merger->vCrossSection = vCrossSection;
  merger->funPathTFIn = funPathTFIn;
  merger->dirTFTemp = dirCondorPackCurrent;
  merger->funNameTFTemp = funNameTFTemp;
  merger->dirTFOut = dirCondorPackCurrent;
  merger->funNameTFOut = funNameTFOut;
  merger->seperatorPath = seperatorPath;
  merger->adjustHistSettingPerLeafTreeExtra = adjustHistSetting;
  merger->toRecreateOutFile = toRecreateOutFile;
  merger->debug = debug;
  merger->allowMissing = allowMissing;
  merger->nLeavesToUseCorrectedTempFileMin = 0; // Merging error (limit inconsistent when using tfCorrectedHist)
  merger->funIsToVetoLeaf = [](TString nameTT, TString nameLeafModified)->Bool_t{
    if (containsOfTString(nameLeafModified, "jEntr")) {
      return true;
    }
    return false;
  };
  merger->funTitleLeaf = [](TLeaf *leaf)->TString{ return leaf->GetBranch()->GetTitle(); };
  merger->Run();
}
