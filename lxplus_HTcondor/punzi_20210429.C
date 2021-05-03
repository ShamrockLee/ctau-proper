#include <Rtypes.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TError.h>

#include <vector>
#include <array>
#include <iostream>
#include <utility>

void punzi_20210429(const Bool_t recreate=true, const Bool_t debug=false) {
  const TString namesLepton[] = {"Electron", "Muon", "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon", "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};

  const TString pathCondorPacks = ".";
  const TString nameCondorPack = "preselect";

  const TString nameDate = "20210429";

  auto funOpenFile = [pathCondorPacks, nameCondorPack](const TString nameDatagroup, const TString nameClusterID, const TString nameTT)->TFile*{
    return TFile::Open(pathCondorPacks + "/" + nameCondorPack + "/" + "output_" + nameDatagroup + "_" + nameTT + "_" + nameClusterID + "_hist.root");
  };

  const std::vector<std::pair<TString, TString>> vNameDatagroupAndClusterIdSignal({
      std::pair<TString, TString>("signal_Mx2-150_Mv-500_Mx1-1_ctau-1", nameDate),
      std::pair<TString, TString>("signal_Mx2-1_Mv-500_Mx1-0p1_ctau-1", nameDate),
      std::pair<TString, TString>("signal_Mx2-50_Mv-500_Mx1-1_ctau-10", nameDate),
      });

  const Size_t nDatagroupsSignal = vNameDatagroupAndClusterIdSignal.size();

  const std::vector<std::pair<TString, TString>> vNameDatagroupAndClusterIdBackground({
      std::pair<TString, TString>("DYJets", "1400584"),
      std::pair<TString, TString>("TT", "1375491"),
      });

  const Size_t nDatagroupsBackground = vNameDatagroupAndClusterIdBackground.size();

  auto funRatioFromBoolHist = [](const TH1 *const hist)->Double_t {
    const Double_t valFalse =  static_cast<Double_t>(hist->GetBinContent(1));
    const Double_t valTrue = static_cast<Double_t>(hist->GetBinContent(2));
    if (valFalse == 0 && valTrue == 0) {
      return -1;
    }
    return valTrue / (valFalse + valTrue);
  };

  std::vector<std::array<std::array<TFile *, 2>, 2>> vaaFilesSignal(nDatagroupsSignal, std::array<std::array<TFile *, 2>, 2>());

  std::vector<std::array<std::array<TFile *, 2>, 2>> vaaFilesBackground(nDatagroupsBackground, std::array<std::array<TFile *, 2>, 2>());

  for (Size_t iDatagroup = 0; iDatagroup < nDatagroupsSignal; ++iDatagroup) {
    const TString& nameDatagroup = vNameDatagroupAndClusterIdSignal[iDatagroup].first;
    const TString& nameClusterID = vNameDatagroupAndClusterIdSignal[iDatagroup].second;
    for (Byte_t iLepton=0; iLepton<2; iLepton++) {
      for (Byte_t iAK=0; iAK<2; iAK++) {
        vaaFilesSignal[iDatagroup][iLepton][iAK] = funOpenFile(
            nameDatagroup,
            nameClusterID,
            (TString)"PreselectedMatching" + (iAK ? "FatJet" : "Jet") + namesLepton[iLepton]
            );
      }
    }
  }

  for (Size_t iDatagroup = 0; iDatagroup < nDatagroupsBackground; ++iDatagroup) {
    const TString& nameDatagroup = vNameDatagroupAndClusterIdBackground[iDatagroup].first;
    const TString& nameClusterID = vNameDatagroupAndClusterIdBackground[iDatagroup].second;
    for (Byte_t iLepton=0; iLepton<2; iLepton++) {
      for (Byte_t iAK=0; iAK<2; iAK++) {
        vaaFilesBackground[iDatagroup][iLepton][iAK] = funOpenFile(
            nameDatagroup,
            nameClusterID,
            (TString)"PreselectedMatching" + (iAK ? "FatJet" : "Jet") + namesLepton[iLepton]
            );
      }
    }
  }

  for (Byte_t iLepton = 0; iLepton < 2; iLepton++) {
    for (Byte_t iAK = 0; iAK < 2; iAK++) {
      const TString nameTT = TString("Preselected") + (iAK ? "FatJet" : "Jet") + namesLepton[iLepton];
      if (recreate) {
        for (Size_t iDatagroupSignal = 0; iDatagroupSignal < nDatagroupsSignal; ++iDatagroupSignal) {
          const TString & nameDatagroupSignal = vNameDatagroupAndClusterIdSignal[iDatagroupSignal].first;
          TFile *tfOut = TFile::Open(
                pathCondorPacks + "/" + nameCondorPack
                + "/" + "output_punzi" + "_" + nameDatagroupSignal + "_" + nameTT + "_" + nameDate + ".root",
                "RECREATE");
          tfOut->Close();
        }
      }
      for (TObject *&& keyRawBackgroundFirst: *vaaFilesBackground[0][iLepton][iAK]->GetListOfKeys()) {
        if (!keyRawBackgroundFirst || keyRawBackgroundFirst->IsZombie()) {
          continue;
        }
        TObject *&& histRawBackgroundFirst = static_cast<TKey*>(keyRawBackgroundFirst)->ReadObj();
        if (!histRawBackgroundFirst
            || histRawBackgroundFirst->IsZombie()
            || !histRawBackgroundFirst->IsA()->InheritsFrom("TH1")) {
          if (debug) std::cout << keyRawBackgroundFirst->GetName() << " is not a readable histogram in " << vNameDatagroupAndClusterIdBackground[0].first << std::endl;
          continue;
        }
        const TString nameHist = keyRawBackgroundFirst->GetName();
        const TString nameLeaf(nameHist(1, nameHist.Length()));
        std::vector<TH1 *> vHistBackgrounds(nDatagroupsBackground, nullptr);
        vHistBackgrounds[0] = static_cast<TH1 *>(histRawBackgroundFirst);
        Bool_t areKeysInAllBackgroundFiles = true;
        for (Ssiz_t iDatagroupBackground = 1; iDatagroupBackground < nDatagroupsBackground; ++iDatagroupBackground) {
          TObject *&& keyRaw = vaaFilesBackground[iDatagroupBackground][iLepton][iAK]->GetKey(nameHist);
          if (!keyRaw || keyRaw->IsZombie()) {
            areKeysInAllBackgroundFiles = false;
            if (debug) std::cout << keyRawBackgroundFirst->GetName() << " not available in " << vNameDatagroupAndClusterIdBackground[iDatagroupBackground].first << std::endl;
            break;
          }
          TObject *&& histRaw = static_cast<TKey*>(keyRaw)->ReadObj();
          if (!histRaw || histRaw->IsZombie()
              || !histRaw->IsA()->InheritsFrom("TH1")) {
          if (debug) std::cout << keyRawBackgroundFirst->GetName() << " is not a readable histogram in " << vNameDatagroupAndClusterIdBackground[iDatagroupBackground].first << std::endl;
            areKeysInAllBackgroundFiles = false;
            break;
          }
          vHistBackgrounds[iDatagroupBackground] = static_cast<TH1 *>(histRaw);
        }
        if (!areKeysInAllBackgroundFiles) {
          continue;
        }
        Double_t binDensity = 0;
        if (nameLeaf.BeginsWith("n")) {
          if (debug) std::cout << "Found n*" << std::endl;
          binDensity = 1;
        } else if (nameLeaf.BeginsWith("is") || nameLeaf.BeginsWith("are") || nameLeaf.BeginsWith("have") || nameLeaf.BeginsWith("has")) {
          continue;
        } else if (nameLeaf.Contains("MET_") || nameLeaf.BeginsWith("SoftActivityJetHT")) {
          if (debug) std::cout << "Found *MET_*" << std::endl;
          binDensity = TMath::Power(10, TMath::Nint(-TMath::Log10(vHistBackgrounds[0]->GetBinWidth(1))));
        } else {
          continue;
        }
        std::vector<TH1 *> vHistCumuBackground(nDatagroupsBackground, nullptr);
        for (Size_t iDatagroupBackground = 0; iDatagroupBackground < nDatagroupsBackground; ++iDatagroupBackground) {
          vHistCumuBackground[iDatagroupBackground] = vHistBackgrounds[iDatagroupBackground]->GetCumulative(false);
        }
        for (Size_t iDatagroupSignal = 0; iDatagroupSignal < nDatagroupsSignal; ++iDatagroupSignal) {
          const TString nameDatagroupSignal = vNameDatagroupAndClusterIdSignal[iDatagroupSignal].first;
          TObject *&& keyRawSignal = vaaFilesSignal[iDatagroupSignal][iLepton][iAK]->GetKey(nameHist);
          if (!keyRawSignal || keyRawSignal->IsZombie()) {
            if (debug) std::cout << keyRawBackgroundFirst->GetName() << " not available in " << vNameDatagroupAndClusterIdSignal[iDatagroupSignal].first << std::endl;
            continue;
          }
          TObject *&& histRawSignal = static_cast<TKey *>(keyRawSignal)->ReadObj();
          if (!histRawSignal || histRawSignal->IsZombie()
              || !histRawSignal->IsA()->InheritsFrom("TH1")) {
            if (debug) std::cout << keyRawBackgroundFirst->GetName() << " is not a readable histogram in " << vNameDatagroupAndClusterIdSignal[iDatagroupSignal].first << std::endl;
            continue;
          }
          if (debug) std::cout << "Computing " << nameHist << std::endl;
          TFile *tfOut = TFile::Open(
              pathCondorPacks + "/" + nameCondorPack
              + "/" + "output_punzi" + "_" + nameDatagroupSignal + "_" + nameTT + "_" + nameDate + ".root",
              "UPDATE");
          TH1 *histSignal = static_cast<TH1 *>(histRawSignal);
          Double_t lowerSignal = histSignal->GetBinLowEdge(1);
          Ssiz_t nBinsSignal = histSignal->GetNbinsX();
          std::vector<Long64_t> vI0BinBackgrounds(nDatagroupsBackground, 0);
          for (Size_t iDatagroupBackground = 0; iDatagroupBackground < nDatagroupsBackground; ++iDatagroupBackground) {
            const Double_t lowerBackground = vHistBackgrounds[iDatagroupBackground]->GetBinLowEdge(1);
            vI0BinBackgrounds[iDatagroupSignal] = static_cast<Long64_t>(TMath::Nint((lowerSignal - lowerBackground) * binDensity));
          }
          // Use hist1Cumu = hist1->GetCumulative(false) to get unnormalized reversely cumulative histogram about hist1
          // Use hist1Cumu->Scale(1./(hist1->GetBinContent(1))) to normalize
          auto histPunzi = histSignal->GetCumulative(false);
          histPunzi->SetName("hPunzi_" + nameLeaf);
          histPunzi->SetTitle("Punzi_" + nameLeaf);
          if (debug && histPunzi->GetNbinsX() != nBinsSignal) {
            Fatal("punzi", "Binnumber between histPunzi and histSignal (%s) does not match", nameHist.Data());
          }
          histPunzi->SetXTitle((TString)"Punzi significance of " + histSignal->GetXaxis()->GetTitle());
          Double_t sumSignal = histPunzi->GetBinContent(1);
          if (sumSignal <= __FLT_EPSILON__ * 1024) {
            sumSignal = 1;
          }
          for (Ssiz_t iBinSignal = 1; iBinSignal <= nBinsSignal; ++iBinSignal) {
            Double_t sumBackgrounds = 0.;
            for (Ssiz_t iDatagroupBackground = 0; iDatagroupBackground < nDatagroupsBackground; ++iDatagroupBackground) {
              const Ssiz_t nBinsBackground = vHistBackgrounds[iDatagroupBackground]->GetNbinsX();
              const Ssiz_t iBinBackground = vI0BinBackgrounds[iDatagroupBackground] + iBinSignal;
              sumBackgrounds += (iBinBackground >= 1 && iBinBackground <= nBinsBackground)
                                ? vHistCumuBackground[iDatagroupBackground]->GetBinContent(iBinBackground)
                                : 0;
            }
            histPunzi->SetBinContent(iBinSignal, histPunzi->GetBinContent(iBinSignal) / (sumSignal * (TMath::Sqrt(sumBackgrounds) + 1)));
          }
          histPunzi->Write();
          tfOut->Close();
        }
      }
    }
  }

}


int main(int argc, char* argv[]) {
  punzi_20210429();
}
