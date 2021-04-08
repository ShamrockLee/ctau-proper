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

#include <vector>
#include <array>
#include <iostream>

void punzi_20210408(const Bool_t debug=false) {
  const TString namesLepton[] = {"Electron", "Muon", "Tau"};  //< the name of the leptons (Xxx)
  const TString namesLeptonLower[] = {"electron", "muon", "tau"};  //< the name of the leptons (xxx)
  const TString namesLeptonNota[] = {"e", "mu", "tau"};

  const TString pathCondorPacks = "/home/shamrock/Documents/Project_HEP/ROOT_assignment_20191210/ctau-proper/lxplus_HTcondor";
  const TString nameCondorPack = "preselect";

  const TString pathOutImages = "/home/shamrock/Documents/Project_HEP/ROOT_assignment_20191210/out_images";
  const TString nameFolderToday = "output_" + nameCondorPack + "_20210408";

  auto funOpenFile = [pathCondorPacks, nameCondorPack](const TString nameDatagroup, const TString nameClusterID, const TString nameTT)->TFile*{
    return TFile::Open(pathCondorPacks + "/" + nameCondorPack + "/" + "output_" + nameDatagroup + "_" + nameTT + "_" + nameClusterID + "_hist.root");
  };

  const TString nameDatagroupSignalHeavy = "signal_Mx2-150_Mv-500_Mx1-1_ctau-1";
  const TString nameClusterIDSignalHeavy = "20210401";

  const TString nameDatagroupDYJets = "DYJets";
  const TString nameClusterIDDYJets = "1325718";

  const TString nameDatagroupTT = "TT";
  const TString nameClusterIDTT = "1345864";

  const TString suffixImage = "svg";

  auto funRatioFromBoolHist = [](const TH1 *const hist)->Double_t {
    const Double_t valFalse =  static_cast<Double_t>(hist->GetBinContent(1));
    const Double_t valTrue = static_cast<Double_t>(hist->GetBinContent(2));
    if (valFalse == 0 && valTrue == 0) {
      return -1;
    }
    return valTrue / (valFalse + valTrue);
  };

  std::array<std::array<TFile *, 2>, 2> aaFileSignalHeavyPreselected;
  {
    const TString nameDatagroup = nameDatagroupSignalHeavy;
    const TString nameClusterID = nameClusterIDSignalHeavy;
    for (Byte_t iLepton=0; iLepton<2; iLepton++) {
      for (Byte_t iAK=0; iAK<2; iAK++) {
        aaFileSignalHeavyPreselected[iLepton][iAK] = funOpenFile(nameDatagroup, nameClusterID, (TString)"PreselectedMatching" + (iAK ? "FatJet" : "Jet") + namesLepton[iLepton]);
      }
    }
  }
  std::array<std::array<TFile *, 2>, 2> aaFileDYJetsPreselected;
  {
    const TString nameDatagroup = nameDatagroupDYJets;
    const TString nameClusterID = nameClusterIDDYJets;
    for (Byte_t iLepton=0; iLepton<2; iLepton++) {
      for (Byte_t iAK=0; iAK<2; iAK++) {
        aaFileDYJetsPreselected[iLepton][iAK] = funOpenFile(nameDatagroup, nameClusterID, (TString)"PreselectedMatching" + (iAK ? "FatJet" : "Jet") + namesLepton[iLepton]);
      }
    }
  }
  std::array<std::array<TFile *, 2>, 2> aaFileTTPreselected;
  {
    const TString nameDatagroup = nameDatagroupTT;
    const TString nameClusterID = nameClusterIDTT;
    for (Byte_t iLepton=0; iLepton<2; iLepton++) {
      for (Byte_t iAK=0; iAK<2; iAK++) {
        aaFileTTPreselected[iLepton][iAK] = funOpenFile(nameDatagroup, nameClusterID, (TString)"PreselectedMatching" + (iAK ? "FatJet" : "Jet") + namesLepton[iLepton]);
      }
    }
  }

  gSystem->mkdir(pathOutImages + "/" + nameFolderToday);

  for (Byte_t iLepton = 0; iLepton < 2; iLepton++) {
    for (Byte_t iAK = 0; iAK < 2; iAK++) {
      for (TObject *keyRawSignal: *aaFileSignalHeavyPreselected[iLepton][iAK]->GetListOfKeys()) {
        TObject* histRawSignal = static_cast<TKey*>(keyRawSignal)->ReadObj();
        if (!histRawSignal || histRawSignal->IsZombie()) {
          continue;
        }
        if (!histRawSignal->IsA()->InheritsFrom("TH1")) {
          continue;
        }
        TString nameHist = histRawSignal->GetName();
        TString nameLeaf = nameHist(1, nameHist.Length() - 1);
        TH1 *histSignal = static_cast<TH1*>(histRawSignal);
        TH1 *histDYJets = static_cast<TH1*>(aaFileDYJetsPreselected[iLepton][iAK]->Get(nameHist));
        TH1 *histTT = static_cast<TH1*>(aaFileTTPreselected[iLepton][iAK]->Get(nameHist));
        if (histSignal->IsZombie() ||
            // !histDYJets || histDYJets->IsZombie() ||
            !histTT || histTT->IsZombie()) {
          continue;
        }
        const Double_t sumSignal = histSignal->Integral();
        if (sumSignal <= 0) {
          continue;
        }
        const Double_t lowerCorrectSignal = histSignal->GetBinCenter(1) - 0.5;
        const Int_t nBinsSignal = histSignal->GetNbinsX();
        Double_t binDensity = 0;
        if (nameLeaf.BeginsWith("n")) {
          binDensity = 1;
        } else if (nameLeaf.BeginsWith("is") || nameLeaf.BeginsWith("are") || nameLeaf.BeginsWith("have") || nameLeaf.BeginsWith("has")) {
          continue;
        } else if (nameLeaf.Contains("MET_") || nameLeaf.BeginsWith("SoftActivityJetHT")) {
          binDensity = TMath::Power(10, TMath::Nint(-TMath::Log10(histSignal->GetBinWidth(1))));
        } else {
          continue;
        }
        if (debug) std::cout << "nameLeaf: " << nameLeaf << std::endl;
        if (debug) std::cout << "nBinsSignal: " << nBinsSignal << std::endl;
        // TH2F *hhPunzi = new TH2F("Punzi_" + nameLeaf, (TString)"Punzi significance of " + histSignal->GetTitle(),
        // nBinsSignal, lowerCorrectSignal, lowerCorrectSignal + static_cast<Double_t>(nBinsSignal) / binDensity,
        // nBinsSignal, lowerCorrectSignal, lowerCorrectSignal + static_cast<Double_t>(nBinsSignal) / binDensity);
        const TString nameTT = TString("Preselected") + (iAK ? "FatJet" : "Jet") + namesLepton[iLepton];
        const TString nameHistPunzi = nameTT + "_" + nameLeaf;
        auto hPunzi = new TH1F(nameHistPunzi,  (TString)"Punzi significance of " + histSignal->GetTitle(),
            nBinsSignal, lowerCorrectSignal, lowerCorrectSignal + static_cast<Double_t>(nBinsSignal) / binDensity);
        hPunzi->SetXTitle( (TString)"Punzi significance of " + histSignal->GetXaxis()->GetTitle());
        Int_t jTT = TMath::Min(histTT->GetNbinsX(), TMath::Nint((histSignal->GetBinCenter(nBinsSignal) - histTT->GetBinCenter(1)) * binDensity) + 1);
        for (Int_t i = 1; i <= nBinsSignal; i++) {
          // for (Int_t j = i; j <= nBinsSignal; j++) {
          //   Double_t efficiencySignal = histSignal->Integral(i, j) / sumSignal;
          //   Double_t nBackground = 0.;
          //   {
          //     Int_t iTT = TMath::Nint((histSignal->GetBinCenter(i) - histTT->GetBinCenter(1)) * binDensity) + 1;
          //     Int_t jTT = TMath::Nint((histSignal->GetBinCenter(j) - histTT->GetBinCenter(1)) * binDensity) + 1;
          //     if (iTT >= 1 && jTT >= iTT) {
          //       nBackground += histTT->Integral(iTT, jTT);
          //     }
          //   }
          //   hhPunzi->SetBinContent(i, j, efficiencySignal / (1. + TMath::Sqrt(nBackground)));
          // }
          Double_t efficiencySignal = histSignal->Integral(i, nBinsSignal) / sumSignal;
          Double_t nBackground = 0.;
          Int_t iTT = TMath::Nint((histSignal->GetBinCenter(1) - histTT->GetBinCenter(1)) * binDensity) + i;
          // if (debug) std::cout << iTT << "\t" << jTT << std::endl;
          if (iTT >= 1 && iTT <= jTT)
            nBackground += histTT->Integral(iTT, jTT);
          hPunzi->SetBinContent(i, efficiencySignal / (1. + TMath::Sqrt(nBackground)));
        }
        if (!gPad) gROOT->MakeDefCanvas();
        gPad->SetLogy(false);
        // hhPunzi->Draw();
        hPunzi->Draw();
        gPad->SetTitle(nameHistPunzi);
        gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + nameHistPunzi + "." + suffixImage);
        // gPad->SetLogy(true);
        // gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "Punzi_" + nameLeaf + "_logy." + suffixImage);
      }
    }
  }

  // for (Byte_t iAK = 0; iAK<2; iAK++) {
  //   {
  //     if (debug) Info("plotting", "Getting names from file %s", arrFileSignalHeavyPreselectedElectron[iAK]->GetName());
  //     std::vector<TString> vName = {};
  //     std::vector<UInt_t> vIdx = {};
  //     refgetIdxHistOpenableRecent(arrFileSignalHeavyPreselectedElectron[iAK], vName, vIdx);
  //     for (const auto& nameHist: vName) {
  //       const TString nameLeaf = nameHist(1, nameHist.Length());
  //       if (nameLeaf.BeginsWith(iAK ? "FatJet_" : "Jet_")) {
  //         if (debug) Info("plotting", "Plotting %s", nameLeaf.Data());
  //         const TString tstrUnderscoreAndSuffix = nameLeaf(iAK ? 6 : 3, nameLeaf.Length());
  //         std::vector<Bool_t> vIsHist = {};
  //         std::vector<TH1*> vHist = {};
  //         std::vector<TString> vNameLeg = {};
  //         TString tstrXTitle = "";
  //         {
  //           const TString nameDatagroupShorter = "SignalHeavy";
  //           {
  //             const TString nameLeafToGet = TString::Format("%sMatched%s", iAK ? "FatJet" : "Jet", tstrUnderscoreAndSuffix.Data());
  //             vHist.push_back(static_cast<TH1*>(arrFileSignalHeavyAllMatchedElectron[iAK]->Get("h" + nameLeafToGet)));
  //             if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
  //             vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
  //             vNameLeg.push_back(Form("AllMatched%s_%s", iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
  //           }
  //           // for (Int_t nFirst=1; nFirst <= (iAK ? 6 : 8); nFirst++) {
  //           {
  //             Int_t nFirst = (iAK ? 2 : 6);
  //             const TString nameLeafToGet = TString::Format("%sFirst%d%s", iAK ? "FatJet" : "Jet", nFirst, tstrUnderscoreAndSuffix.Data());
  //             vHist.push_back(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[iAK]->Get("h" + nameLeafToGet)));
  //             if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
  //             vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
  //             vNameLeg.push_back(Form("Preselected%s_%s", iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
  //           }

  //           {
  //             const TString& nameLeafToGet = nameLeaf;
  //             vHist.push_back(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[iAK]->Get("h" + nameLeafToGet)));
  //             if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
  //             vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
  //             vNameLeg.push_back(Form("Preselected%s_%s", iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
  //             if (vIsHist.back()) {
  //               tstrXTitle = vHist.back()->GetXaxis()->GetTitle();
  //             }
  //           }
  //         }
  //         if (gPad) gPad->SetLogy(false);
  //         drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeaf, nameLeaf, tstrXTitle);
  //         gPad->SetTitle(nameLeaf);
  //         gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedAndAllMatched_" + nameLeaf + "." + suffixImage);
  //         gPad->SetLogy(true);
  //         gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedAndAllMatched_" + nameLeaf + "_logy." + suffixImage);
  //         {
  //           const TString nameDatagroupShorter = nameDatagroupDYJets;
            
  //           // for (Int_t nFirst=1; nFirst <= (iAK ? 6 : 8); nFirst++) {
  //           {
  //             Int_t nFirst = (iAK ? 2 : 6);
  //             const TString nameLeafToGet = TString::Format("%sFirst%d%s", iAK ? "FatJet" : "Jet", nFirst, tstrUnderscoreAndSuffix.Data());
  //             vHist.push_back(static_cast<TH1*>(arrFileDYJetsPreselectedEletron[iAK]->Get("h" + nameLeafToGet)));
  //             if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
  //             vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
  //             vNameLeg.push_back(Form("%s_Preselected%s_%s", nameDatagroupShorter.Data(), iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
  //           }

  //           {
  //             const TString nameLeafToGet = nameLeaf;
  //             vHist.push_back(static_cast<TH1*>(arrFileDYJetsPreselectedEletron[iAK]->Get("h" + nameLeafToGet)));
  //             if (debug && vHist.back() != nullptr) Info("plotting", "Plotting hist (%p) name %s", vHist.back(), vHist.back()->GetName());
  //             vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
  //             vNameLeg.push_back(Form("%s_Preselected%s_%s", nameDatagroupShorter.Data(), iAK ? "FatJet" : "Jet", nameLeafToGet.Data()));
  //           }
  //         }
  //         if (gPad) gPad->SetLogy(false);
  //         drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeaf, nameLeaf, tstrXTitle);
  //         gPad->SetTitle(nameLeaf);
  //         gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedAndAllMatched_" + nameLeaf + "vsDYJets" + "." + suffixImage);
  //         gPad->SetLogy(true);
  //         gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedAndAllMatched_" + nameLeaf + "vsDYJets" + "_logy." + suffixImage);
  //       }
  //       if (nameLeaf.EqualTo("SoftActivityJetHT")) {
  //         for (auto cstrNumber: {"", "2", "5", "10"}) {
  //           TString nameLeafToGet = nameLeaf + cstrNumber;
  //           std::vector<Bool_t> vIsHist = {};
  //           std::vector<TH1 *> vHist = {};
  //           std::vector<TString> vNameLeg = {};
  //           vHist.push_back(static_cast<TH1 *>(arrFileSignalHeavyPreselectedElectron[iAK]->Get("h" + nameLeafToGet)));
  //           vIsHist.push_back(vHist.back() != nullptr &&
  //                             !vHist.back()->IsZombie());
  //           vNameLeg.push_back(TString::Format("SignalHeavy_Electron_%s",
  //                                              iAK ? "FatJet" : "Jet"));
  //           vHist.push_back(static_cast<TH1 *>(
  //               arrFileDYJetsPreselectedEletron[iAK]->Get("h" + nameLeafToGet)));
  //           vIsHist.push_back(vHist.back() != nullptr &&
  //                             !vHist.back()->IsZombie());
  //           vNameLeg.push_back(TString::Format("DYJets_Electron_%s",
  //                                              iAK ? "FatJet" : "Jet"));
  //           if (gPad) gPad->SetLogy(false);
  //           drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeafToGet,
  //                             nameLeafToGet + " (Preselected)",
  //                             vHist.back()->GetXaxis()->GetTitle());
  //           gPad->SetTitle(nameLeaf);
  //           gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
  //                       "PreselectedJetElectron_" + nameLeaf + "vsDYJets" +
  //                       "." + suffixImage);
  //           gPad->SetLogy(false);
  //           gPad->Print(pathOutImages + "/" + nameFolderToday + "/" +
  //                       "PreselectedJetElectron_" + nameLeaf + "vsDYJets" +
  //                       "_logy." + suffixImage);
  //         }
  //       }
  //       if (nameLeaf.Contains("MET_")) {
  //         std::vector<Bool_t> vIsHist = {};
  //         std::vector<TH1*> vHist = {};
  //         std::vector<TString> vNameLeg = {};
  //         vHist.push_back(static_cast<TH1*>(arrFileSignalHeavyPreselectedElectron[iAK]->Get(nameHist)));
  //         vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
  //         vNameLeg.push_back(TString::Format("SignalHeavy_Electron_%s", iAK ? "FatJet" : "Jet"));
  //         vHist.push_back(static_cast<TH1*>(arrFileDYJetsPreselectedEletron[iAK]->Get(nameHist)));
  //         vIsHist.push_back(vHist.back() != nullptr && !vHist.back()->IsZombie());
  //         vNameLeg.push_back(TString::Format("DYJets_Electron_%s", iAK ? "FatJet" : "Jet"));
  //         if (gPad) gPad->SetLogy(false);
  //         drawHistsMultiple(vIsHist, vHist, vNameLeg, nameLeaf, nameLeaf + " (Preselected)", TString(vHist[0]->GetXaxis()->GetTitle()) + " (Pt > #)");
  //         gPad->SetTitle(nameLeaf);
  //         gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedJetElectron_" + nameLeaf + "vsDYJets" + "." + suffixImage);
  //         gPad->SetLogy(false);
  //         gPad->Print(pathOutImages + "/" + nameFolderToday + "/" + "PreselectedJetElectron_" + nameLeaf + "vsDYJets" + "_logy." + suffixImage);
  //       }
  //     }
  //   }
  // }
}


int main(int argc, char* argv[]) {
  punzi_20210408();
}
