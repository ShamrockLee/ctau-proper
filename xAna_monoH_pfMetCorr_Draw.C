#define DEBUGGING true

#include <vector>
#include <iostream>
#include <string>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer_sh.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TStyle.h>

Double_t fSumHist(TH1F* phist, Int_t iBinStart, Int_t iBinEnd) {
    Stat_t sum=0;
    for (Int_t i=iBinStart;i<iBinEnd;i++) sum+=phist->GetBinContent(i);
    if (DEBUGGING) {
        using namespace std;
        cout << "Summing " << phist->GetName() << ":" << endl;
        cout << "from bin " << iBinStart << " to " << iBinEnd << endl;
        cout << "sum: " << sum << endl << endl;
    }
    return sum;
}

Double_t fSumHist(TH1F* phist, Int_t iBinStart=1) {
    return fSumHist(phist, iBinStart, phist->GetNbinsX()+1);
}

Double_t doubleDecayingExpWithTau(Double_t *x, Double_t *p) {
    return p[0] * TMath::Exp(-*x/p[1]);
}

void FitExp(TH1F* hist, TCanvas* ptc, const char* strFitOptions, const char* strFitDrawOptions, const char* nameTF, Double_t xmin, Double_t xmax) {
    TF1 *tfExpWithTau = new TF1(nameTF, &doubleDecayingExpWithTau, xmin, xmax, 2);
    tfExpWithTau->SetParameter(0, hist->GetMaximum());
    tfExpWithTau->SetParameter(1, hist->GetBinWidth(1) * fSumHist(hist) / tfExpWithTau->GetParameter(0));
    tfExpWithTau->SetParLimits(0, 0, tfExpWithTau->GetParameter(0)*1.5);
    tfExpWithTau->SetParLimits(1, 0, tfExpWithTau->GetParameter(1)*1.5);
    if (DEBUGGING) {
        std::cout << nameTF << " " << "initial parameters:" << std::endl;
        std::cout << "p1: " << tfExpWithTau->GetParameter(0) << std::endl;
        std::cout << "p2: " << tfExpWithTau->GetParameter(1) << std::endl;
    }
    ptc->Clear();
    hist->Fit(tfExpWithTau, strFitOptions, strFitDrawOptions);
}

void FitExp(TH1F* hist, TCanvas* ptc, const char* strFitOptions="L", const char* strFitDrawOptions="", const char* nameTF="tfExpWithTau") {
    return FitExp(hist, ptc, strFitOptions, strFitDrawOptions, nameTF, 0, hist->GetBinWidth(1)*hist->GetNbinsX());
}

void Main(const char* nameOutImageDir, TCanvas* c1 = nullptr, bool deleteCanvas=true) {
    bool reallyDeleteCanvas = false;
    if (c1 == nullptr) {
        c1 = new TCanvas;
        reallyDeleteCanvas = deleteCanvas;
    }
    TString nameFileHead = "pfMetCorr";
    TString nameFileTail = "20200224";
    TString namesFileVariable[] = {"Mx-1", "Mx2-1_Mx1-0p1_ctau-1", "Mx-150", "Mx2-150_Mx1-1_ctau-1"};
    TString namesHistogramSelected[] = {"hPfMetCorrPt", "hPfMetCorrPhi", "hPfMetCorrSig", 
        "hNGoodTHINJet", "hNGoodTHINBJet", 
        "hEleP4M", "hMuP4M", "hHPSTau_4MomentumM"};
    for (int i=0; i<4; i++) {
        TString nameFileVariable = namesFileVariable[i];
        TFile* ptfile = new TFile((TString)"output_" + nameFileHead + "_selected_" + nameFileVariable + "_" + nameFileTail + ".root");
        TString nameHistogram;
        for (int j=0;j<8;j++) {
            c1->Clear();
            nameHistogram = namesHistogramSelected[j];
            ((TH1F *)ptfile->Get(nameHistogram))->Draw();
            c1->Print((TString)nameOutImageDir + "/" + nameFileHead + "_selected_" + nameFileVariable + "_" +  nameHistogram + ".svg");
        }
        gStyle->SetOptFit(1111);
        c1->Clear();
        FitExp((TH1F *)ptfile->Get("hEleP4M"), c1);
        c1->Print((TString)nameOutImageDir + "/" + nameFileHead + "_selected_" + nameFileVariable + "_" +  "hEleP4M" + "_fit" + ".svg");
    }
}
