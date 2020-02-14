#define DEBUGGING true

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <TString.h>
#include <TH1D.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TPavesText.h>
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

void fitProperHistInFile(const char* nameFileIn, const char* nameOutImageDir, TCanvas* c1 = nullptr, bool deleteCanvas=true) {
    TFile* fileIn = new TFile(nameFileIn);
    if (DEBUGGING) {
        fileIn->ls();
    }
    TH1F* hCtauProper;
    fileIn->GetObject("h_ctau_proper", hCtauProper);
    Double_t xlowCtauProper = hCtauProper->GetBinLowEdge(1);
    Double_t binWidthCtauProper = hCtauProper->GetBinWidth(1);
    Int_t maximumCtauProper = hCtauProper->GetMaximum();
    Double_t xLimitsCtauProper[2];
    xLimitsCtauProper[0] =  xlowCtauProper;
    xLimitsCtauProper[1] = xLimitsCtauProper[0] + hCtauProper->GetBinWidth(1) * hCtauProper->GetNbinsX();
    const int nPars = 2;
    TF1 *tfCtauProper = new TF1("f_ctau_proper", [](Double_t *x, Double_t *p){return p[0]*TMath::Exp(-*x/p[1]);}, xLimitsCtauProper[0], xLimitsCtauProper[1], nPars);
    Double_t parametersPredefined[nPars];
    parametersPredefined[0] = maximumCtauProper;
    parametersPredefined[1] = binWidthCtauProper * fSumHist(hCtauProper) / maximumCtauProper;
    tfCtauProper->SetParameters(parametersPredefined);
    tfCtauProper->SetParLimits(0, tfCtauProper->GetParameter(0)*0.5, tfCtauProper->GetParameter(0)*1.5);
    if (DEBUGGING) {
        std::cout << "Initial parameters of " << tfCtauProper->GetName() << ":" << std::endl;
        std::cout << tfCtauProper->GetParameter(0) << ", " << tfCtauProper->GetParameter(1) << std::endl;
        std::cout << std::endl;
    }
    bool reallyDeleteCanvas = false;
    if (c1 == nullptr) {
        c1 = new TCanvas;
        reallyDeleteCanvas = deleteCanvas;
    }
    gStyle->SetOptFit(1111);
    const int nStrFitOptions = 2;
    std::string strFitOptions[nStrFitOptions] = {"", "L"};
    for (int i=0; i<nStrFitOptions; i++) {
        c1->Clear();
        tfCtauProper->SetParameters(parametersPredefined);
        hCtauProper->Fit(tfCtauProper, (TString) strFitOptions[i]);
        c1->Print((TString)nameOutImageDir + "/fit_ctau_proper_" + ((strFitOptions[i].empty()) ? "none" : strFitOptions[i]) + ".svg");
    }
    if (reallyDeleteCanvas) {
        c1->Close();
        delete c1;
        delete tfCtauProper;
        fileIn->Close();
        delete fileIn;
        delete hCtauProper;
    }
}
