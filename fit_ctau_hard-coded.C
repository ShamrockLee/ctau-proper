#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <TString.h>
#include <TH1D.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TPavesText.h>

using namespace std;

template<typename TypeNumeric>
TypeNumeric myHeaviside(TypeNumeric x, TypeNumeric xShift=0, TypeNumeric yCentral=1) {
    return x == xShift ? yCentral : (x > xShift ? 1 : 0);
}

/*
char* strConcat(const char *s1, const char *s2) {
    // Adapted from David Heffernan's work at https://stackoverflow.com/a/8465083
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}
*/

TString fStrExpAndHillExpPart(bool inverseHorizotnal=false, bool inverseVerticle=false, const char* strParameter1st="[1]", const char* strParameter2nd="[2]") {
    TString formula = inverseVerticle ? "-" : "";
    TString strXDisp = (TString)"(x-" + strParameter1st + ")";
    formula += strParameter2nd;
    formula += "*(";
    formula += (TString)"exp(" + (!inverseHorizotnal ? "-" : "") + strXDisp + "/[4])";
    formula += ")";
    return formula;
}

TString fStrExpAndHill(bool inverseHorizotnal=false, bool inverseVerticle=false) {
    // myHeaviside<Double_t>((x - [1])*signHorizontal)*[2]*exp(([1] - x)*signHorizontal/[4]) 
    // + myHeaviside<Double_t>(([1] - x)*signHorizontal)*[2]*(1 + (x - [1])*signHorizontal/[3])
    // std::string* formula = new string();
    // std::string* strXMinusArg1 = (string*) "(x - [1])";
    // std::string* strXMinusArg1TimesSignHorizontal = signHorizontal > 0 ? strcat("-", strXMinusArg1) : strXMinusArg1;
    // std::string* strXMinusArg1TimesSignHorizontalNeg = signHorizontal < 0 ? strcat("-", strXMinusArg1) : strXMinusArg1;
    TString formula = "";
    TString strXDisp = "(x-[1])";
    formula += (inverseVerticle ? "-" : "");
    formula += "[2]*(";
    formula += (TString)"exp(" + (!inverseHorizotnal ? "-" : "") + strXDisp + "/[4])";
    formula += (TString)"*myHeaviside<Double_t>(";
    formula += (inverseHorizotnal ? (TString)"-" + strXDisp + ", 0, 0" : strXDisp);
    formula += ")";
    formula += (TString)"+(1+" + (inverseHorizotnal ? "-" : "") + strXDisp + "/[3])";
    formula += (TString)"*myHeaviside<Double_t>(";
    formula += (!inverseHorizotnal ? (TString)"-" + strXDisp + ", 0, 0" : strXDisp);
    formula += ")";
    formula += ")";
    /*
    if (signHorizontal > 0) {
        *formula = "myHeaviside<Double_t>((x - [1]))*[2]*exp(([1] - x)/[4]) "\
        "+ myHeaviside<Double_t>(([1] - x))*[2]*(1 + (x - [1])/[3])";
    } else {
        *formula = "myHeaviside<Double_t>(-(x - [1]))*[2]*exp(-([1] - x)/[4]) "\
        "+ myHeaviside<Double_t>(-([1] - x))*[2]*(1 - (x - [1])/[3])";
    }
    */
    // std::cout << *formula << std::endl;
    return formula;

}

/*
Double_t fSumHist(TH1F* hist, Int_t *pBinStart=NULL, Int_t *pBinEnd=NULL) {
    Int_t* pBinStartGot = (pBinStart == nullptr) ? (Int_t*)0 : pBinStart;
    Int_t* pBinEndGot = pBinEnd;
    if (pBinEnd == nullptr) { 
        *pBinEndGot = hist->GetNbinsX();
    }
    Stat_t sum=0;
    for (Int_t i=*pBinStartGot;i<*pBinEndGot;i++) sum+=hist->GetBin(i);
    return sum;
}
*/

Double_t fSumHist(TH1F* phist, Int_t iBinStart, Int_t iBinEnd) {
    Stat_t sum=0;
    for (Int_t i=iBinStart;i<iBinEnd;i++) sum+=phist->GetBin(i);
    return sum;
}

Double_t fSumHist(TH1F* phist, Int_t iBinStart=0) {
    return fSumHist(phist, iBinStart, phist->GetNbinsX());
}

void setParametersExpAndHill(TH1F* phist, Double_t parametersTF[], bool inverseHorizotnal, Double_t xmin, Double_t xmax) { //Double_t *pxmin=NULL, Double_t *pxmax=NULL) {
    Int_t maximumBin = phist->GetMaximumBin();
    Double_t binWidth = phist->GetBinWidth(maximumBin);
    parametersTF[0] = binWidth * maximumBin;
    parametersTF[1] = phist->GetBin(maximumBin);
    parametersTF[2] = inverseHorizotnal ? 1 - parametersTF[0] : parametersTF[0];
    /*
    Double_t xmin;
    Double_t xmax;
    if (pxmin == nullptr) {
        xmin = hist->GetBarOffset();
    } else {
        xmin = *pxmin;
    }
    if (pxmax == nullptr) {
        xmax = xmin + hist->GetBarWidth();
    } else {
        xmax = *pxmax;
    }
    */
    // TF1 *tfExpPart = new TF1("tfExpPart", fStrExpAndHillExpPart(inverseHorizotnal), xmin, xmax);
    parametersTF[3] = binWidth * fSumHist(phist, xmin, xmax) / parametersTF[1];
}

void setParametersExpAndHill(TH1F* phist, TF1* ptf, bool inverseHorizotnal, Double_t xmin, Double_t xmax) {
    Double_t parametersTF[4];
    setParametersExpAndHill(phist, parametersTF, inverseHorizotnal, xmin, xmax);
    ptf->SetParameters(parametersTF[0], parametersTF[1], parametersTF[2], parametersTF[3]);
}

void setParametersExpAndHill(TH1F* phist, Double_t parametersTF[], bool inverseHorizotnal=false, Double_t xmin=0.f) {
    return setParametersExpAndHill(phist, parametersTF, inverseHorizotnal, xmin, phist->GetBarOffset()+phist->GetBarWidth());
}

void setParametersExpAndHill(TH1F* phist, TF1* ptf, bool inverseHorizotnal=false) {
    // return fitParametersExpAndHill(hist, tf, inverseHorizotnal, xmin, hist->GetBarOffset()+hist->GetBarWidth());
    return setParametersExpAndHill(phist, ptf, inverseHorizotnal, ptf->GetXmin(), ptf->GetXmax());
}

void fitHistInFile(const char* nameFileIn) {
    TFile* fileIn = new TFile(nameFileIn);
    TH1F* h_betatau_lab;
    TH1F* h_chi2beta;
    TH1F* h_chi2gamma;
    TH1F* h_ctau_lab;
    TH1F* h_ctau_proper;
    fileIn->GetObject("h_betatau_lab", h_betatau_lab);
    fileIn->GetObject("h_chi2_beta", h_chi2beta);
    fileIn->GetObject("h_chi2_gamma", h_chi2gamma);
    fileIn->GetObject("h_ctau_lab", h_ctau_lab);
    fileIn->GetObject("h_ctau_proper", h_ctau_proper);
    // Suitable when nameFileIn is "output_proper_Mx2-150_Mx1-1_20191224.root"
    TF1 *tfBetatauLab = new TF1("f_betatau_lab", fStrExpAndHill(), 0, 2);
    TF1 *tfChi2Beta = new TF1("f_chi2_beta", fStrExpAndHill(true), 0, 1);
    TF1 *tfChi2Gamma = new TF1("f_chi2_gamma", fStrExpAndHill(), 1, 20);
    TF1 *tfCtauLab = new TF1("f_ctau_lab", fStrExpAndHill(), 0, 5);
    TF1 *tfCtauProper = new TF1("f_ctau_proper", fStrExpAndHill(), 0, 0.6);
    setParametersExpAndHill(h_betatau_lab, tfBetatauLab);
    setParametersExpAndHill(h_chi2beta, tfChi2Beta, true);
    setParametersExpAndHill(h_chi2gamma, tfChi2Gamma);
    setParametersExpAndHill(h_ctau_lab, tfCtauLab);
    setParametersExpAndHill(h_ctau_proper, tfCtauProper);
    h_betatau_lab->Fit(tfBetatauLab, "L");
    h_chi2beta->Fit(tfChi2Beta, "L");
    h_chi2gamma->Fit(tfChi2Gamma, "L");
    h_ctau_lab->Fit(tfCtauLab, "L");
    h_ctau_proper->Fit(tfCtauProper, "L");
    // TPaveText* fitlabel = new TPaveText(0.6,0.4,0.9,0.75,"NDC");
}
