#define DEBUGGING true

#include <vector>
#include <iostream>
#include <fstream>
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
#include <cmath>

using namespace std;

template<typename TypeNumeric>
TypeNumeric myHeaviside(TypeNumeric x, TypeNumeric xShift=0, TypeNumeric yCentral=1) {
    return x == xShift ? yCentral : (x > xShift ? 1 : 0);
}

/*
TString fStrFormulaFormat(const char* strFormulaContent, const char* nameX="x", const char* nameP="p", const char* nameTypeNumeric="double") {
    TString formula = "[&](";
    formula += (TString)nameTypeNumeric + " *" + nameX + ", ";
    formula += (TString)nameTypeNumeric + " *" + nameP + ")";
    formula += "{return ";
    formula += strFormulaContent;
    formula += ";}";
    return formula;
}
*/
/*
TString fStrExpAndHillExpPart(bool inverseHorizotnal=false, bool inverseVerticle=false, const char* strParameter1st="[1]", const char* strParameter2nd="[2]") {
    TString formula = inverseVerticle ? "-" : "";
    TString strXDisp = (TString)"(x-" + strParameter1st + ")";
    formula += strParameter2nd;
    formula += "*(";
    formula += (TString)"exp(" + (!inverseHorizotnal ? "-" : "") + strXDisp + "/p[3])";
    formula += ")";
    return formula;
}
*/
/*
TString fStrExpAndHill(bool inverseHorizotnal=false, bool inverseVerticle=false) {
    // myHeaviside<Double_t>((x - [1])*signHorizontal)*[2]*exp(([1] - x)*signHorizontal/[4]) 
    // + myHeaviside<Double_t>(([1] - x)*signHorizontal)*[2]*(1 + (x - [1])*signHorizontal/[3])
    
    TString formula = "";
    TString strXDisp = "(*x-p[0])";
    formula += (inverseVerticle ? "-" : "");
    formula += "p[1]*(";
    formula += (TString)"exp(" + (!inverseHorizotnal ? "-" : "") + strXDisp + "/p[3])";
    formula += (TString)"*myHeaviside<double>(";
    formula += (inverseHorizotnal ? (TString)"-" + strXDisp + ", 0, 0" : strXDisp);
    formula += ")";
    formula += (TString)"+(1+" + (inverseHorizotnal ? "-" : "") + strXDisp + "/p[2])";
    formula += (TString)"*myHeaviside<double>(";
    formula += (!inverseHorizotnal ? (TString)"-" + strXDisp + ", 0, 0" : strXDisp);
    formula += ")";
    formula += ")";
    // std::cout << *formula << std::endl;
    return fStrFormulaFormat(formula);
}
*/
/*
Double_t fSumHist(TH1F* hist, Int_t *pBinStart=NULL, Int_t *pBinEnd=NULL) {
    Int_t* pBinStartGot = pBinStart;
    if (pBinStart == nullptr) {
        pBinStart = new Int_t(1);
    }
    Int_t* pBinEndGot = pBinEnd;
    if (pBinEnd == nullptr) { 
        pBinEndGot = new Int_t(hist->GetNbinsX()+1);
    }
    Stat_t sum=0;
    for (Int_t i=*pBinStartGot;i<*pBinEndGot;i++) sum+=hist->GetBin(i);
    return sum;
}
*/

Double_t fSumHist(TH1F* phist, Int_t iBinStart, Int_t iBinEnd) {
    Stat_t sum=0;
    for (Int_t i=iBinStart;i<iBinEnd;i++) sum+=phist->GetBinContent(i);
    if (DEBUGGING) {
        cout << "Summing " << phist->GetName() << ":" << endl;
        cout << "from bin " << iBinStart << " to " << iBinEnd << endl;
        cout << "sum: " << sum << endl << endl;
    }
    return sum;
}

Double_t fSumHist(TH1F* phist, Int_t iBinStart=1) {
    return fSumHist(phist, iBinStart, phist->GetNbinsX()+1);
}

Double_t fMaxHist(TH1F* phist, Int_t iBinStart, Int_t iBinEnd) {
    Stat_t maximum = 0;
    Double_t contentNow;
    for (Int_t i=iBinStart; i<iBinEnd; i++) {
        contentNow = phist->GetBinContent(i);
        if (maximum < contentNow) {
            maximum = contentNow;
        }
    }
    return maximum;
}

Double_t fMaxHist(TH1F* phist, Int_t iBinStart=1) {
    return fMaxHist(phist, iBinStart, phist->GetNbinsX()+1);
}

/*
class clDoubleExpAndHill {
public:
    bool inverseHorizotnal=false;
    bool inverseVerticle=false;
    clDoubleExpAndHill(bool inverseHorizotnal=false, bool inverseVerticle=false) {
        clDoubleExpAndHill::inverseHorizotnal = inverseHorizotnal;
        clDoubleExpAndHill::inverseVerticle = inverseVerticle;
    }
    Double_t fDoubleExpAndHill(Double_t *x, Double_t *p) {
        Double_t xMinusP1Adjusted = inverseHorizotnal ? p[0] - *x : *x - p[0];
        return myHeaviside<Double_t>(xMinusP1Adjusted)*p[1]*exp(-xMinusP1Adjusted/p[3]) + myHeaviside<Double_t>(xMinusP1Adjusted)*p[1]*(1 + xMinusP1Adjusted/p[2]);
    }
};
*/

inline Double_t fDoubleExpAndHillMeta(Double_t *x, Double_t *p, bool inverseHorizotnal=false) {
    Double_t xMinusP1Adjusted = inverseHorizotnal ? p[0] - *x : *x - p[0];
    // Double_t valMyHeavisideAdjusted = inverseHorizotnal ? myHeaviside<Double_t>(xMinusP1Adjusted, 0,  0) : myHeaviside<Double_t>(xMinusP1Adjusted);
    // return valMyHeavisideAdjusted*p[1]*exp(-xMinusP1Adjusted/p[3]) + valMyHeavisideAdjusted*p[1]*(1 + xMinusP1Adjusted/p[2]);
    return (xMinusP1Adjusted >= 0) ?  p[1]*exp(-xMinusP1Adjusted/p[3]) : p[1]*(1 + xMinusP1Adjusted/p[2]);
}

/*
auto fLambdaExpAndHill(bool inverseHorizotnal=false, bool inverseVerticle=false) {
    auto lambdaXMinusAdjusted = inverseHorizotnal ? [&](Double_t *x, Double_t *p){p[0]-*x;} : [&](Double_t *x, Double_t *p){return *x - p[0];};
    return [&](Double_t *x, Double_t *p){return p[1]*(exp(-lambdaXMinusAdjusted(x, p)/p[3])*myHeaviside<Double_t>(lambdaXMinusAdjusted(x, p))+(1+lambdaXMinusAdjusted(x, p)/p[2])*myHeaviside<Double_t>(-lambdaXMinusAdjusted(x, p), 0, 0));};
}
*/

void setParametersExpAndHill(TH1F* phist, Double_t parametersTF[4], bool inverseHorizotnal, Double_t xmin, Double_t xmax) { //Double_t *pxmin=NULL, Double_t *pxmax=NULL) {
    Int_t maximumBin = phist->GetMaximumBin();
    Double_t binWidth = phist->GetBinWidth(maximumBin);
    parametersTF[0] = binWidth * maximumBin;
    parametersTF[1] = phist->GetBinContent(maximumBin);
    parametersTF[2] = inverseHorizotnal ? 1 - parametersTF[0] : parametersTF[0];
    Double_t xlow = phist->GetBinLowEdge(1);
    Int_t numBins = phist->GetNbinsX();
    // Double_t xup = xlow + binWidth * binNumber;
    Int_t iMin, iMax;
    // iMin = (xmin - xlow >= binWidth/2) ? (xmin - xlow)//binWidth : 1;
    Double_t areaExpPart;
    Int_t deltaiDetectRange = numBins/100;
    iMin = 1+TMath::Max(0, (Int_t)((xmin - xlow)/binWidth));
    iMax = TMath::Min(numBins, (Int_t)(xmax/binWidth));
    if (not inverseHorizotnal) {
        if (iMax >= numBins-deltaiDetectRange) {
            areaExpPart = binWidth * fSumHist(phist, maximumBin);
        } else {
            areaExpPart = binWidth * fSumHist(phist, maximumBin, iMax) / (1 - fMaxHist(phist, iMax-deltaiDetectRange, iMax+deltaiDetectRange+1)/(deltaiDetectRange*2+1));
        }
    } else {
        if (iMin <= TMath::Min(deltaiDetectRange, maximumBin-deltaiDetectRange)) {
            areaExpPart = binWidth * fSumHist(phist, iMin, maximumBin+1);
        } else {
            areaExpPart = binWidth * fSumHist(phist, iMin, maximumBin+1) / (1 - fMaxHist(phist, iMin-deltaiDetectRange, iMin+deltaiDetectRange+1)/(deltaiDetectRange*2+1));
        }
    }
    parametersTF[3] = areaExpPart / parametersTF[1];
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
    if (DEBUGGING) {
        string strSepValues = "\n";
        string strEndValues = "\n";
        cout << "hist: " << phist->GetName() << endl;
        cout << "maximumBin: " << maximumBin << endl;
        cout << "maximum: " << phist->GetMaximum() << endl;
        cout << "binWidth: " << binWidth << endl;
        cout << "parametersTF:" << endl;
        cout << "xlow: " << xlow << endl;
        cout << "numBins: " << numBins << endl;
        cout << "xmin: " << xmin << " xmax: " << xmax << endl;
        cout << "iMin: " << iMin << " iMax: " << iMax << endl;
        cout << "deltaiDetectRange: " << deltaiDetectRange << endl;
        cout << "areaExpPart: " << areaExpPart << endl;
        for (int i=0; i<4-1; i++) {
            cout << parametersTF[i] << strSepValues;
        }
        cout << parametersTF[3] << strEndValues << endl;
        }
    
}

void setParametersExpAndHill(TH1F* phist, TF1* ptf, bool inverseHorizotnal, Double_t xmin, Double_t xmax) {
    Double_t parametersTF[4];
    setParametersExpAndHill(phist, parametersTF, inverseHorizotnal, xmin, xmax);
    // ptf->SetParameters(parametersTF[0], parametersTF[1], parametersTF[2], parametersTF[3]);
    ptf->SetParameters(parametersTF);
    ptf->SetParLimits(0, TMath::Max(xmin, parametersTF[0]-(xmax-xmin)/20), TMath::Min(xmax, parametersTF[0]+(xmax-xmin)/20));
    ptf->SetParLimits(1, parametersTF[2]*0.8, parametersTF[1]*1.5);
    ptf->SetParLimits(2, 0, parametersTF[2]*2);
    ptf->SetParLimits(3, 0, inverseHorizotnal ? TMath::Max(parametersTF[0], parametersTF[3]*0.5) : parametersTF[3]*2);
}
/*
void setParametersExpAndHill(TH1F* phist, Double_t parametersTF[], bool inverseHorizotnal=false, Double_t xmin=0.f) {
    return setParametersExpAndHill(phist, parametersTF, inverseHorizotnal, xmin, phist->GetBinLowEdge(1)+phist->GetBinWidth(1)*phist->GetNbinsX());
}
*/
void setParametersExpAndHill(TH1F* phist, TF1* ptf, bool inverseHorizotnal=false) {
    // return fitParametersExpAndHill(hist, tf, inverseHorizotnal, xmin, hist->GetBarOffset()+hist->GetBarWidth());
    if (DEBUGGING) {
        cout << "tfunction " << ptf->GetName() << ": " << endl;
        cout << "(xmin, xmax, npar): " << ptf->GetXmin() << ", " << ptf->GetXmax() << ", " << ptf->GetNpar() << endl;
        cout << endl;
    }
    return setParametersExpAndHill(phist, ptf, inverseHorizotnal, ptf->GetXmin(), ptf->GetXmax());
}

void fitHistInFile(const char* nameFileIn) {
    TFile* fileIn = new TFile(nameFileIn);
    if (DEBUGGING) {
        fileIn->ls();
    }
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
    auto doubleExpAndHill = [&](Double_t *x, Double_t *p){return fDoubleExpAndHillMeta(x, p);};
    auto doubleExpAndHillInversed = [&](Double_t *x, Double_t *p){return fDoubleExpAndHillMeta(x, p, true);};
    /*
    // Suitable when nameFileIn is "output_proper_Mx2-150_Mx1-1_20191224.root"
    TF1 *tfBetatauLab = new TF1("f_betatau_lab", doubleExpAndHill, 0, 2, 4);
    TF1 *tfChi2Beta = new TF1("f_chi2_beta", doubleExpAndHillInversed, 0, 1, 4);
    TF1 *tfChi2Gamma = new TF1("f_chi2_gamma", doubleExpAndHill, 1, 20, 4);
    TF1 *tfCtauLab = new TF1("f_ctau_lab", doubleExpAndHill, 0, 5, 4);
    TF1 *tfCtauProper = new TF1("f_ctau_proper", doubleExpAndHill, 0, 0.6, 4);
    */
    Double_t xLimitsBetaTauLab[2];
    Double_t xLimitsChi2Beta[2];
    Double_t xLimitsChi2Gamma[2];
    Double_t xLimitsCtauLab[2];
    Double_t xLimitsCtauProper[2];
    xLimitsBetaTauLab[0] = h_betatau_lab->GetBinLowEdge(1);
    xLimitsChi2Beta[0] = h_chi2beta->GetBinLowEdge(1);
    xLimitsChi2Gamma[0] = h_chi2gamma->GetBinLowEdge(1);
    xLimitsCtauLab[0] = h_ctau_lab->GetBinLowEdge(1);
    xLimitsCtauProper[0] = h_ctau_proper->GetBinLowEdge(1);
    xLimitsBetaTauLab[1] = xLimitsBetaTauLab[0] + h_betatau_lab->GetBinWidth(1)*h_betatau_lab->GetNbinsX();
    xLimitsChi2Beta[1] = xLimitsChi2Beta[0] + h_chi2beta->GetBinWidth(1)*h_chi2beta->GetNbinsX();
    xLimitsChi2Gamma[1] = xLimitsChi2Gamma[0] + h_chi2gamma->GetBinWidth(1)*h_chi2gamma->GetNbinsX();
    xLimitsCtauLab[1] = xLimitsCtauLab[0] + h_ctau_lab->GetBinWidth(1)*h_ctau_lab->GetNbinsX();
    xLimitsCtauProper[1] = xLimitsCtauProper[0] + h_ctau_proper->GetBinWidth(1)*h_ctau_proper->GetNbinsX();
    xLimitsChi2Beta[1] = 1.f;
    xLimitsChi2Gamma[0] = 1.f;
    TF1 *tfBetatauLab = new TF1("f_betatau_lab", doubleExpAndHill, xLimitsBetaTauLab[0], xLimitsBetaTauLab[1], 4);
    TF1 *tfChi2Beta = new TF1("f_chi2_beta", doubleExpAndHillInversed, xLimitsChi2Beta[0], xLimitsChi2Beta[1], 4);
    TF1 *tfChi2Gamma = new TF1("f_chi2_gamma", doubleExpAndHill, xLimitsChi2Gamma[0], xLimitsChi2Gamma[1], 4);
    TF1 *tfCtauLab = new TF1("f_ctau_lab", doubleExpAndHill, xLimitsCtauLab[0], xLimitsCtauLab[1], 4);
    TF1 *tfCtauProper = new TF1("f_ctau_proper", doubleExpAndHill, xLimitsCtauProper[0], xLimitsCtauProper[1], 4);
    setParametersExpAndHill(h_betatau_lab, tfBetatauLab);
    setParametersExpAndHill(h_chi2beta, tfChi2Beta, true);
    setParametersExpAndHill(h_chi2gamma, tfChi2Gamma);
    setParametersExpAndHill(h_ctau_lab, tfCtauLab);
    setParametersExpAndHill(h_ctau_proper, tfCtauProper);
    // h_betatau_lab->Draw();
    // tfBetatauLab->Draw("SAME");
    // h_chi2beta->Draw();
    // tfChi2Beta->Draw("SAME");
    // h_chi2gamma->Draw();
    // tfChi2Gamma->Draw("SAME");
    // h_ctau_lab->Draw();
    // tfCtauLab->Draw("SAME");
    // h_ctau_proper->Draw();
    // tfCtauProper->Draw("SAME");
    // TFitResultPtr tfitrpBetatauLab = h_betatau_lab->Fit(tfBetatauLab, "", "");
    // TFitResultPtr tfitrpChi2Beta = h_chi2beta->Fit(tfChi2Beta, "", "");
    // TFitResultPtr tfitrpChi2Gamma = h_chi2gamma->Fit(tfChi2Gamma, "", "");
    // TFitResultPtr tfitrpCtauLab = h_ctau_lab->Fit(tfCtauLab, "", "");
    TFitResultPtr tfitrpCtauProper = h_ctau_proper->Fit(tfCtauProper, "", "");
    TString pathImageFolder = "../out_images/output_proper_Mx2-150_Mx1-1_20200130";
    // TPaveText* fitlabel = new TPaveText(0.6,0.4,0.9,0.75,"NDC");
    // fileIn->Close();
}
