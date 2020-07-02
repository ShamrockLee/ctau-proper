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
#include <TLegend.h>

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
Double_t doubleGrowingExpWithTau(Double_t *x, Double_t *p) {
    return p[0] * TMath::Exp(*x/p[1]);
}

void FitExp(TH1F* hist, TCanvas* ptc, bool isDecaying, const char* strFitOptions, const char* strFitDrawOptions, const char* nameTF, Double_t xmin, Double_t xmax) {
    Double_t (*doubleFunctionPointer)(Double_t*, Double_t*);
    doubleFunctionPointer = isDecaying ? doubleDecayingExpWithTau : doubleGrowingExpWithTau;
    TF1 *tfExpWithTau = new TF1(nameTF, doubleFunctionPointer, xmin, xmax, 2);
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

void FitExp(TH1F* hist, TCanvas* ptc, bool isDecaying=true, const char* strFitOptions="L", const char* strFitDrawOptions="", const char* nameTF="tfExpWithTau") {
    return FitExp(hist, ptc, isDecaying, strFitOptions, strFitDrawOptions, nameTF, 0, hist->GetBinWidth(1)*hist->GetNbinsX());
}

void chi2NbinsCompare(const TH1 *h1, const TH1 *h2, double &chi2, int &nbins, int binLo = -1, int binHi = -1) {

    printf("Chi2 Calculation: \n");

    chi2 = 0;
    nbins = 0;

    int nBins = 0;
    double evtThresh = 1e-6;

    double h1Err, h1Value, h2Err, h2Value, lowEdge, highEdge, binChi2;

    double binWidth = h1->GetXaxis()->GetBinWidth(1);

    if (binLo < 0 || binHi < 0) {
        binLo = 1;
        binHi = h1->GetNbinsX();
    }

    for (int i = binLo; i <= binHi; i++) {

        h1Value = h1->GetBinContent(i);
        h1Err = h1->GetBinError(i);

        h2Value = h2->GetBinContent(i);
        h2Err = h2->GetBinError(i);

        if (h1Value < evtThresh && h2Value < evtThresh)
            continue;
        if (h1Err < 1e-6 && h2Err < 1e-6)
            continue;
        lowEdge = h1->GetXaxis()->GetBinLowEdge(i);
        highEdge = h1->GetXaxis()->GetBinUpEdge(i);

        binChi2 = (h1Value - h2Value) / sqrt(h1Err * h1Err + h2Err * h2Err);
        binChi2 *= binChi2;

        chi2 += binChi2;

        // printf( "%d) [%d, %d] [%f, %f] h1: %f h2: %f h1Err: %f h2Err: %f chi2: %f total: %f \n ",
        // 	    nBins, i, i, lowEdge, highEdge,h1Value, h2Value, h1Err, h2Err, binChi2, chi2);

        nBins++;
    }
    int NDF = nBins;

    printf("Fit chi2/NDF = %f/%d, prob: %f\n", chi2, NDF, TMath::Prob(chi2, NDF) * 100);
    nbins = NDF;
}

void Main(std::string nameOutImageDir, TCanvas* c1 = nullptr, bool deleteCanvas=true) {
    bool reallyDeleteCanvas = false;
    if (c1 == nullptr) {
        c1 = new TCanvas;
        reallyDeleteCanvas = deleteCanvas;
    }
    gStyle->SetOptStat(111111);
    TString nameFileHead = "jets";
    TString nameFileTail = "20200702";
    const Int_t nFVarComparing = 2;
    const Int_t nFVarableGroups = 2;
    TString arrsNameFileVariable[nFVarableGroups][nFVarComparing];
    arrsNameFileVariable[0][0] = "Mx-1";
    arrsNameFileVariable[0][1] = "Mx2-1_Mx1-0p1_ctau-1";
    arrsNameFileVariable[1][0] = "Mx-150";
    arrsNameFileVariable[1][1] = "Mx2-150_Mx1-1_ctau-1";
    std::vector<TString> namesHistogram({
        "hNumDWanted",
        "hIsDsFoundSignCorrect",
        "hnTHINMatched",
        "hnTHINMatchedUnique",
        "hDeltaRTHINjetPairsFromChi2ordi",
        "hDeltaRTHINjetPairsFromChi2bar",
        "hDeltaRBetweenTwoTHINjetPairs",
        "hTHINjetP4Pt00",
        "hTHINjetPairP4Pt0",
        "hTHINjetP4Eta00",
        "hTHINjetPairP4Eta0",
        "hTHINjetP4Pt00Mismatched"
    });
    // , "hElePairP4MMax", "hElePairP4MtMax", "hMuPairP4MMax", "hMuPairP4MtMax"
    if (DEBUGGING) std::cout << "namesHistogram size: " << namesHistogram.size() << std::endl;
    bool normalize = true;

    TLegend* leg = new TLegend(0.3333,0.7027,0.8333,0.9023);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    for (int iFGroup=0; iFGroup<nFVarableGroups; iFGroup++) {
        TFile *plTFile[nFVarComparing];
        TString nameFileHeadAndVars = nameFileHead;
        for (int iFVar=0; iFVar<nFVarComparing; iFVar++) {
            TString nameFileVariable =  arrsNameFileVariable[iFGroup][iFVar];
            TString nameFile = (TString)"output_" + nameFileHead + "_" + nameFileVariable + "_" + nameFileTail + ".root";
            if (DEBUGGING) std::cout << "Opening file: " << nameFile << std::endl;
            plTFile[iFVar] = new TFile(nameFile);
            if (DEBUGGING) plTFile[iFVar]->ls();
            nameFileHeadAndVars += (TString)"_" + nameFileVariable;

            /*
            gStyle->SetOptFit(1111);
            TString nameHistogramToFit;
            c1->Clear();
            nameHistogramToFit = "hPfMetCorrPt";
            ((TH1F *)plTFile[iFVar]->Get(nameHistogramToFit))->Fit("landau");
            c1->Print((TString)nameOutImageDir + "/" + nameFileHead + "_" + nameFileVariable + "_" +  nameHistogramToFit + "_fit" + ".svg");
            c1->Clear();
            nameHistogramToFit = "hPfMetCorrSig";
            ((TH1F *)plTFile[iFVar]->Get(nameHistogramToFit))->Fit("landau");
            c1->Print((TString)nameOutImageDir + "/" + nameFileHead + "_" + nameFileVariable + "_" +  nameHistogramToFit + "_fit" + ".svg");
            c1->Clear();
            nameHistogramToFit = "hTHINnJet";
            ((TH1F *)plTFile[iFVar]->Get(nameHistogramToFit))->Fit("landau");
            c1->Print((TString)nameOutImageDir + "/" + nameFileHead + "_" + nameFileVariable + "_" +  nameHistogramToFit + "_fit" + ".svg");
            // c1->Clear();
            // FitExp((TH1F *)ptfile->Get("hEleP4M"), c1);
            // c1->Print((TString)nameOutImageDir + "/" + nameFileHead + "_" + nameFileVariable + "_" +  "hEleP4M" + "_fit" + ".svg");
            */
        }
        for (int iHist=0;iHist<(Int_t)(namesHistogram.size());iHist++) {
            TH1 *plTHist[nFVarComparing];
            for (int iFVar=0; iFVar<nFVarComparing; iFVar++) {
                TString nameFileVariable = arrsNameFileVariable[iFGroup][iFVar];
                TString nameHistogram = namesHistogram[iHist];
                if (DEBUGGING) std::cout<< "Extracting " << nameHistogram << " from " << plTFile[iFVar]->GetName() << std::endl;
                plTHist[iFVar] = ((TH1F *)plTFile[iFVar]->Get(nameHistogram));
                if (DEBUGGING) std::cout<< "Extracted " << nameHistogram << std::endl;
                // plTFile[iFVar]->Close(); // DONT CLOSE!! The histogram need the data. Closing cause segment violation when drawing.
                // if (DEBUGGING) std::cout << "Closed file " << plTFile[iFVar]->GetName() << std::endl;
                c1->Clear();
                if (DEBUGGING) std::cout << "Drawing histogram " << plTHist[iFVar]->GetName() << std::endl;
                gStyle->SetOptStat(111111);
                plTHist[iFVar]->Draw();
                TString pathImage = (TString)nameOutImageDir + "/" + nameFileHead + "_" + nameFileVariable + "_" +  nameHistogram + ".svg";
                if (DEBUGGING) std::cout << "Printing canvas to " << pathImage << std::endl;
                c1->Print(pathImage);
            }
            if (DEBUGGING) std::cout << "Decorating histogram " << plTHist[0]->GetName() << " and " << plTHist[1]->GetName() << std::endl;
            plTHist[0]->SetMarkerStyle(8);
            plTHist[0]->SetMarkerSize(1);
            plTHist[0]->SetLineWidth(3);
            plTHist[0]->SetLineColor(4);
            plTHist[1]->SetMarkerStyle(25);
            plTHist[1]->SetMarkerSize(1);
            plTHist[1]->SetLineWidth(3);
            plTHist[1]->SetLineColor(2);
		    for (Int_t iFVar=0;iFVar<nFVarComparing;iFVar++) {
                gStyle->SetOptStat(0);
                plTHist[iFVar]->Sumw2();
                if(normalize){
                    if (DEBUGGING) std::cout << "Normalizing histogram " << plTHist[iFVar]->GetName() << std::endl;
                    plTHist[iFVar]->Scale(1.0/plTHist[iFVar]->Integral());
	            }
            }
            if (DEBUGGING) std::cout << "Setting maximum for histogram " << plTHist[0]->GetName() << " and " << plTHist[1]->GetName() << std::endl;
            float max = TMath::Min(0xFFFFFF.0p0, TMath::Max(plTHist[0]->GetMaximum(), plTHist[1]->GetMaximum()));
            plTHist[0]->SetMaximum(1.1*max);
            plTHist[1]->SetMaximum(1.1*max);

            double chi2=0.;
            int nbins=0;
            if (DEBUGGING) std::cout << "Chi2 comparing histogram " << plTHist[0]->GetName() << " and " << plTHist[1]->GetName() << std::endl;
            chi2NbinsCompare(plTHist[0],plTHist[1],chi2,nbins,1,plTHist[0]->GetNbinsX());

            if (DEBUGGING) std::cout << "Drawing histogram " << plTHist[0]->GetName() << " and " << plTHist[1]->GetName() << std::endl;
            plTHist[0]->Draw("hist");
            plTHist[1]->Draw("histsame");
            leg->Clear();
            leg->SetHeader("");
            leg->AddEntry((TObject*)0, Form("#chi^{2}/NDF=%.1f/%d",chi2,nbins), "");
            leg->AddEntry((TObject*)0, Form("#chi^{2} Prob=%.2f",TMath::Prob(chi2,nbins)), "");
            leg->AddEntry((TObject*)0, "", "");
            // leg->AddEntry(h1,endfix1.Data(),"l");
            // leg->AddEntry(h2,endfix2.Data(),"l");
            leg->Draw("same");
            TString pathImage = (TString)nameOutImageDir + "/" + nameFileHeadAndVars + "_" +  namesHistogram[iHist] + ".svg";
            if (DEBUGGING) std::cout << "Printing canvas to " << pathImage << std::endl;
            c1->Print(pathImage);

        }
    }
}
