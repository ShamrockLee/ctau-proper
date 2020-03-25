/** Compare variables to determine if long-lived particle exists
 * 
 * This file compares 
 * * `pfMetCorrPt`
 * * `pfMetCorrPhi`
 * * `pfMetCorrSig`
 * * Number of AK4 jets passing basic jet IDs
 * * Invariant mass of a pair of electrons (or muons)
 * to differ the ROOT file with long-lived particles 
 * from the other.
 */

#define DEBUGGING true

#include <vector>
#include <iostream>
#include <string>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TStyle.h>

void efferr(float nsig,float ntotal,float factor=1) {
    float eff = nsig/ntotal;
    float err = sqrt( (1-eff)*eff/ntotal);
    std::cout << "efficiency = " << eff*factor << " +- " << err*factor << std::endl;
}

template<typename TypeTArray=TClonesArray, TypeTArray* TArray=nullptr, typename TypeIndex=Int_t>
bool compareTArrayElementsAt(TypeIndex i, TypeIndex j) {
    return TArray->At(i) < TArray->At(j);
}

template<typename TypeArray=std::vector<int>, TypeArray array=nullptr, typename TypeIndex=Int_t>
bool compareArrayElementsAt(TypeIndex i, TypeIndex j) {
    return array[i] < array[j];
}

bool fcnIsM2CloseEnough(Double_t M2_a, Double_t M2_b) {
    return TMath::Max(M2_a, M2_b) <= TMath::Min(M2_a, M2_b) * 0x1.004p0;
}

void selectAndFillElectronPairs(Int_t nEle, Int_t* eleCharge, TClonesArray* eleP4, TH1F* hElePairP4MMax, TH1F* hElePairP4MtMax) {
    Int_t nEleNeg = 0;
    std::vector<Int_t> vectorIElesNeg(nEle);
    std::vector<Int_t> vectorIElesPos(vectorIElesNeg);
    for (Int_t iEle=0; iEle<nEle; iEle++) {
        if ((eleCharge[iEle]) < 0) {
            vectorIElesNeg[nEleNeg] = iEle;
            nEleNeg++;
        } else {
            vectorIElesPos[iEle - nEleNeg] = iEle;
        }
        vectorIElesNeg.resize(nEleNeg);
        vectorIElesPos.resize(nEle-nEleNeg);
    }

    std::vector<Double_t> vectorEleP4M2(nEle);
    for (Int_t i=0; i<nEle; i++) {
        vectorEleP4M2[i] = ((TLorentzVector*)eleP4->At(i))->M2();
    }
    // auto compareM2 = [eleP4](Int_t i, Int_t j){return ((TLorentzVector*)eleP4->At(i))->M2() < ((TLorentzVector*)eleP4->At(j))->M2();};
    // auto compareM2 = [vectorEleP4M2](Int_t i, Int_t j, bool (*comparor)(Int_t, Int_t) = nullptr) {
    //     return (comparor == nullptr) ? (*comparor)(vectorEleP4M2[i], vectorEleP4M2[j]) : vectorEleP4M2[i] < vectorEleP4M2[j];};
    auto compareEleM2Inversed = [vectorEleP4M2](Int_t i, Int_t j) {return vectorEleP4M2[i] > vectorEleP4M2[j];};
    std::sort(vectorIElesNeg.begin(), vectorIElesNeg.end(), compareEleM2Inversed);
    std::sort(vectorIElesPos.begin(), vectorIElesPos.end(), compareEleM2Inversed);
    // std::vector<Int_t*> vectorElePairsValid(nEle <= nEleNeg*2 ? nEleNeg : nEle - nEleNeg);
    std::vector<TLorentzVector> vectorP4ElePairsValid(nEle <= nEleNeg*2 ? nEleNeg : nEle - nEleNeg);
    // for (Int_t jEleNeg=0, jElePos=0; jEleNeg<nEleNeg && jElePos<nEle-nEleNeg; (vectorEleP4M2.at(jEleNeg) <= vectorEleP4M2.at(jElePos)) ? jEleNeg++ : jElePos++)
    //     if (fcnIsM2CloseEnough(vectorEleP4M2.at(jEleNeg), vectorEleP4M2.at(jElePos))) vectorElePairsValid
    Int_t jEleNeg = 0, jElePos = 0, jElePairValid = 0;
    // Int_t **currentIElePairs;
    while (jEleNeg < nEleNeg && jElePos < nEle - nEleNeg) {
        if (fcnIsM2CloseEnough(vectorEleP4M2.at(jEleNeg), vectorEleP4M2.at(jElePos))) {
            // currentIElePairs = new Int_t*;
            // (*currentIElePairs)[0] = vectorIElesNeg[jEleNeg];
            // (*currentIElePairs)[1] = vectorIElesPos[jElePos];
            // vectorElePairsValid[jElePairValid] = *currentIElePairs;
            vectorP4ElePairsValid[jElePairValid] = *((TLorentzVector*)eleP4->At(vectorIElesNeg[jEleNeg])) + *((TLorentzVector*)(eleP4->At(vectorIElesPos[jElePos])));
            jEleNeg++, jElePos++, jElePairValid++;
        } else if (jEleNeg <= jElePos) jEleNeg++; else jElePos++;
    }
    vectorP4ElePairsValid.resize(jElePairValid);
    vectorP4ElePairsValid.shrink_to_fit();
    Double_t p4MPairCurrent, p4MPairMax, p4MPairMaxMt, p4MtPairCurrent, p4MtPairMax;
    if (jElePairValid > 0) {
        p4MPairMax = p4MPairMaxMt = (vectorP4ElePairsValid.at(0)).M();
        p4MtPairMax = vectorP4ElePairsValid.at(0).Mt();
        for(Int_t iP4 = 0; iP4 < jElePairValid; iP4++) {
            p4MPairCurrent = vectorP4ElePairsValid[iP4].M();
            p4MtPairCurrent = vectorP4ElePairsValid[iP4].Mt();
            if (p4MPairCurrent > p4MPairMax) {
                p4MPairMax = p4MPairCurrent;
            }
            if (p4MtPairCurrent > p4MtPairMax) {
                p4MtPairMax = p4MtPairCurrent;
                // p4MPairMaxMt = p4MPairCurrent;
            }
        }
        hElePairP4MMax->Fill(p4MPairMax);
        // hElePairP4MMaxMt->Fill(p4MPairMaxMt);
        hElePairP4MtMax->Fill(p4MtPairMax);
    }
}

void xAna_monoH_PfMetCorrs(std::string nameInputFile, std::string nameOutputFile, bool toRecreateOutFile=true, bool debug=false) {
    TreeReader data(nameInputFile.data());
    ULong64_t nTotal=0;
    // ULong64_t nPass[20]={0};
    TH1F* hPfMetCorrPt = new TH1F("hPfMetCorrPt", "PfMetCorrPt", 512, 0, 512);
    TH1F* hPfMetCorrPhi = new TH1F("hPfMetCorrPhi", "PfMetCorrPhi", 1024, -8, 8);
    TH1F* hPfMetCorrSig = new TH1F("hPfMetCorrSig", "PfMetCorrSig", 512, 0, 125);
    TH1F* hFATnJet = new TH1F("hFATnJet", "FATnJet", 6, 0 - 0.5f, 6 -0.5f);
    TH1F* hTHINnJet = new TH1F("hTHINnJet", "THINnJet", 64, 0 - 0.5f, 64 - 0.5f);
    // TH1C* hNGoodTHINJet = new TH1C("hNGoodTHINJet", "nGoodTHINJet", 64, 0, 64);
    // TH1C* hNGoodTHINBJet = new TH1C("hNGoodTHINBJet", "nGoodTHINBJet", 64, 0, 64);
    TH1F* hElePairP4MMax = new TH1F("hEleP4MMax", "Max mass of e-pairs", 512, 0, 0.5);
    // TH1F* hElePairP4MMaxMt = (TH1F*) hElePairP4MMax->Clone("hEleP4MMaxMt");
    // hElePairP4MMaxMt->SetTitle("Mass of e-pairs with max Mt");
    TH1F* hElePairP4MtMax = new TH1F("hEleP4MtMax", "Max m_t of e-pairs", 512, 0, 0.5);
    TH1F* hMuPairP4MMax = new TH1F("hMuP4MMax", "Max mass of mu-pairs", 512, 0.2112, 0.2115);
    TH1F* hMuPairP4MtMax = new TH1F("hMuP4MtMax", "Max m_t of mu-pairs", 512, 0, 0.5);
    // TH1F* hHPSTau_4MomentumM = new TH1F("hHPSTau_4MomentumM", "HPSTau_4Momentum->M()", 512, 0, 0.25);
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++) {
        /*
        if (jEntry % 50000 == 0)
            fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
        nPass[0]++;
        */
        data.GetEntry(jEntry);
        nTotal ++;
        /*
        //0. has a good vertex
        int nVtx        = data.GetInt("nVtx");
        if(nVtx<1)continue;
        nPass[0]++;
        */

        /*
        //1. trigger 
        std::string* trigName = data.GetPtrString("hlt_trigName");
        std::vector<bool> &trigResult = *((std::vector<bool>*) data.GetPtr("hlt_trigResult"));

        bool passTrigger=false;
        for(unsigned int it=0; it< trigResult.size(); it++) {
            std::string thisTrig= trigName[it];
            bool results = trigResult[it];

            if( (thisTrig.find("HLT_PFMET90_PFMHT90_")!= 
                std::string::npos && results==1) || 
                (thisTrig.find("HLT_PFMET170_NoiseCleaned")!= 
                std::string::npos && results==1)) {
                // cout << thisTrig << endl;
                passTrigger=true;
                break;
            }
        }

        if(!passTrigger)continue;
        nPass[1]++;

        // apply noise filters if it is data
        bool isData = data.GetBool("isData");
        std::string* filterName = data.GetPtrString("hlt_filterName");
        std::vector<bool> &filterResult = *((std::vector<bool>*) data.GetPtr("hlt_filterResult"));
        bool passFilter=false;
        for(unsigned int it=0; it< filterResult.size(); it++){
            std::string thisFilter= filterName[it];
            bool results = filterResult[it];

            if( (thisFilter.find("Flag_CSCTightHaloFilter")!= 
                std::string::npos && results==1) &&
                (thisFilter.find("Flag_eeBadScFilter")!= 
                std::string::npos && results==1) &&
                (thisFilter.find("Flag_HBHENoiseFilter")!= 
                std::string::npos && results==1) &&
                (thisFilter.find("Flag_HBHENoiseIsoFilter")!= 
                std::string::npos && results==1) ) {
                passFilter=true;
                break;
            }	
        }
        if( isData && !passFilter )continue;
        nPass[2]++;
        */
        Double_t pfMetCorrPt = data.GetFloat("pfMetCorrPt");
        Double_t pfMetCorrPhi = data.GetFloat("pfMetCorrPhi");
        /*
        if(pfMetCorrPt<170.)continue;
        nPass[3]++;
        */
        /*
        // (NOT) // veto extra electrons
        int nEle = data.GetInt("nEle");
        TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");
        std::vector<bool>& eleIsPassLoose = *((std::vector<bool>*) data.GetPtr("eleIsPassLoose"));
        std::vector<int> myEles;
        myEles.clear();
        for(int ie = 0; ie < nEle; ie++) {    
            TLorentzVector* myEle = (TLorentzVector*)eleP4->At(ie);
            if( myEle->Pt()<10 )continue;
            if( fabs(myEle->Eta())>2.5 )continue;
            if( !eleIsPassLoose[ie] )continue;
            myEles.push_back(ie);
            hEleP4M->Fill(myEle->M());
        } 
        // if(myEles.size()>0)continue;
        nPass[4]++;
        */
        /*
        // (NOT) //veto extra muons
        int nMu = data.GetInt("nMu");
        TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");
        std::vector<bool>& isLooseMuon = *((std::vector<bool>*) data.GetPtr("isLooseMuon"));
        float* muChHadIso = data.GetPtrFloat("muChHadIso");
        float* muNeHadIso = data.GetPtrFloat("muNeHadIso");
        float* muGamIso   = data.GetPtrFloat("muGamIso");
        float* muPUPt     = data.GetPtrFloat("muPUPt");
        */
        
        /*
        std::vector<int> myMuos;
        myMuos.clear();
        for(int im = 0; im < nMu; im++) {
        TLorentzVector* myMu = (TLorentzVector*)muP4->At(im);
        if( myMu->Pt()<10 )continue;
        if( fabs(myMu->Eta())>2.4 ) continue;
        if( !isLooseMuon[im] )continue;
        
        float relPFIso = (muChHadIso[im]+ 
                TMath::Max(0., muNeHadIso[im] + muGamIso[im] - 0.5*muPUPt[im]))/myMu->Pt();
        if(relPFIso>0.4)continue;
        myMuos.push_back(im);
        hMuP4M->Fill(myMu->M());
        }
        // if(myMuos.size()>0)continue;
        nPass[5]++;
        */
        
        /*
        // (NOT) //veto extra taus
        int nTau = data.GetInt("HPSTau_n");
        TClonesArray* tauP4 = (TClonesArray*) data.GetPtrTObject("HPSTau_4Momentum");
        std::vector<bool>& isDecayModeFinding = *((std::vector<bool>*) data.GetPtr("disc_decayModeFinding"));
        std::vector<bool>& passLooseTauIso = *((std::vector<bool>*) data.GetPtr("disc_byLooseIsolationMVA3oldDMwLT"));
    
        std::vector<int> myTaus;
        for(int it=0; it < nTau; it++) {
        TLorentzVector* myTau = (TLorentzVector*)tauP4->At(it);
        if( myTau->Pt()<20 )continue;
        if( fabs(myTau->Eta())>2.3 )continue;
        if( !isDecayModeFinding[it] )continue;
        if( !passLooseTauIso[it] )continue;
        myTaus.push_back(it);
        hHPSTau_4MomentumM->Fill(myTau->M());
        }
        //if(myTaus.size()>0)continue;
        nPass[6]++;
        */
        
        /*
        //find a pair of b-jets that could be a Higgs candidate
        const int nTHINJets     = data.GetInt("THINnJet");
        TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
        float* thinJetCSV =  data.GetPtrFloat("THINjetCISVV2");
        std::vector<bool>& passThinJetLooseID = *((std::vector<bool>*) data.GetPtr("THINjetPassIDLoose"));
        std::vector<bool>& passThinJetPUID = *((std::vector<bool>*) data.GetPtr("THINisPUJetIDTight"));    
        
        float maxHpt=-999;
        int Hindex[2]={-1,-1};

        for(int ij=0; ij < nTHINJets; ij++){
            TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
            if(thisJet->Pt()<30)continue;
            if(fabs(thisJet->Eta())>2.4)continue;
            if(!passThinJetLooseID[ij])continue;
            if(!passThinJetPUID[ij])continue;
            
            // for b-jet (medium ID)
            if(thinJetCSV[ij]<0.80)continue;

            for(int jj=0; jj < ij ; jj++){
                TLorentzVector* thatJet = (TLorentzVector*)thinjetP4->At(jj);
                if(thatJet->Pt()<30)continue;
                if(fabs(thatJet->Eta())>2.4)continue;
                if(!passThinJetLooseID[jj])continue;
                if(!passThinJetPUID[jj])continue;
                
                // for b-jet (medium ID)
                if(thinJetCSV[jj]<0.80)continue;

                float thisHpt = (*thisJet + *thatJet).Pt();
                float thisHmass = (*thisJet + *thatJet).M();
                if(thisHpt < 150.)continue;

                if(thisHmass < 100)continue;
                if(thisHmass > 150)continue;
                
                if(thisHpt>maxHpt){
                    Hindex[0]=jj;
                    Hindex[1]=ij;
                    maxHpt   =thisHpt;
                }
            } // end of inner loop jet
        } // end of outer loop jet
        */

        /*
        if(Hindex[0]<0 || Hindex[1]<0)continue;
        nPass[7]++;
        */

        /*
        TLorentzVector  bjet[2];
        for(int ib=0; ib<2;ib++)bjet[ib] = 
                    *((TLorentzVector*)thinjetP4->At(Hindex[ib]));
        TLorentzVector  higgsJet = bjet[0]+bjet[1];

        // cout << "Hindex [0] = " << Hindex[0] << endl;
        // cout << "Hindex [1] = " << Hindex[1] << endl;

        // (NOT) //veto n>=1 AK4 jets and extra AK4 b jet
        
        unsigned int nGoodTHINBJets=0;
        unsigned int nGoodTHINJets=0;
        int jetIndex=-1; // extra AK4 jets if there is one
        std::vector<int> indexForDPhi;
        indexForDPhi.clear();

        for(int ij=0; ij < nTHINJets; ij++){
        TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
        if(thisJet->Pt()<30)continue;
        if(fabs(thisJet->Eta())>4.5)continue;
        if(!passThinJetLooseID[ij])continue;
        if(!passThinJetPUID[ij])continue;

        indexForDPhi.push_back(ij);

        if(thisJet->DeltaR(bjet[0])<0.4)continue;
        if(thisJet->DeltaR(bjet[1])<0.4)continue;


        nGoodTHINJets++;
        jetIndex=ij;

        // for b-jet
        if(fabs(thisJet->Eta())>2.4)continue;
        if(thinJetCSV[ij]<0.46)continue;
        nGoodTHINBJets++;

        } // end of loop
        */
        
        Int_t nGenPar = data.GetInt("nGenPar");
        Int_t nEle = data.GetInt("nEle");
        Int_t* eleCharge = (Int_t*) data.GetPtrInt("eleCharge");
        TClonesArray* eleP4 = (TClonesArray*)data.GetPtrTObject("eleP4");
        // TClonesArray* eleCharge = (TClonesArray*) data.GetPtrTObject("eleCharge");

        // Electron pairs

        selectAndFillElectronPairs(nEle, eleCharge, eleP4, hElePairP4MMax, hElePairP4MtMax);

        Int_t nMu = data.GetInt("nMu");
        Int_t* muCharge = (Int_t*) data.GetPtrInt("muCharge");
        TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");
        
        // Muon pairs
        selectAndFillElectronPairs(nMu, muCharge, muP4, hMuPairP4MMax, hMuPairP4MtMax);

        hPfMetCorrPt->Fill(pfMetCorrPt);
        hPfMetCorrPhi->Fill(pfMetCorrPhi);
        hPfMetCorrSig->Fill(data.GetFloat("pfMetCorrSig"));
        // hNGoodTHINJet->Fill(nGoodTHINJets);
        // hNGoodTHINBJet->Fill(nGoodTHINBJets);
        hFATnJet->Fill(data.GetInt("FATnJet"));
        hTHINnJet->Fill(data.GetInt("THINnJet"));
        // hEleP4M->Fill(((TLorentzVector *)data.GetPtrTObject("eleP4"))->M());
        // hMuP4M->Fill(((TLorentzVector *)data.GetPtrTObject("muP4"))->M());
    }
    TFile* outFile = new TFile(nameOutputFile.data(), toRecreateOutFile ? "recreate" : "update");
    gStyle->SetOptStat(111111);
    hPfMetCorrPt->Write();
    hPfMetCorrPhi->Write();
    hPfMetCorrSig->Write();
    hFATnJet->Write();
    hTHINnJet->Write();
    //hNGoodTHINJet->Write();
    //hNGoodTHINBJet->Write();
    hElePairP4MMax->Write();
    // hElePairP4MMaxMt->Write();
    hElePairP4MtMax->Write();
    hMuPairP4MMax->Write();
    hMuPairP4MtMax->Write();
    //hHPSTau_4MomentumM->Write();
    if (DEBUGGING) outFile->ls();
    outFile->Close();
}
