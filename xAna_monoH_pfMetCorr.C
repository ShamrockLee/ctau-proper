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

#define DEBUG_LEVEL 1

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
#include <TMath.h>
#include "compare_utils_11.h" // For vectorFirstN
#include <algorithm> // for `std::max_element`
#include <numeric> // for `std::iota`

/*
void efferr(float nsig,float ntotal,float factor=1) {
    float eff = nsig/ntotal;
    float err = sqrt( (1-eff)*eff/ntotal);
    std::cout << "efficiency = " << eff*factor << " +- " << err*factor << std::endl;
}
*/

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

struct {
    Int_t nEle;
    std::vector<bool>* eleIsPass;
    Int_t* eleCharge;
    TClonesArray* eleP4;
    TH1F* hElePairP4MMax = nullptr;
    TH1F* hElePairP4MtMax = nullptr;
    TH1F* hElePairP4MMaxPt = nullptr;
    TH1F* hElePairP4MTwoMaxPt = nullptr;
    bool toDistinguishCharge = false;
} typedef SelectAndFillEPairs;

void selectAndFillElectronPairs(SelectAndFillEPairs selectAndFillEPairs) {
    Int_t nEle = selectAndFillEPairs.nEle;
    std::vector<bool>* eleIsPass = selectAndFillEPairs.eleIsPass;
    Int_t* eleCharge = selectAndFillEPairs.eleCharge;
    TClonesArray* eleP4 = selectAndFillEPairs.eleP4;
    TH1F* hElePairP4MMax = selectAndFillEPairs.hElePairP4MMax;
    TH1F* hElePairP4MtMax = selectAndFillEPairs.hElePairP4MtMax;
    TH1F* hElePairP4MMaxPt = selectAndFillEPairs.hElePairP4MMaxPt;
    TH1F* hElePairP4MWithTwoMaxPt = selectAndFillEPairs.hElePairP4MTwoMaxPt;
    bool toDistinguishCharge = selectAndFillEPairs.toDistinguishCharge;

    Int_t nEleNeg = 0;
    std::vector<Int_t> vectorIElePassed;
    vectorIElePassed.clear();
    Int_t nElePassed = 0;
    // std::vector<bool>& eleIsPassLoose = *((std::vector<bool>*) data.GetPtr("eleIsPassLoose"));
    if (eleIsPass != nullptr) {
        for (Int_t iEle=0; iEle<nEle; iEle++) {
            if (!eleIsPass->at(iEle)) continue;
            vectorIElePassed.push_back(iEle);
            nElePassed++;
        }
        vectorIElePassed.resize(nElePassed);
    } else {
        nElePassed = nEle;
        // vectorIElePassed.resize(nEle);
        // std::iota(vectorIElePassed.begin(), vectorIElePassed.begin() + nEle, 0);
        for (Int_t i=0; i<nEle; vectorIElePassed.push_back(i++));
    }
    if (DEBUG_LEVEL > 1) std::cout << "vectorIElePassed got, nEle: " << nEle << ", nElePassed: " << nElePassed << std::endl;
    if (DEBUG_LEVEL > 1) assert(vectorIElePassed.size() == nElePassed);
    // // inversely sort according to Pt
    // std::sort(vectorIElePassed.begin(), vectorIElePassed.end(), [eleP4](Int_t i, Int_t j) {return ((TLorentzVector*)eleP4->At(i))->Pt() > ((TLorentzVector*)eleP4->At(j))->Pt();});
    
    
    /*
    std::vector<Double_t> vectorEleP4M2(nElePassed);
    for (Int_t i=0; i<nElePassed; i++) {
        vectorEleP4M2[i] = ((TLorentzVector*)eleP4->At(i))->M2();
    }
    */
    if (nElePassed < 2) return;
    TLorentzVector* p4ElePairWithTwoMaxPt = new TLorentzVector();
    std::vector<TLorentzVector> vectorP4ElePairsValid;
    vectorP4ElePairsValid.clear();
    Int_t nElePairsValid = 0;
    if (toDistinguishCharge) {
        std::vector<Int_t> vectorIElesNeg(nElePassed);
        std::vector<Int_t> vectorIElesPos(vectorIElesNeg);
        std::vector<Int_t>::iterator iterIElesNeg = vectorIElesNeg.begin();
        std::vector<Int_t>::iterator iterIElesPos = vectorIElesPos.begin();
        TLorentzVector *eleP4MaxPtNeg, *eleP4MaxPtPos;
        eleP4MaxPtNeg = eleP4MaxPtPos = nullptr;
        for (std::vector<Int_t>::iterator iterIElePassed = vectorIElePassed.begin(); iterIElePassed != vectorIElePassed.end(); iterIElePassed++) {
            TLorentzVector* eleP4Current = (TLorentzVector*)eleP4->At(*iterIElePassed);
            if ((eleCharge[*iterIElePassed]) < 0) {
                // vectorIElesNeg[nEleNeg] = iEle;
                *(iterIElesNeg++) = *iterIElePassed;
                if (eleP4MaxPtNeg == nullptr) {
                    eleP4MaxPtNeg = eleP4Current;
                } else if (eleP4Current->Pt() > eleP4MaxPtNeg->Pt()) {
                    eleP4MaxPtNeg = eleP4Current;
                }
                nEleNeg++;
            } else {
                // vectorIElesPos[iEle - nEleNeg] = iEle;
                *(iterIElesPos++) = *iterIElePassed;
                if (eleP4MaxPtPos == nullptr) {
                    eleP4MaxPtPos = eleP4Current;
                } else if (eleP4Current->Pt() > eleP4MaxPtPos->Pt()) {
                    eleP4MaxPtPos = eleP4Current;
                }
            }
            vectorIElesNeg.resize(nEleNeg);
            vectorIElesPos.resize(nElePassed-nEleNeg);
        }
        if (eleP4MaxPtNeg != nullptr && eleP4MaxPtPos != nullptr) {
            *p4ElePairWithTwoMaxPt = *eleP4MaxPtNeg + *eleP4MaxPtPos;
        } else {
            delete p4ElePairWithTwoMaxPt;
            p4ElePairWithTwoMaxPt = nullptr;
        }
        for (std::vector<Int_t>::iterator iterIElesNeg = vectorIElesNeg.begin(); 
             iterIElesNeg != vectorIElesNeg.end(); iterIElesNeg++) {
            for (std::vector<Int_t>::iterator iterIElesPos = vectorIElesPos.begin(); 
                iterIElesPos != vectorIElesPos.end(); iterIElesPos++) {
                vectorP4ElePairsValid.push_back(*((TLorentzVector*)eleP4->At(*iterIElesNeg)) + *((TLorentzVector*)(eleP4->At(*iterIElesPos))));
                nElePairsValid++;
            }
        }
    } else {
        // auto lambdaMapperPt = ([eleP4](Int_t i){return (Double_t)(((TLorentzVector*)eleP4->At(i))->Pt());});
    // Double_t (*mapperPt)(Int_t) = Lambda::ptr1<Double_t, Int_t>(lambdaMapperPt);
        std::vector<Int_t> vectorITwoElesWithMaxPt = vectorFirstNWithMapper<Int_t, Double_t>(
            vectorIElePassed, 2, false, true, 
            [eleP4](Int_t i){return (Double_t)(((TLorentzVector*)eleP4->At(i))->Pt());}
        );
        *p4ElePairWithTwoMaxPt = *(TLorentzVector*)eleP4->At(vectorITwoElesWithMaxPt[0]) 
         + *(TLorentzVector*)eleP4->At(vectorITwoElesWithMaxPt[1]);
         for (std::vector<Int_t>::iterator iterIElePassedLeft = vectorIElePassed.begin(); 
              iterIElePassedLeft != vectorIElePassed.end(); 
              iterIElePassedLeft++) {
                  for (std::vector<Int_t>::iterator iterIElePassedRight = iterIElePassedLeft + 1;
                      iterIElePassedRight != vectorIElePassed.end();
                      iterIElePassedRight++) {
                      vectorP4ElePairsValid.push_back(*((TLorentzVector*)eleP4->At(*iterIElePassedLeft)) + *((TLorentzVector*)(eleP4->At(*iterIElePassedRight))));
                      nElePairsValid++;
                }
            }
    }
    if (hElePairP4MWithTwoMaxPt != nullptr && p4ElePairWithTwoMaxPt != nullptr) {
        hElePairP4MWithTwoMaxPt->Fill(p4ElePairWithTwoMaxPt->M());
    }
    // auto compareEleM2Inversed = [vectorEleP4M2](Int_t i, Int_t j) {return vectorEleP4M2[i] > vectorEleP4M2[j];};
    // std::sort(vectorIElesNeg.begin(), vectorIElesNeg.end(), compareEleM2Inversed);
    // std::sort(vectorIElesPos.begin(), vectorIElesPos.end(), compareEleM2Inversed);
    vectorP4ElePairsValid.resize(nElePairsValid);
    vectorP4ElePairsValid.shrink_to_fit();
    Double_t p4MPairCurrent, p4MPairMax, p4MPairMaxMt, p4MPairMaxPt, p4MtPairCurrent, p4MtPairMax, p4PtPairCurrent, p4PtPairMax;
    if (nElePairsValid > 0) {
        p4MPairMax = p4MPairMaxMt = p4MPairMaxPt = (vectorP4ElePairsValid.at(0)).M();
        p4MtPairMax = vectorP4ElePairsValid.at(0).Mt();
        p4PtPairMax = vectorP4ElePairsValid.at(0).Pt();
        for(std::vector<TLorentzVector>::iterator iterP4ElePairsValid = vectorP4ElePairsValid.begin(); 
            iterP4ElePairsValid != vectorP4ElePairsValid.end(); 
            iterP4ElePairsValid++) {
            p4MPairCurrent = iterP4ElePairsValid->M();
            p4MtPairCurrent = iterP4ElePairsValid->Mt();
            p4PtPairCurrent = iterP4ElePairsValid->Pt();
            if (p4MPairCurrent > p4MPairMax) {
                p4MPairMax = p4MPairCurrent;
            }
            if (p4MtPairCurrent > p4MtPairMax) {
                p4MtPairMax = p4MtPairCurrent;
                // p4MPairMaxMt = p4MPairCurrent;
            }
            if (p4PtPairCurrent > p4PtPairMax) {
                p4PtPairMax = p4PtPairCurrent;
                p4MPairMaxPt = p4MPairCurrent;
            }
        }
        if (hElePairP4MMax != nullptr) hElePairP4MMax->Fill(p4MPairMax);
        // hElePairP4MMaxMt->Fill(p4MPairMaxMt);
        if (hElePairP4MtMax != nullptr) hElePairP4MtMax->Fill(p4MtPairMax);
        if (hElePairP4MMaxPt != nullptr) hElePairP4MMaxPt->Fill(p4MPairMaxPt);
    }
}

void xAna_monoH_PfMetCorrs(std::string nameInputFile, std::string nameOutputFile, bool toRecreateOutFile=true, bool debug=false) {
    TreeReader data(nameInputFile.data());
    ULong64_t nTotal=0;
    // ULong64_t nPass[20]={0};
    TH1F* hPfMetCorrPt = new TH1F("hPfMetCorrPt", "PfMetCorrPt", 250, 0, 500);
    TH1F* hPfMetCorrPhi = new TH1F("hPfMetCorrPhi", "PfMetCorrPhi", 100, -10, 10);
    TH1F* hPfMetCorrSig = new TH1F("hPfMetCorrSig", "PfMetCorrSig", 100, 0, 100);
    TH1F* hFATnJet = new TH1F("hFATnJet", "FATnJet", 6, 0 - 0.5f, 6 -0.5f);
    TH1F* hTHINnJet = new TH1F("hTHINnJet", "THINnJet", 50, 0 - 0.5f, 50 - 0.5f);
    // TH1C* hNGoodTHINJet = new TH1C("hNGoodTHINJet", "nGoodTHINJet", 50, 0, 50);
    // TH1C* hNGoodTHINBJet = new TH1C("hNGoodTHINBJet", "nGoodTHINBJet", 50, 0, 50);
    TH1F* hElePairP4MMax = new TH1F("hElePairP4MMax", "Max mass of e-pairs", 1000*1, 0, 1000);
    // TH1F* hElePairP4MMaxMt = (TH1F*) hElePairP4MMax->Clone("hElePairP4MMaxMt");
    // hElePairP4MMaxMt->SetTitle("Mass of e-pairs with max Mt");
    TH1F* hElePairP4MMaxPt = (TH1F*) hElePairP4MMax->Clone("hElePairP4MMaxPt");
    hElePairP4MMaxPt->SetTitle("Mass of e-pairs with max Pt");
    TH1F* hElePairP4MTwoMaxPt = (TH1F*) hElePairP4MMax->Clone("hElePairP4MTwoMaxPt");
    hElePairP4MTwoMaxPt->SetTitle("Mass of e-pairs constructed from the two with max Pt");
    TH1F* hElePairP4MtMax = new TH1F("hElePairP4MtMax", "Max m_t of e-pairs", 200, 0, 5000);
    // TH1F* hMuPairP4MMax = new TH1F("hMuPairP4MMax", "Max mass of mu-pairs", 512, 0.2112, 0.2115);
    TH1F* hMuPairP4MMax = new TH1F("hMuPairP4MMax", "Max mass of mu-pairs", 1000*1, 0, 1000);
    TH1F* hMuPairP4MMaxPt = (TH1F*) hMuPairP4MMax->Clone("hMuPairP4MMaxPt");
    hMuPairP4MMaxPt->SetTitle("Mass of mu-pairs with max Pt");
    TH1F* hMuPairP4MTwoMaxPt = (TH1F*) hMuPairP4MMax->Clone("hMuPairP4MTwoMaxPt");
    hMuPairP4MTwoMaxPt->SetTitle("Mass of mu-pairs constructed from the two with max Pt");
    TH1F* hMuPairP4MtMax = new TH1F("hMuPairP4MtMax", "Max m_t of mu-pairs", 200, 0, 1000);
    // TH1F* hMuPairP4MMaxPt = 
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
        
        // TClonesArray* eleCharge = (TClonesArray*) data.GetPtrTObject("eleCharge");

        // Electron pairs
        SelectAndFillEPairs inputParsElectron;
        inputParsElectron.nEle = data.GetInt("nEle");
        inputParsElectron.eleCharge = (Int_t*) data.GetPtrInt("eleCharge");
        inputParsElectron.eleP4 = (TClonesArray*)data.GetPtrTObject("eleP4");
        inputParsElectron.eleIsPass = ((std::vector<bool>*) data.GetPtr("eleIsPassLoose"));
        inputParsElectron.hElePairP4MMax = hElePairP4MMax;
        inputParsElectron.hElePairP4MMaxPt = hElePairP4MMaxPt;
        inputParsElectron.hElePairP4MTwoMaxPt = hElePairP4MTwoMaxPt;
        inputParsElectron.hElePairP4MtMax = hElePairP4MtMax;
        inputParsElectron.toDistinguishCharge = false;
        selectAndFillElectronPairs(inputParsElectron);

        
        // Muon pairs
        SelectAndFillEPairs inputParsMuon;
        inputParsMuon.nEle = data.GetInt("nMu");
        inputParsMuon.eleCharge = (Int_t*) data.GetPtrInt("muCharge");
        inputParsMuon.eleP4 = (TClonesArray*) data.GetPtrTObject("muP4");
        inputParsMuon.eleIsPass = nullptr;//((std::vector<bool>*) data.GetPtr("muIsPassLoose"));
        inputParsMuon.hElePairP4MMax = hMuPairP4MMax;
        inputParsMuon.hElePairP4MMaxPt = hMuPairP4MMaxPt;
        inputParsMuon.hElePairP4MTwoMaxPt = hMuPairP4MTwoMaxPt;
        inputParsMuon.hElePairP4MtMax = hMuPairP4MtMax;
        inputParsMuon.toDistinguishCharge = true;
        selectAndFillElectronPairs(inputParsMuon);

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
    // gStyle->SetOptStat(111111);
    hPfMetCorrPt->Write();
    hPfMetCorrPhi->Write();
    hPfMetCorrSig->Write();
    hFATnJet->Write();
    hTHINnJet->Write();
    //hNGoodTHINJet->Write();
    //hNGoodTHINBJet->Write();
    hElePairP4MMax->Write();
    // hElePairP4MMaxMt->Write();
    hElePairP4MMaxPt->Write();
    hElePairP4MTwoMaxPt->Write();
    hElePairP4MtMax->Write();
    hMuPairP4MMax->Write();
    hMuPairP4MMaxPt->Write();
    hMuPairP4MTwoMaxPt->Write();
    hMuPairP4MtMax->Write();
    //hHPSTau_4MomentumM->Write();
    if (DEBUG_LEVEL > 0) outFile->ls();
    outFile->Close();
}
