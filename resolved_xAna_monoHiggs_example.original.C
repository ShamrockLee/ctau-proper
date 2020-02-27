// example code to run 2015 mono-Higgs resolved selections on signal (EXO-16-012)

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

void efferr(float nsig,float ntotal,float factor=1)
{
  float eff = nsig/ntotal;
  float err = sqrt( (1-eff)*eff/ntotal);
  cout << "efficiency = " << eff*factor << " +- " << err*factor << endl;

}


using namespace std;
void resolved_xAna_monoHiggs(std::string inputFile){

  //get TTree from file ...
  TreeReader data(inputFile.data());

  Long64_t nTotal=0;
  Long64_t nPass[20]={0};


  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal++;
    //0. has a good vertex
    int nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[0]++;

    //1. trigger 
    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));

    bool passTrigger=false;
    for(unsigned int it=0; it< trigResult.size(); it++)
      {
	std::string thisTrig= trigName[it];
	bool results = trigResult[it];

	if( (thisTrig.find("HLT_PFMET90_PFMHT90_")!= 
	     std::string::npos && results==1) || 
	    (thisTrig.find("HLT_PFMET170_NoiseCleaned")!= 
	     std::string::npos && results==1))
	  {
	    //	    cout << thisTrig << endl;
	    passTrigger=true;
	    break;
	  }


      }

    if(!passTrigger)continue;
    nPass[1]++;

    // apply noise filters if it is data
    bool isData = data.GetBool("isData");
    std::string* filterName = data.GetPtrString("hlt_filterName");
    vector<bool> &filterResult = *((vector<bool>*) data.GetPtr("hlt_filterResult"));
    bool passFilter=false;
    for(unsigned int it=0; it< filterResult.size(); it++)
      {
	std::string thisFilter= filterName[it];
	bool results = filterResult[it];

	if( (thisFilter.find("Flag_CSCTightHaloFilter")!= 
	     std::string::npos && results==1) &&
	    (thisFilter.find("Flag_eeBadScFilter")!= 
	     std::string::npos && results==1) &&
	    (thisFilter.find("Flag_HBHENoiseFilter")!= 
	     std::string::npos && results==1) &&
	    (thisFilter.find("Flag_HBHENoiseIsoFilter")!= 
	     std::string::npos && results==1) )
	  {
	    passFilter=true;
	    break;
	  }	
      }
    if( isData && !passFilter )continue;
    nPass[2]++;
	



    float pfMet = data.GetFloat("pfMetCorrPt");
    float pfMetPhi = data.GetFloat("pfMetCorrPhi");
    
    if(pfMet<170.)continue;
    nPass[3]++;



    //veto extra electrons
    int    nEle       = data.GetInt("nEle");
    TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");
    vector<bool>& eleIsPassLoose = *((vector<bool>*) data.GetPtr("eleIsPassLoose"));
    vector<int> myEles;
    myEles.clear();
    for(int ie = 0; ie < nEle; ie++){    
      TLorentzVector* myEle = (TLorentzVector*)eleP4->At(ie);
      if( myEle->Pt()<10 )continue;
      if( fabs(myEle->Eta())>2.5 )continue;
      if( !eleIsPassLoose[ie] )continue;
      myEles.push_back(ie);
    } 
    if(myEles.size()>0)continue;
    nPass[4]++;

    //veto extra muons
    int    nMu       = data.GetInt("nMu");
    TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");
    vector<bool>& isLooseMuon = *((vector<bool>*) data.GetPtr("isLooseMuon"));
    float* muChHadIso = data.GetPtrFloat("muChHadIso");
    float* muNeHadIso = data.GetPtrFloat("muNeHadIso");
    float* muGamIso   = data.GetPtrFloat("muGamIso");
    float* muPUPt     = data.GetPtrFloat("muPUPt");
       
    vector<int> myMuos;
    myMuos.clear();
    for(int im = 0; im < nMu; im++){
      TLorentzVector* myMu = (TLorentzVector*)muP4->At(im);
      if( myMu->Pt()<10 )continue;
      if( fabs(myMu->Eta())>2.4 ) continue;
      if( !isLooseMuon[im] )continue;
      
      float relPFIso = (muChHadIso[im]+ 
			TMath::Max(0., muNeHadIso[im] + muGamIso[im] - 0.5*muPUPt[im]))/myMu->Pt();
      if(relPFIso>0.4)continue;
      myMuos.push_back(im);
    }
    if(myMuos.size()>0)continue;
    nPass[5]++;

    //veto extra taus
    int    nTau      = data.GetInt("HPSTau_n");
    TClonesArray* tauP4 = (TClonesArray*) data.GetPtrTObject("HPSTau_4Momentum");
    vector<bool>& isDecayModeFinding = *((vector<bool>*) data.GetPtr("disc_decayModeFinding"));
    vector<bool>& passLooseTauIso = *((vector<bool>*) data.GetPtr("disc_byLooseIsolationMVA3oldDMwLT"));
   
    vector<int> myTaus;
    for(int it=0; it < nTau; it++)
      {
	TLorentzVector* myTau = (TLorentzVector*)tauP4->At(it);
	if( myTau->Pt()<20 )continue;
	if( fabs(myTau->Eta())>2.3 )continue;
	if( !isDecayModeFinding[it] )continue;
	if( !passLooseTauIso[it] )continue;
	myTaus.push_back(it);
      }
    if(myTaus.size()>0)continue;
    nPass[6]++;
				      

    //find a pair of b-jets that could be a Higgs candidate
    const int nTHINJets     = data.GetInt("THINnJet");
    TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
    float* thinJetCSV =  data.GetPtrFloat("THINjetCISVV2");
    vector<bool>& passThinJetLooseID = *((vector<bool>*) data.GetPtr("THINjetPassIDLoose"));
    vector<bool>& passThinJetPUID = *((vector<bool>*) data.GetPtr("THINisPUJetIDTight"));    
    
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
	
	if(thisHpt>maxHpt)
	  {
	    Hindex[0]=jj;
	    Hindex[1]=ij;
	    maxHpt   =thisHpt;
	  }

      } // end of inner loop jet


    } // end of outer loop jet


    if(Hindex[0]<0 || Hindex[1]<0)continue;
    nPass[7]++;

    TLorentzVector  bjet[2];
    for(int ib=0; ib<2;ib++)bjet[ib] = 
			      *((TLorentzVector*)thinjetP4->At(Hindex[ib]));
    TLorentzVector  higgsJet = bjet[0]+bjet[1];

    // cout << "Hindex [0] = " << Hindex[0] << endl;
    // cout << "Hindex [1] = " << Hindex[1] << endl;

    //veto n>=1 AK4 jets and extra AK4 b jet
    
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

    if(nGoodTHINJets>1)continue;
    nPass[8]++;

    if(nGoodTHINBJets>0)continue;
    nPass[9]++;


    bool passDphi=true;
    for(unsigned int i=0; i<indexForDPhi.size(); i++)
      {
	int jetIndex=indexForDPhi[i];
	TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(jetIndex);
	double dphi=TVector2::Phi_mpi_pi(pfMetPhi-thisJet->Phi());

	if(fabs(dphi)<0.4)
	  {
	    passDphi=false;
	    break;
	  }	
      }
    
    if(!passDphi)continue;
    nPass[10]++;
    

  } // end of loop over entries


  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;


  efferr(nPass[10],nTotal);

  // TFile* outFile = new TFile("test.root","recreate");

  // h_hh->Write();
  // h_hh2->Write();
  // h_SD1->Write();
  // h_SD2->Write();

  // outFile->Close();



}

