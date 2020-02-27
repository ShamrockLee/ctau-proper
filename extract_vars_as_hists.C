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



void extractVarsAsTH1s(const std::string inputFile, const std::string arrNamesVariable[], const TreeReader::ETypes arrKtypes[], TObjArray *arrHists, const Int_t nHists, bool debug=false) {
    //get TTree from InputFile
    TreeReader data(inputFile.data());
    Long64_t nTotal=0;
    Long64_t nPass[20]={0};
    // TH1F* hOutput = new TH1F((TString)nameVariable,"",100000,0,10000);
    for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
        data.GetEntry(jEntry);
        nTotal ++;
        Long64_t nPass[20] = {0};
        Int_t nGenPar = data.GetInt("nGenPar");
        for(int ig=0; ig < nGenPar; ig++) {
            for(int ihist = 0; ihist < nHists; ihist++) {
                if (arrKtypes[ihist] == TreeReader::ETypes::kBool) {
                    ((TH1 *)((*arrHists)[ihist]))->Fill(data.GetBool((TString)arrNamesVariable[ihist]));
                } else if (arrKtypes[ihist] == TreeReader::ETypes::kChar) {
                    ((TH1 *)((*arrHists)[ihist]))->Fill(data.GetChar((TString)arrNamesVariable[ihist]));
                } else if (arrKtypes[ihist] == TreeReader::ETypes::kShort) {
                    ((TH1 *)((*arrHists)[ihist]))->Fill(data.GetShort((TString)arrNamesVariable[ihist]));
                } else if (arrKtypes[ihist] == TreeReader::ETypes::kInt) {
                    ((TH1 *)((*arrHists)[ihist]))->Fill(data.GetInt((TString)arrNamesVariable[ihist]));
                } else if (arrKtypes[ihist] == TreeReader::ETypes::kFloat) {
                    ((TH1 *)((*arrHists)[ihist]))->Fill(data.GetFloat((TString)arrNamesVariable[ihist]));
                } else if (arrKtypes[ihist] == TreeReader::ETypes::kDouble) {
                    ((TH1 *)((*arrHists)[ihist]))->Fill(data.GetDouble((TString)arrNamesVariable[ihist]));
                } else {
                }
            }
        }
    }
}

