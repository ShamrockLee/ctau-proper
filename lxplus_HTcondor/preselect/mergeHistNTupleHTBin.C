#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <array>

#include <TH1.h>
#include <TFile.h>
#include <TKey.h>
#include <TString.h>

void mergeHistNTupleHTBin(const char* const clusterId) {
  constexpr size_t nRanges = 8;
  const std::array<std::string, nRanges> aStrHTRanges {"70to100", "100to200", "200to400", "400to600", "600to800", "800to1200", "1200to2500", "2500toInf"};
  const std::array<double, nRanges> aXSecHT {175300., 147400., 41040., 5674., 1358., 622.9, 151.2, 3.659};
  // /eos/user/y/yuehshun/ntuple_filelist_outputmerged.tmp/2016BkgMC/DY/DYJets_HTBin/DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_merged_2806993.root
  const std::string outputMergedSpace = "/eos/user/y/yuehshun/ntuple_filelist_outputmerged.tmp";
  const std::string subDirRootHT = "2016BkgMC/DY/DYJets_HTBin";
  const std::string outputPath = Form(
      "%s/%s/DYJetsToLL_M-50_HT-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_allmerged.root",
  TFile *tfOut = TFile::Open(outputPath.c_str(), "recreate");
  for (size_t i = 0; i < nRanges; ++i) {
    const std::string inputPath = Form(
        "%s/%s/DYJetsToLL_M-50_HT-%s_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_merged_%s.root",
        outputSpace, subDirRootHT, strHTRange, clusterId);
        outputSpace, subDirRootHT, strHTRange, clusterId);
    TFile *tfIn = TFile::Open(inputPath.c_str());
    if (! i) {
      tfOut->cd();
    }
  }
}

int main(int argc, char* argv[]) {
  return 0;
}
