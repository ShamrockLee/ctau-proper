#include <Rtypes.h>
#include <TTree.h>
#include <TBranch.h>

#include <iostream>

template<size_t n>
TBranch *BranchOnArrayElement(TTree *tree, Int_t (&testArray)[n]) {
  return tree->Branch("testBranchOnArrayElement", &(testArray[0]));
}

TBranch *BranchOnVectorElement(TTree *tree, std::vector<Int_t> &testVector) {
  return tree->Branch("testBranchOnVectorElement", &(testVector[0]));
}

TBranch *BranchOnSingleVar(TTree *tree, Int_t &testInt) {
  return tree->Branch("testBranchOnSingleVar", &testInt);
}

void test_branching() {
  TTree *tree = new TTree;
  Int_t testInt = 0;
  std::vector<Int_t> testVector = {4};
  Int_t testArray[] = {5};
  std::cout << BranchOnSingleVar(tree, testInt) << std::endl;
  std::cout << BranchOnVectorElement(tree, testVector) << std::endl;
  std::cout << BranchOnArrayElement(tree, testArray) << std::endl;
}
