#include <iostream>
#include <vector>
#include <algorithm>
#include <Rtypes.h>

typedef Int_t TypeE;

std::vector<TypeE> vec123({1, 2, 3});

void printVector(std::vector<TypeE>& vec, std::string sep=", ") {
  std::cout << "{";
  typename std::vector<TypeE>::iterator iter = vec.begin();
  if (vec.size() > 0) {
    std::cout << *iter;
    iter++;
  }
  for(;iter != vec.end(); iter++) {
    std::cout << sep << *iter;
  }
  std::cout << "}" << std::endl;
}

void testNextPermutation(std::vector<TypeE>& vec) {
  printVector(vec);
  std::cout << std::next_permutation(vec.begin(), vec.end()) << std::endl;
  printVector(vec);;
}

void testDoWhileNextPermutation(std::vector<TypeE>& vec) {
  std::vector<TypeE> vecSwapped;
  vecSwapped = vec;
  std::sort(vecSwapped.begin(), vecSwapped.end());
  do {
    printVector(vecSwapped);
  } while (std::next_permutation(vecSwapped.begin(), vecSwapped.end()));
}
