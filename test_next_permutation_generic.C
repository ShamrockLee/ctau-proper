#include <iostream>
#include <vector>
#include <algorithm>
#include <Rtypes.h>
std::vector<int> vec123({1, 2, 3});

template<typename E>
void printVector(std::vector<E>& vec, std::string sep=", ") {
  std::cout << "{";
  typename std::vector<E>::iterator iter = vec.begin();
  if (vec.size() > 0) {
    std::cout << *iter;
    iter++;
  }
  for(;iter != vec.end(); iter++) {
    std::cout << sep << *iter;
  }
  std::cout << "}" << std::endl;
}

template<typename E>
void testNextPermutation(std::vector<E>& vec) {
  printVector(vec);
  std::cout << std::next_permutation(vec.begin(), vec.end()) << std::endl;
  printVector(vec);;
}

template<typename E>
void testDoWhileNextPermutation(std::vector<E>& vec) {
  std::vector<E> vecSwapped(vec);
  std::sort(vecSwapped.begin(), vecSwapped.end());
  do {
    printVector(vecSwapped);
  } while (std::next_permutation(vecSwapped.begin(), vecSwapped.end()));
}
