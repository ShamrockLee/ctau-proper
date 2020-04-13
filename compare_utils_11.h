#ifndef COMPARATOR_UTILS_H
#define COMPARATOR_UTILS_H

#define COMPARATOR_UTILS_H_DEBUG true

#include <vector>
#include <functional> // for `std::function` and `std::less`, available for c++11
template <typename E>
class LinkedListUnit {
public:
    E element;
    LinkedListUnit* ptrPrevious;
    LinkedListUnit* ptrNext;
LinkedListUnit(E element=nullptr, LinkedListUnit* ptrPrevious=nullptr, LinkedListUnit* ptrNext=nullptr) {
    this->element = element;
    this->ptrPrevious = ptrPrevious;
    this->ptrNext = ptrNext;
}

};

template <typename E>
std::vector<E> vectorFirstN(std::vector<E>& vectorToCompare, int nRanks, bool isUnique, bool isGreatest, std::function<bool (E, E)> comparator) {
    typedef typename std::vector<E>::iterator iteratorE;
    std::vector<E> vectorOutput(nRanks);
    if (vectorToCompare.size() == 0) {
        return vectorOutput;
    }
    bool isNotFull = true;
    iteratorE iterVectorToCompare = vectorToCompare.begin();
    LinkedListUnit<E>* ptrFirst = new LinkedListUnit<E>(*(iterVectorToCompare++));
    LinkedListUnit<E>* ptrLast = ptrFirst;
    int nRanksCurrent=1;
    while (iterVectorToCompare != vectorToCompare.end()) {
        LinkedListUnit<E>* ptrCurrent = ptrFirst;
        while (true) {
            if (isUnique && *iterVectorToCompare==ptrCurrent->element) break;
            if (isGreatest xor comparator(*iterVectorToCompare, ptrCurrent->element)) {
                LinkedListUnit<E>* ptrInserted = new LinkedListUnit<E>(*iterVectorToCompare, ptrCurrent->ptrPrevious, ptrCurrent);
                    if (ptrCurrent->ptrPrevious != nullptr) {
                        ptrCurrent->ptrPrevious->ptrNext = ptrInserted;
                    } else {
                        ptrFirst = ptrInserted;
                    }
                    ptrCurrent->ptrPrevious = ptrInserted;
                    // ptrCurrent=ptrCurrent->ptrNext;
                    if (isNotFull) {
                        nRanksCurrent++;
                        isNotFull = (nRanksCurrent < nRanks);
                    } else {
                        ptrLast = ptrLast->ptrPrevious;
                        delete ptrLast->ptrNext;
                        ptrLast->ptrNext = nullptr;
                    }
                    break;
            }
            if (ptrCurrent == ptrLast) {
                if (isNotFull && !(isUnique && *iterVectorToCompare==ptrLast->element)) {
                    ptrLast->ptrNext = new LinkedListUnit<E>(*iterVectorToCompare, ptrLast);
                    ptrLast->ptrNext->ptrPrevious = ptrLast;
                    ptrLast = ptrLast->ptrNext;
                    nRanksCurrent++;
                    isNotFull = (nRanksCurrent < nRanks);
                }
                break;
            }
            ptrCurrent = ptrCurrent->ptrNext;
        }
        iterVectorToCompare++;
    }
    iteratorE iterVectorOutput = vectorOutput.begin();
    LinkedListUnit<E>* ptrCurrent = ptrFirst;
    for (int i=0; i<nRanks-1; i++) {
        *iterVectorOutput = ptrCurrent->element;
        iterVectorOutput++;
        ptrCurrent=ptrCurrent->ptrNext;
        delete ptrCurrent->ptrPrevious;
        ptrCurrent->ptrPrevious = nullptr;
    }
    *iterVectorOutput = ptrCurrent->element;
    delete ptrCurrent->ptrPrevious;
    delete ptrCurrent;
    return vectorOutput;
}

template<typename E, typename M, typename Tmapper = std::function<M(E)>>
std::vector<E> vectorFirstNWithMapper(std::vector<E>& vectorToCompare, int nRanks, bool isUnique, bool isGreatest, Tmapper mapper) {
    std::function<M(E)> mapperCasted = std::function<M(E)>(mapper);
    auto lambdaMapper = [mapperCasted](E var1, E var2){return mapperCasted(var1) < mapperCasted(var2);};
    return vectorFirstN(vectorToCompare, nRanks, isUnique, isGreatest, std::function<bool(E, E)>(lambdaMapper));
}

template <typename E>
std::vector<E> vectorFirstN(std::vector<E>& vectorToCompare, int nRanks, bool isUnique=false, bool isGreatest=false) {
    return vectorFirstN(vectorToCompare, nRanks, isUnique, isGreatest, static_cast<std::function<bool(E, E)>>(std::less<E>()));
}

#endif
