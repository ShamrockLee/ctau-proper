#ifndef BRANCHMOUNTER_H
#define BRANCHMOUNTER_H

#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <TTree.h>

#include <functional>
#include <iterator>
#include <vector>


class BranchDescription {
 public:
  TString name;
  TString title;
  TString typeName;
  BranchDescription() {}
  BranchDescription(const char *name, const char *title, const char *typeName) {
    this->name = name;
    this->title = title;
    this->typeName = typeName;
  }
  explicit BranchDescription(TString name, TString title, TString typeName) {
    this->name = name;
    this->title = title;
    this->typeName = typeName;
  }
};

template <typename TypeSize = size_t>
class VirtualBranchMounter {
 public:
  virtual Bool_t BranchFor(BranchDescription description) = 0;
  virtual void PrepareForBranching() = 0;
  virtual void BranchOn(TTree *treeOut) = 0;
  virtual void PrepareForPushing(TypeSize n) = 0;
  virtual void PushOrSkip(Bool_t isAcceptable) = 0;
};

namespace BranchMounterHelper {
template <class BranchMounter, typename TypeSize = size_t>
void RunOnFunNext(BranchMounter &mounter, std::function<Bool_t()> funNext,
                  TypeSize n) {
  mounter.PrepareForPushing(n);
  for (TypeSize i = 0; i < n; i++) {
    mounter.PushOrSkip(funNext());
  }
}
template <class BranchMounter, typename TypeSize = size_t>
void RunOnFunAt(BranchMounter &mounter, std::function<Bool_t(TypeSize i)> funAt,
                TypeSize n) {
  mounter.PrepareForPushing(n);
  for (TypeSize i = 0; i < n; i++) {
    mounter.PushOrSkip(funAt(i));
  }
}
template <class BranchMounter, typename TypeIterBool,
          typename TypeSize = size_t>
void RunForEachN(BranchMounter &mounter, TypeIterBool begin, TypeSize n) {
  mounter.PrepareForPushing(n);
  std::for_each(begin, begin + n, mounter.PushOrSkip);
}
template <class BranchMounter, typename TypeIterBool,
          typename TypeSize = size_t>
void RunForEachN(BranchMounter &mounter, TypeIterBool begin, TypeIterBool end) {
  mounter.PrepareForPushing(std::distance(begin, end));
  std::for_each(begin, end, mounter.PushOrSkip);
}
}  // namespace BranchMounterHelper

template <typename E, typename TypeIterE=E*>
class VirtualBranchMounterVector : public VirtualBranchMounter<size_t> {
  typedef size_t TypeSize;
 public:
  virtual Bool_t GetIsAcceptable(BranchDescription description) = 0;
  virtual Bool_t BranchFor(BranchDescription description) {
    if (GetIsAcceptable(description)) {
      vecDescription.push_back(description);
      return true;
    }
    return false;
  }
  virtual void PrepareForBranching() {
    vecIterE.reserve(vecDescription.size());
  }
  virtual void BranchOn(TTree *treeOut) {
    for (auto description: vecDescription) {
      treeOut->Branch(description.name, &vecE)->SetTitle(description.title);
    }
  }
  virtual void SetFunIterE(std::function<TypeIterE (BranchDescription description)> funIterE) {
    this->funIterE = funIterE;
  }
  virtual void PrepareForPushing(TypeSize n) {
    vecIterE.clear();
    for (auto description: vecDescription){
      vecIterE.emplace_back(funIterE(description));
    }
  }
  virtual void PushOrSkip(Bool_t isAcceptable) {
    for (auto iterE: vecIterE) {
      if (isAcceptable) vecE.emplace_back(static_cast<E>(*iterE));
      iterE++;
    }
  }
  VirtualBranchMounterVector() {
    vecDescription.clear();
  }
  VirtualBranchMounterVector(std::function<TypeIterE (BranchDescription description)> funIterE):
  VirtualBranchMounterVector() {
    this->SetFunIterE(funIterE);
  }
 protected:
  std::vector<BranchDescription> vecDescription;
  std::function<TypeIterE (BranchDescription description)> funIterE;
  std::vector<TypeIterE> vecIterE;
  std::vector<E> vecE;
};

template<typename TypeIterE=LongDouble_t*>
class BranchMounterLongDouble : public VirtualBranchMounterVector<LongDouble_t, TypeIterE> {
  typedef LongDouble_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return description.typeName.Contains("LongDouble") ||
           description.typeName.Contains("long double") ||
           description.typeName.Contains("Float128");
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description)
  } 
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Double_t*>
class BranchMounterDouble : public VirtualBranchMounterVector<Double_t, TypeIterE> {
  typedef Double_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return (description.typeName.Contains("Double") ||
           description.typeName.Contains("double") ||
           description.typeName.Contains("Float64")) &&
           !BranchMounterLongDouble::GetIsAcceptableStatic(description);
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description)
  } 
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Float16_t*>
class BranchMounterFloat16 : public VirtualBranchMounterVector<Float16_t, TypeIterE> {
  typedef Float16_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return description.typeName.Contains("Float16") ||
           description.typeName.Contains("float16");
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Float_t*>
class BranchMounterFloat : public VirtualBranchMounterVector<Float_t, TypeIterE> {
  typedef Float_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return (description.typeName.Contains("Float") ||
           description.typeName.Contains("float")) &&
           !BranchMounterDouble::GetIsAcceptableStatic(description) &&
           !BranchMounterFloat16::GetIsAcceptableStatic(description);
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=UInt_t*>
class BranchMounterUInt : public VirtualBranchMounterVector<UInt_t, TypeIterE> {
  typedef UInt_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return description.typeName.Contains("UInt") ||
           description.typeName.Contains("unsigned int");
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Long64_t*>
class BranchMounterLong64 : public VirtualBranchMounterVector<Long64_t, TypeIterE> {
  typedef Long64_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return description.typeName.Contains("Long64") ||
           description.typeName.Contains("long long") ||
           description.typeName.Contains("Int64");
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Long_t*>
class BranchMounterLong : public VirtualBranchMounterVector<Long_t, TypeIterE> {
  typedef Long_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return (description.typeName.Contains("Long") ||
           description.typeName.Contains("long")) &&
           !BranchMounterLong64::GetIsAcceptableStatic(description);
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Short_t*>
class BranchMounterShort : public VirtualBranchMounterVector<Short_t, TypeIterE> {
  typedef Short_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return description.typeName.Contains("Short") ||
           description.typeName.Contains("short") ||
           description.typeName.Contains("Int16");
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Char_t*>
class BranchMounterChar : public VirtualBranchMounterVector<Char_t, TypeIterE> {
  typedef Char_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return description.typeName.Contains("Char") ||
           description.typeName.Contains("char");
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Int_t*>
class BranchMounterInt : public VirtualBranchMounterVector<Int_t, TypeIterE> {
  typedef Int_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return (description.typeName.Contains("Int") ||
           description.typeName.Contains("int")) &&
           !BranchMounterUInt::GetIsAcceptableStatic(description) &&
           !BranchMounterLong64::GetIsAcceptableStatic(description) &&
           !BranchMounterShort::GetIsAcceptableStatic(description);
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};

template<typename TypeIterE=Bool_t*>
class BranchMounterBool : public VirtualBranchMounterVector<Bool_t, TypeIterE> {
  typedef Bool_t E;
 public:
  static inline Bool_t GetIsAcceptableStatic(BranchDescription description) {
    return description.typeName.Contains("Bool") ||
           description.typeName.Contains("bool");
  }
  Bool_t GetIsAcceptable(BranchDescription description) {
    return GetIsAcceptableStatic(description);
  }
  using VirtualBranchMounterVector<E, TypeIterE>::VirtualBranchMounterVector;
};
#endif