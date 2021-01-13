#ifndef BRANCH_MOUNTER_H
#define BRANCH_MOUNTER_H

#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <TTree.h>

#include <array>
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

class VirtualBranchMounterBase {
 public:
  virtual void BranchOn(TTree *treeOut) {};
  virtual ~VirtualBranchMounterBase() {};
};
class VirtualBranchMounterScalar : public VirtualBranchMounterBase {
 public:
  virtual void Push() = 0;
};
template <typename TypeSize = size_t>
class VirtualBranchMounterIterable : public VirtualBranchMounterBase {
 public:
  virtual void PrepareForPushing(const TypeSize n) = 0;
  virtual void PushOrSkip(const Bool_t isAcceptable) = 0;
};

class VirtualDescriptionCollector {
 public:
  virtual Bool_t BranchFor(BranchDescription description) = 0;
  virtual void Prepare() {};
  virtual ~VirtualDescriptionCollector() {};
};

template<typename E>
class VirtualBranchMounterScalarSingle: public VirtualBranchMounterScalar {
 public:
  void SetFunPush(std::function<E(BranchDescription description)> funPush) {
    this->funPush = funPush;
  }
  // ~VirtualBranchMounterScalarSingle() {};
 protected:
  std::function<E(BranchDescription description)> funPush;
};

template<typename E, typename TypeIterEIn = E*>
class VirtualBranchMounterIterableSingle
    : public VirtualBranchMounterIterable<size_t> {
 protected:
  std::function<TypeIterEIn(BranchDescription description)> funIterEIn;
 public:
  void SetFunIterEIn(std::function<TypeIterEIn(BranchDescription description)> funIterEIn) {
    this->funIterEIn = funIterEIn;
  }
};

class VirtualDescriptionCollectorSingle: public VirtualDescriptionCollector{
 public:
  std::vector<BranchDescription> vDescriptions;
  virtual Bool_t BranchFor(BranchDescription description) {
    if (GetIsAcceptable(description)){
      vDescriptions.emplace_back(description);
      return true;
    }
    return false;
  }
  virtual void Prepare() {
    vDescriptions.shrink_to_fit();
  }
  VirtualDescriptionCollectorSingle() {
    vDescriptions.clear();
  }
 protected:
  virtual Bool_t GetIsAcceptable(BranchDescription description)=0;
};

template <typename E>
// typedef Int_t E
class BranchMounterScalarSingle
    : public VirtualBranchMounterScalarSingle<E> {
 protected:
  typedef size_t TypeSize;
  const std::vector<BranchDescription> vDescriptions;
  TypeSize nDescriptions;
  std::function<E(BranchDescription)> funPush;
  std::function<BranchDescription(BranchDescription)> funDescriptionModified;
 public:
  // typedef VirtualDescriptionCollectorSingle Setting;
  typedef std::vector<E> Binding;
  Binding vE;
  void SetFunPush(std::function<E(BranchDescription)> funPush) {
    this->funPush = funPush;
  }
  void BranchOn(TTree *treeOut) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = funDescriptionModified == nullptr ? vDescriptions[iDescription] : funDescriptionModified(vDescriptions[iDescription]);
      treeOut->Branch(description.name, static_cast<E*>(&(vE[iDescription])))->SetTitle(description.title);
    }
  }
  void Push() {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      vE[iDescription] = funPush(vDescriptions[iDescription]);
    }
  }
  BranchMounterScalarSingle(const std::vector<BranchDescription> vDescriptions,
  std::function<E(BranchDescription)> funPush = nullptr,
  std::function<BranchDescription(BranchDescription)> funDescriptionModified = nullptr)
  : vDescriptions(vDescriptions), nDescriptions(vDescriptions.size()), funPush(funPush), funDescriptionModified(funDescriptionModified) {
    vE.resize(nDescriptions);
  }
  BranchMounterScalarSingle(const VirtualDescriptionCollectorSingle &setting,
  std::function<E(BranchDescription)> funPush = nullptr,
  std::function<BranchDescription(BranchDescription)> funDescriptionModified = nullptr)
  : BranchMounterScalarSingle(setting.vDescriptions, funPush, funDescriptionModified) {}
};

template<>
class BranchMounterScalarSingle<Bool_t>: public BranchMounterScalarSingle<Char_t> {
 protected:
  using BranchMounterScalarSingle<Char_t>::TypeSize;
  using BranchMounterScalarSingle<Char_t>::vDescriptions;
 public:
  using BranchMounterScalarSingle<Char_t>::vE;
  // using BranchMounterScalarSingle<Char_t>::Binding;
  void BranchOn(TTree *treeOut) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = funDescriptionModified == nullptr ? vDescriptions[iDescription] : funDescriptionModified(vDescriptions[iDescription]);
      treeOut->Branch(description.name, (Bool_t*)&(vE[iDescription]))->SetTitle(description.title);
    }
  }
  void SetFunPush(std::function<Bool_t(BranchDescription)> funPush) {
    BranchMounterScalarSingle<Char_t>::SetFunPush([funPush](BranchDescription description)->Char_t{
      return static_cast<Char_t>(funPush(description));
    });
  }
  BranchMounterScalarSingle(const std::vector<BranchDescription> vDescriptions,
    std::function<Bool_t(BranchDescription)> funPush = nullptr,
    std::function<BranchDescription(BranchDescription)> funDescriptionModified = nullptr)
  : BranchMounterScalarSingle<Char_t>(vDescriptions, funPush, funDescriptionModified) {}
  BranchMounterScalarSingle(const VirtualDescriptionCollectorSingle &setting,
    std::function<Bool_t(BranchDescription)> funPush = nullptr,
    std::function<BranchDescription(BranchDescription)> funDescriptionModified = nullptr)
  : BranchMounterScalarSingle<Char_t>(setting, funPush, funDescriptionModified) {}
};

template <typename E, typename TypeIterEIn>
// typedef Int_t E;
// typedef Int_t* TypeIterEIn;
class BranchMounterVectorSingle: public VirtualBranchMounterIterableSingle<E, TypeIterEIn> {
 protected:
  typedef size_t TypeSize;
  const std::vector<BranchDescription> vDescriptions;
  const size_t nDescriptions;
  std::function<TypeIterEIn(BranchDescription)> funIterEIn;
  std::vector<TypeIterEIn> vIterEIn;
  std::function<BranchDescription(BranchDescription)> funDescriptionModified;
 public:
  // typedef VirtualDescriptionCollectorSingle Setting;
  typedef std::vector<std::vector<E>> Binding;
  Binding vvE;
  void BranchOn(TTree *treeOut) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = funDescriptionModified == nullptr ? vDescriptions[iDescription] : funDescriptionModified(vDescriptions[iDescription]);
      treeOut->Branch(description.name, &(vvE[iDescription]))->SetTitle(description.title);
    }
  }
  void SetFunIterEIn(std::function<TypeIterEIn(BranchDescription)> funIterEIn) {
    this->funIterEIn = funIterEIn;
  }
  void PrepareForPushing(TypeSize n) {
    vIterEIn.clear();
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      vvE[iDescription].clear();
      vvE[iDescription].reserve(n);
      vIterEIn.emplace_back(funIterEIn(vDescriptions[iDescription]));
    }
  }
  void PushOrSkip(const Bool_t isAcceptable) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      if (isAcceptable) {
        vvE[iDescription].emplace_back(*(vIterEIn[iDescription]));
      }
      vIterEIn[iDescription]++;
    }
  }
  BranchMounterVectorSingle(const std::vector<BranchDescription> vDescriptions,
    std::function<TypeIterEIn(BranchDescription)> funIterEIn = nullptr,
    std::function<BranchDescription(BranchDescription)> funDescriptionModified = nullptr)
  : vDescriptions(vDescriptions), nDescriptions(vDescriptions.size()), funIterEIn(funIterEIn), funDescriptionModified(funDescriptionModified) {
    vvE.resize(nDescriptions);
    vIterEIn.reserve(nDescriptions);
  }
  BranchMounterVectorSingle(const VirtualDescriptionCollectorSingle &setting,
    std::function<TypeIterEIn(BranchDescription)> funIterEIn = nullptr,
    std::function<BranchDescription(BranchDescription)> funDescriptionModified = nullptr)
  : BranchMounterVectorSingle(setting.vDescriptions, funIterEIn, funDescriptionModified) {}
};

class DescriptionCollectorFromFun: public VirtualDescriptionCollectorSingle {
 public:
  DescriptionCollectorFromFun(std::function<Bool_t(BranchDescription description)> funIsAcceptable): 
    funIsAcceptable(funIsAcceptable){
  }
 protected:
  const std::function<Bool_t(BranchDescription description)> funIsAcceptable;
  Bool_t GetIsAcceptable(BranchDescription description) {
    return funIsAcceptable(description);
  }
};

template<class ClassCollector=VirtualDescriptionCollector, Bool_t deletePointers=false>
class DescriptionCollectorChained: public VirtualDescriptionCollector {
 public:
  std::vector<ClassCollector *> vCollectors;
  virtual Bool_t BranchFor(BranchDescription description) {
    for (ClassCollector *collector: vCollectors) {
      if (collector->BranchFor(description)) {
        return true;
      }
    }
    return false;
  }
  virtual void Prepare() {
    for (ClassCollector *collector: vCollectors) {
      collector->Prepare();
    }
  }
  DescriptionCollectorChained(const std::vector<ClassCollector *> vCollectors):
    vCollectors(vCollectors) {
    \
  }
  virtual ~DescriptionCollectorChained() {
    if (deletePointers) {
      for (ClassCollector *collector: vCollectors) {
        delete collector;
      }
    }
  }
};

template<class ClassMounter=VirtualBranchMounterBase, Bool_t deletePointers=false>
class BranchMounterBaseChained: public VirtualBranchMounterBase {
 public:
  std::vector<ClassMounter *> vMounters;
  virtual void BranchOn(TTree *treeOut) {
    for (ClassMounter *mounter: vMounters) {
      mounter->BranchOn(treeOut);
    }
  }
  BranchMounterBaseChained(const std::vector<ClassMounter *> vMounters)
  : vMounters(vMounters) {
    \
  }
  virtual ~BranchMounterBaseChained() {
    if (deletePointers) {
      for (ClassMounter *mounter: vMounters) {
        delete mounter;
      }
    }
  }
};

template<class ClassMounter=VirtualBranchMounterScalar, Bool_t deletePointers=false>
class BranchMounterScalarChained:
  public BranchMounterBaseChained<ClassMounter, deletePointers>, public VirtualBranchMounterScalar {
 public:
  // std::vector<ClassMounter *> vMounters;
  using BranchMounterBaseChained<ClassMounter, deletePointers>::vMounters;
  virtual void BranchOn(TTree *treeOut) {
    BranchMounterBaseChained<ClassMounter, deletePointers>::BranchOn(treeOut);
  }
  virtual void Push() {
    for (ClassMounter *mounter: vMounters) {
      mounter->Push();
    }
  }
  using BranchMounterBaseChained<ClassMounter, deletePointers>::BranchMounterBaseChained;
};

template<class ClassMounter=VirtualBranchMounterIterable<size_t>, typename TypeSize=size_t, Bool_t deletePointers=false>
class BranchMounterIterableChained:
  public BranchMounterBaseChained<ClassMounter, deletePointers>, public VirtualBranchMounterIterable<TypeSize> {
 public:
  // std::vector<ClassMounter *> vMounters;
  using BranchMounterBaseChained<ClassMounter, deletePointers>::vMounters;
  virtual void BranchOn(TTree *treeOut) {
    BranchMounterBaseChained<ClassMounter, deletePointers>::BranchOn(treeOut);
  }
  virtual void PrepareForPushing(TypeSize n) {
    for (ClassMounter *mounter: vMounters) {
      mounter->PrepareForPushing(n);
    }
  }
  virtual void PushOrSkip(Bool_t isAccptable) {
    for (ClassMounter *mounter: vMounters) {
      mounter->PushOrSkip(isAccptable);
    }
  }
  using BranchMounterBaseChained<ClassMounter, deletePointers>::BranchMounterBaseChained;
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
  for (TypeIterBool iterBool = begin; iterBool != begin + n; iterBool++) {
    mounter.PushOrSkip(static_cast<Bool_t>(*iterBool));
  }
}
template <class BranchMounter, typename TypeIterBool,
          typename TypeSize = size_t>
void RunForEach(BranchMounter &mounter, TypeIterBool begin, TypeIterBool end) {
  mounter.PrepareForPushing(std::distance(begin, end));
  for (TypeIterBool iterBool = begin; iterBool != end; iterBool++) {
    mounter.PushOrSkip(static_cast<Bool_t>(*iterBool));
  }
}
template <class BranchMounter, class TypeIterIdx, typename TypeSize = size_t>
TypeSize RunForEachNIdxSorted(BranchMounter &mounter, TypeIterIdx begin, TypeSize n) {
  mounter.PrepareForPushing(n);
  TypeIterIdx iterIdxToPass = begin;
  TypeSize ii = 0;
  TypeSize idxIn = 0;
  while (ii < n) {
    TypeSize idxToPass = static_cast<TypeSize>(*(iterIdxToPass++));
    // std::cout << "ii: " << ii << " idxToPass: " << idxToPass << '\n';
    for (; idxIn < idxToPass; idxIn++) {
      mounter.PushOrSkip(false);
      // std::cout << "idxIn: " << idxIn << ", Skip one." << '\n';
    }
    // std::cout << "idxIn: " << idxIn << ", Push one." << '\n';
    mounter.PushOrSkip(true);
    idxIn++;
    ii++;
  }
  // std::cout << std::flush;
  return ii;
}
template <class BranchMounter, class TypeIterIdx, typename TypeSize = size_t>
TypeSize RunForEachIdxSorted(BranchMounter &mounter, TypeIterIdx begin, TypeIterIdx end) {
  return RunForEachNIdxSorted(mounter, begin, std::distance(begin, end));
}
}  // namespace BranchMounterHelper

namespace CommonAcceptability {
Bool_t GetIsAcceptableLongDouble(BranchDescription description) {
  return description.typeName.Contains("LongDouble") ||
         description.typeName.Contains("long double") ||
         description.typeName.Contains("Float128");
}
Bool_t GetIsAcceptableDouble(BranchDescription description) {
  return (description.typeName.Contains("Double") ||
          description.typeName.Contains("double") ||
          description.typeName.Contains("Float64")) &&
         !GetIsAcceptableLongDouble(description);
}
Bool_t GetIsAcceptableFloat16(BranchDescription description) {
  return description.typeName.Contains("Float16") ||
         description.typeName.Contains("float16");
}
Bool_t GetIsAcceptableFloat(BranchDescription description) {
  return (description.typeName.Contains("Float") ||
          description.typeName.Contains("float")) &&
         !GetIsAcceptableLongDouble(description) &&
         !GetIsAcceptableDouble(description) &&
         !GetIsAcceptableFloat16(description);
}
Bool_t GetIsAcceptableULong64(BranchDescription description) {
  return description.typeName.Contains("ULong64") ||
         description.typeName.Contains("unsigned long long");
}
Bool_t GetIsAcceptableLong64(BranchDescription description) {
  return (description.typeName.Contains("Long64") ||
         description.typeName.Contains("long long")) &&
         !GetIsAcceptableULong64(description);
}
Bool_t GetIsAcceptableULong(BranchDescription description) {
  return (description.typeName.Contains("ULong") ||
         description.typeName.Contains("unsigned long") ||
         description.typeName.Contains("ulong")) &&
         !GetIsAcceptableULong64(description) &&
         !description.typeName.Contains("double");
}
Bool_t GetIsAcceptableLong(BranchDescription description) {
  return (description.typeName.Contains("Long") ||
          description.typeName.Contains("long")) &&
         !GetIsAcceptableULong64(description) &&
         !GetIsAcceptableLong64(description) &&
         !GetIsAcceptableULong(description) &&
         !description.typeName.Contains("double");
}
Bool_t GetIsAcceptableUInt(BranchDescription description) {
  return description.typeName.Contains("UInt") ||
         description.typeName.Contains("unsigned int") ||
         description.typeName.Contains("uint");
}
Bool_t GetIsAcceptableInt(BranchDescription description) {
  return (description.typeName.Contains("Int") ||
          description.typeName.Contains("int")) &&
         !GetIsAcceptableUInt(description);
}
Bool_t GetIsAcceptableUShort(BranchDescription description) {
  return description.typeName.Contains("UShort") ||
         description.typeName.Contains("unsigned short") ||
         description.typeName.Contains("ushort");
}
Bool_t GetIsAcceptableShort(BranchDescription description) {
  return (description.typeName.Contains("Short") ||
         description.typeName.Contains("short")) &&
         !GetIsAcceptableUShort(description);
}
Bool_t GetIsAcceptableUChar(BranchDescription description) {
  return description.typeName.Contains("UChar") ||
         description.typeName.Contains("Byte") ||
         description.typeName.Contains("unsigned char");
}
Bool_t GetIsAcceptableChar(BranchDescription description) {
  return (description.typeName.Contains("Char") ||
         description.typeName.Contains("char")) &&
         !GetIsAcceptableUChar(description);
}
Bool_t GetIsAcceptableBool(BranchDescription description) {
  return description.typeName.Contains("Bool") ||
         description.typeName.Contains("bool");
}
}  // namespace CommonAcceptability

#if false
class DescriptionCollectorDouble: public VirtualDescriptionCollectorSingle {
 protected:
  // Bool_t GetIsAcceptable(BranchDescription description) {
  //   return CommonAcceptability::GetIsAcceptableDouble(description);
  // }
};

class BranchMounterScalarDouble: public BranchMounterScalarSingle<Double_t> {
 public:
  using BranchMounterScalarSingle<Double_t>::BranchOn;
//  protected:
//   Bool_t GetIsAccdptable(BranchDescription description) {
//     return CommonAcceptability::GetIsAcceptableDouble(description);
//   }
};

template <typename TypeIterEIn = Double_t*>
class BranchMounterVectorDouble: public BranchMounterVectorSingle<Double_t, TypeIterEIn> {
 public:
  using BranchMounterVectorSingle<Double_t, TypeIterEIn>::BranchOn;
 protected:
  Bool_t GetIsAccdptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableDouble(description);
  }
};

class DescriptionCollectorBool: public DescriptionCollectorFromFun {
 public:
  DescriptionCollectorBool()
  : DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableBool) {}
};

typedef BranchMounterScalarSingle<Bool_t> BranchMounterScalarBool;

class BranchMounterScalarChainedTest: public BranchMounterScalarChained<VirtualBranchMounterScalar> {
 public:
  BranchMounterScalarChainedTest()
  : BranchMounterScalarChained({(VirtualBranchMounterScalar*) new BranchMounterScalarSingle<Double_t>(std::vector<BranchDescription>())}) {}
};
#endif
#endif
