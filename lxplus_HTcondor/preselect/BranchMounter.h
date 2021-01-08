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
};

template<typename E>
class VirtualBranchMounterScalarSingle: public VirtualBranchMounterScalar {
 public:
  void SetFunPush(std::function<E(BranchDescription description)> funPush) {
    this->funPush = funPush;
  }
 protected:
  std::function<E(BranchDescription description)> funPush;
};

template<typename E, typename TypeIterEIn = E*>
class VirtualBranchMounterIterableSingle
    : public VirtualBranchMounterIterable<size_t> {
 public:
  void SetFunIterEIn(std::function<TypeIterEIn(BranchDescription description)> funIterEIn) {
    this->funIterEIn = funIterEIn;
  }
 protected:
  std::function<TypeIterEIn(BranchDescription description)> funIterEIn;
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
  typedef size_t TypeSize;
 public:
  // typedef VirtualDescriptionCollectorSingle Setting;
  typedef std::vector<E> Binding;
  Binding &vE;
  BranchMounterScalarSingle(Binding &binding, const std::vector<BranchDescription> vDescriptions):
    vE(binding), vDescriptions(vDescriptions), nDescriptions(vDescriptions.size()) {
    vE.resize(nDescriptions);
  }
  BranchMounterScalarSingle(Binding &binding, const VirtualDescriptionCollectorSingle *const setting):
    BranchMounterScalarSingle(binding, setting->vDescriptions) {}
  BranchMounterScalarSingle(Binding &binding, const std::vector<BranchDescription> vDescriptions,
  std::function<E(BranchDescription)> funPush):
    BranchMounterScalarSingle(binding, vDescriptions) {
    SetFunPush(funPush);
  }
  BranchMounterScalarSingle(Binding &binding, const VirtualDescriptionCollectorSingle *const setting,
  std::function<E(BranchDescription)> funPush):
    BranchMounterScalarSingle(binding, setting->vDescriptions, funPush) {}
  void SetFunPush(std::function<E(BranchDescription)> funPush) {
    this->funPush = funPush;
  }
  void BranchOn(TTree *treeOut) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = vDescriptions[iDescription];
      treeOut->Branch(description.name, &(vE[iDescription]))->SetTitle(description.title);
    }
  }
  void BranchOn(TTree *treeOut, Binding &vE) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = vDescriptions[iDescription];
      treeOut->Branch(description.name, &(vE[iDescription]))->SetTitle(description.title);
    }
  }
  void Push() {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      vE[iDescription] = funPush(vDescriptions[iDescription]);
    }
  }
 protected:
  const std::vector<BranchDescription> vDescriptions;
  TypeSize nDescriptions;
  std::function<E(BranchDescription)> funPush;
};

template <typename E, typename TypeIterEIn>
// typedef Int_t E;
// typedef Int_t* TypeIterEIn;
class BranchMounterVectorSingle: public VirtualBranchMounterIterableSingle<E, TypeIterEIn> {
  typedef size_t TypeSize;
 public:
  // typedef VirtualDescriptionCollectorSingle Setting;
  typedef std::vector<std::vector<E>> Binding;
  Binding &vvE;
  BranchMounterVectorSingle(Binding &binding, const std::vector<BranchDescription> vDescriptions):
    vvE(binding), vDescriptions(vDescriptions), nDescriptions(vDescriptions.size()) {
    vvE.resize(nDescriptions);
    vIterEIn.reserve(nDescriptions);
  }
  BranchMounterVectorSingle(Binding &binding, const VirtualDescriptionCollectorSingle *const setting):
    BranchMounterVectorSingle(binding, setting->vDescriptions) {}
  BranchMounterVectorSingle(Binding &binding, const std::vector<BranchDescription> vDescriptions,
  std::function<TypeIterEIn(BranchDescription)> funIterEIn):
    BranchMounterVectorSingle(binding, vDescriptions) {
    SetFunIterEIn(funIterEIn);
  }
  BranchMounterVectorSingle(Binding &binding, const VirtualDescriptionCollectorSingle *const setting,
  std::function<TypeIterEIn(BranchDescription)> funIterEIn):
    BranchMounterVectorSingle(binding, setting->vDescriptions, funIterEIn) {}
  void SetFunIterEIn(std::function<TypeIterEIn(BranchDescription)> funIterEIn) {
    this->funIterEIn = funIterEIn;
  }
  void BranchOn(TTree *treeOut) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = vDescriptions[iDescription];
      treeOut->Branch(description.name, &(vvE[iDescription]))->SetTitle(description.title);
    }
  }
  void BranchOn(TTree *treeOut, Binding &vvE) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = vDescriptions[iDescription];
      treeOut->Branch(description.name, &(vvE[iDescription]))->SetTitle(description.title);
    }
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
 protected:
  const std::vector<BranchDescription> vDescriptions;
  const size_t nDescriptions;
  std::function<TypeIterEIn(BranchDescription)> funIterEIn;
  std::vector<TypeIterEIn> vIterEIn;
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

template<class ClassMounter=VirtualDescriptionCollector>
class DescriptionCollectorChained: public VirtualDescriptionCollector {
 public:
  std::vector<ClassMounter *> vCollectors;
  virtual Bool_t BranchFor(BranchDescription description) {
    for (auto collector: vCollectors) {
      if (collector->BranchFor(description)) {
        return true;
      }
    }
    return false;
  }
  virtual void Prepare() {
    for (auto collector: vCollectors) {
      collector->Prepare();
    }
  }
  DescriptionCollectorChained(const std::vector<ClassMounter *> vCollectors):
    vCollectors(vCollectors) {
    \
  }
};

template<class ClassMounter=VirtualBranchMounterBase>
class BranchMounterBaseChained: public VirtualBranchMounterBase {
 public:
  std::vector<ClassMounter *> vMounters;
  virtual void BranchOn(TTree *treeOut) {
    for (ClassMounter *mounter: vMounters) {
      mounter->BranchOn(treeOut);
    }
  }
  BranchMounterBaseChained(const std::vector<ClassMounter *> vMounters):
    vMounters(vMounters) {
    \
  }
};

template<class ClassMounter=VirtualBranchMounterScalar>
class BranchMounterScalarChained:
  public BranchMounterBaseChained<ClassMounter>, public VirtualBranchMounterScalar {
 public:
  // std::vector<ClassMounter *> vMounters;
  using BranchMounterBaseChained<ClassMounter>::vMounters;
  virtual void BranchOn(TTree *treeOut) {
    BranchMounterBaseChained<ClassMounter>::BranchOn(treeOut);
  }
  virtual void Push() {
    for (ClassMounter *mounter: vMounters) {
      mounter->Push();
    }
  }
  using BranchMounterBaseChained<ClassMounter>::BranchMounterBaseChained;
};

template<class ClassMounter=VirtualBranchMounterIterable<size_t>, typename TypeSize=size_t>
class BranchMounterIterableChained:
  public BranchMounterBaseChained<ClassMounter>, public VirtualBranchMounterIterable<TypeSize> {
 public:
  // std::vector<ClassMounter *> vMounters;
  using BranchMounterBaseChained<ClassMounter>::vMounters;
  virtual void BranchOn(TTree *treeOut) {
    BranchMounterBaseChained<ClassMounter>::BranchOn(treeOut);
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
  using BranchMounterBaseChained<ClassMounter>::BranchMounterBaseChained;
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
void RunForEachN(BranchMounter &mounter, TypeIterBool begin, TypeIterBool end) {
  mounter.PrepareForPushing(std::distance(begin, end));
  for (TypeIterBool iterBool = begin; iterBool != end; iterBool++) {
    mounter.PushOrSkip(static_cast<Bool_t>(*iterBool));
  }
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
Bool_t GetIsAcceptableUInt(BranchDescription description) {
  return description.typeName.Contains("UInt") ||
         description.typeName.Contains("unsigned int");
}
Bool_t GetIsAcceptableLong64(BranchDescription description) {
  return description.typeName.Contains("Long64") ||
         description.typeName.Contains("long long") ||
         description.typeName.Contains("Int64");
}
Bool_t GetIsAcceptableLong(BranchDescription description) {
  return (description.typeName.Contains("Long") ||
          description.typeName.Contains("long")) &&
         !GetIsAcceptableLong64(description);
}
Bool_t GetIsAcceptableShort(BranchDescription description) {
  return description.typeName.Contains("Short") ||
         description.typeName.Contains("short") ||
         description.typeName.Contains("Int16");
}
Bool_t GetIsAcceptableChar(BranchDescription description) {
  return description.typeName.Contains("Char") ||
         description.typeName.Contains("char");
}
Bool_t GetIsAcceptableInt(BranchDescription description) {
  return (description.typeName.Contains("Int") ||
          description.typeName.Contains("int")) &&
         !GetIsAcceptableUInt(description) &&
         !GetIsAcceptableLong64(description) &&
         !GetIsAcceptableLong(description);
}
Bool_t GetIsAcceptableBool(BranchDescription description) {
  return description.typeName.Contains("Bool") ||
         description.typeName.Contains("bool");
}
}  // namespace CommonAcceptability

class DescriptionCollectorLongDouble: public VirtualDescriptionCollectorSingle {
 protected:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
};

class BranchMounterScalarLongDouble: public BranchMounterScalarSingle<LongDouble_t> {
 public:
  using BranchMounterScalarSingle<LongDouble_t>::BranchOn;
 protected:
  Bool_t GetIsAccdptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
};

template <typename TypeIterEIn = LongDouble_t*>
class BranchMounterVectorLongDouble: public BranchMounterVectorSingle<LongDouble_t, TypeIterEIn> {
 public:
  using BranchMounterVectorSingle<LongDouble_t, TypeIterEIn>::BranchOn;
 protected:
  Bool_t GetIsAccdptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
};

#endif
