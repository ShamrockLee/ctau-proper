#ifndef BRANCHMOUNTER_H
#define BRANCHMOUNTER_H

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

template <typename TypeSize = size_t>
class VirtualBranchMounterBase {
 public:
  virtual void BranchOn(TTree *treeOut) = 0;
};
template <typename TypeSize = size_t>
class VirtualBranchMounterScalar : VirtualBranchMounterBase<TypeSize> {
  virtual void Push() = 0;
};
template <typename TypeSize = size_t>
class VirtualBranchMounterIterable : VirtualBranchMounterBase<TypeSize> {
  virtual void PrepareForPushing(const TypeSize n) = 0;
  virtual void PushOrSkip(const Bool_t isAcceptable) = 0;
};

template <class ClassMounter>
class VirtualBranchMounterSetter {
 public:
  virtual Bool_t BranchFor(BranchDescription description) = 0;
  virtual ClassMounter *GetMounter() = 0;
};

template<typename E>
class VirtualBranchMounterScalarSingle: VirtualBranchMounterScalar<size_t> {
 public:
  void SetFunPush(std::function<E(BranchDescription description)> funPush) {
    this->funPush = funPush;
  }
 protected:
  std::function<E(BranchDescription description)> funPush;
};

template <typename E>
class BranchMounterScalar
    : VirtualBranchMounterScalarSingle<E> {
  typedef size_t TypeSize;
 protected:
  const std::vector<BranchDescription> vDescriptions;
  const TypeSize nDescriptions;
  template <size_t capacityDescriptions>
  BranchMounterScalar(E arrE[capacityDescriptions], std::vector<BranchDescription> vDescriptions): arrE(arrE), vDescriptions(vDescriptions), nDescriptions(vDescriptions.size()) {
  }
 public:
  E arrE[];
  void BranchOn(TTree *treeOut) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = vDescriptions[iDescription];
      treeOut->Branch(description.name, &(arrE[iDescription]))->SetTitle(description.title);
    }
  }
  void Push() {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      arrE[iDescription] = funPush(vDescriptions[iDescription]);
    }
  }
};

template<typename E, typename TypeIterEIn = E*>
class VirtualBranchMounterIterableSingle
    : VirtualBranchMounterIterable<size_t> {
 public:
  void SetFunIterEIn(std::function<TypeIterEIn(BranchDescription description)> funIterEIn) {
    this->funIterEIn = funIterEIn;
  }
 protected:
  std::function<TypeIterEIn(BranchDescription description)> funIterEIn;
};

template <typename E, typename TypeIterEIn>
class BranchMounterVector: VirtualBranchMounterIterableSingle<E, TypeIterEIn> {
  typedef size_t TypeSize;
 protected:
  const std::vector<BranchDescription> vDescriptions;
  const size_t nDescriptions;
  std::vector<TypeIterEIn> vIterEIn;
  template<size_t capacityDescriptions>
  BranchMounterVector(std::vector<E> avE[capacityDescriptions], const std::vector<BranchDescription> vDescriptions): avE(avE), vDescriptions(vDescriptions), nDescriptions(0) {
    vIterEIn.reserve(nDescriptions);
  }
 public:
  std::vector<E> avE[capacityDescription];
  void BranchOn(TTree *treeOut) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      const BranchDescription description = vDescriptions[iDescription];
      treeOut->Branch(description.name, &(avE[iDescription]))->SetTitle(description.title);
    }
  }
  void PrepareForPushing(TypeSize n) {
    vIterEIn.clear();
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      avE[iDescription].clear();
      avE[iDescription].reserve(n);
      vIterEIn.emplace_back(funIterEIn(vDescriptions[iDescription]));
    }
  }
  void PushOrSkip(const Bool_t isAcceptable) {
    for (TypeSize iDescription = 0; iDescription < nDescriptions; iDescription++) {
      if (isAcceptable) {
        avE[iDescription].emplace_back(*(vIterEIn[iDescription]));
      }
      vIterEIn[iDescription]++;
    }
  }
};

template<class ClassMounter>
class VirtualBranchMounterSetterSingle: public VirtualBranchMounterSetter<ClassMounter> {
 public:
  virtual Bool_t GetIsAcceptable(BranchDescription description)=0;
  virtual Bool_t BranchFor(BranchDescription description) {
    if (GetIsAcceptable(description)) {
      vDescriptions.emplace_back(description);
      return true;
    }
    return false;
  }
  VirtualBranchMounterSetterSingle() {
    vDescriptions.clear();
  }
 protected:
  std::vector<BranchDescription> vDescriptions;
};

template<typename E>
class BranchMounterSetterScalarSingle: public VirtualBranchMounterSetterSingle<VirtualBranchMounterScalarSingle<E>> {
  typedef VirtualBranchMounterScalarSingle<E> ClassMounter;
  typedef size_t TypeSize;
 public:
  ClassMounter *GetMounter() {
    const TypeSize nDescriptions = vDescriptions.size();
    E arrE[nDescriptions];
    return (ClassMounter *) new BranchMounterScalar<E>(arrE, vDescriptions);
  }
  virtual Bool_t GetIsAcceptable(BranchDescription description)=0;
  virtual Bool_t BranchFor(BranchDescription description) {
    if (GetIsAcceptable(description)) {
      vDescriptions.emplace_back(description);
      return true;
    }
    return false;
  }
  BranchMounterSetterScalarSingle() {
    vDescriptions.clear();
  }
 protected:
  std::vector<BranchDescription> vDescriptions;
};

template<typename E, typename TypeIterEIn = E*>
class BranchMounterSetterVector: public VirtualBranchMounterSetterSingle<VirtualBranchMounterIterableSingle<E, TypeIterEIn>> {
  typedef VirtualBranchMounterIterableSingle<E, TypeIterEIn> ClassMounter;
 public:
  ClassMounter *GetMounter() {
    const TypeSize nDescriptions = vDescriptions.size();
    std::vector<E> avE[];
    return (ClassMounter *) new BranchMounterVector<>(avE, vDescriptions);
  }
  virtual Bool_t GetIsAcceptable(BranchDescription description)=0;
  virtual Bool_t BranchFor(BranchDescription description) {
    if (GetIsAcceptable(description)) {
      vDescriptions.emplace_back(description);
      return true;
    }
    return false;
  }
 protected:
  std::vector<BranchDescription> vDescriptions;
};

// template <typename TypeSize = size_t>
// class VirtualBranchMounterScalarBase
//     : public VirtualBranchMounterBase<TypeSize> {
//   virtual void Push() = 0;
// };
// template <typename E,
//           typename TypeSize = size_t>
// class VirtualBranchMounterScalar
//     : public VirtualBranchMounterScalarBase<TypeSize> {
//   typedef std::function<E(BranchDescription description)> TypeFunPush;
//  public:
//   virtual Bool_t GetIsAcceptable(BranchDescription description) = 0;
//   virtual Bool_t BranchFor(BranchDescription description) {
//     if (GetIsAcceptable(description)) {
//       vDescription.emplace_back(description);
//       return true;
//     }
//     return false;
//   }
//   virtual void PrepareForBranching() {
//     nDescription = vDescription.size();
//     // vE.clear();
//     // vE.resize(nDescription);
//   }
//   virtual void BranchOn(TTree *treeOut) {
//     typename std::vector<BranchDescription>::iterator iterDescription =
//         vDescription.begin();
//     for (TypeSize iDescription = 0; iDescription < nDescription;
//          iDescription++) {
//       treeOut->Branch(iterDescription->name, &(vE[iDescription]))
//           ->SetTitle(iterDescription->title);
//       iterDescription++;
//     }
//   }
//   virtual void SetFunPush(TypeFunPush funPush) { this->funPush = funPush; }
//   virtual void Push() {
//     typename std::vector<BranchDescription>::iterator iterDescription =
//         vDescription.begin();
//     for (TypeSize iDescription = 0; iDescription < nDescription;
//          iDescription++) {
//       vE[iDescription] = funPush(*(iterDescription++));
//     }
//   }

//   // VirtualBranchMounterScalar(std::vector<E> &vE) {
//   //   pvE = &vE;
//   // }

//   std::array<E, capDescription> vE;

//  protected:
//   std::vector<BranchDescription> vDescription;
//   TypeSize nDescription;
//   TypeFunPush funPush;
// };

// template <typename TypeSize = size_t>
// class VirtualBranchMounterIterable : public VirtualBranchMounterBase<TypeSize> {
//  public:
//   virtual void PrepareForPushing(TypeSize n) = 0;
//   virtual void PushOrSkip(Bool_t isAcceptable) = 0;
// };

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

// template <typename E, typename TypeIterEIn = E *>
// class VirtualBranchMounterVector : public VirtualBranchMounterIterable<size_t> {
//   typedef size_t TypeSize;

//  public:
//   virtual Bool_t GetIsAcceptable(BranchDescription description) = 0;
//   virtual Bool_t BranchFor(BranchDescription description) {
//     if (GetIsAcceptable(description)) {
//       vDescription.push_back(description);
//       return true;
//     }
//     return false;
//   }
//   virtual void PrepareForBranching() {
//     nDescriptions = vDescription.size();
//     vIterEIn.reserve(nDescriptions);
//     vvE.clear();
//     vvE.resize(nDescriptions, {});
//   }
//   virtual void BranchOn(TTree *treeOut, std::vector<std::vector<E>> &vvE) {
//     typename std::vector<BranchDescription>::iterator iterDescription =
//         vDescription.begin();
//     for (TypeSize iDescription = 0; iDescription < nDescriptions;
//          iDescription++) {
//       treeOut->Branch(iterDescription->name, &(vvE[iDescription]))
//           ->SetTitle(iterDescription->title);
//     }
//   }
//   virtual void SetFunIterEIn(
//       std::function<TypeIterEIn(BranchDescription description)> funIterEIn) {
//     this->funIterEIn = funIterEIn;
//   }
//   virtual void PrepareForPushing(TypeSize n) {
//     vIterEIn.clear();
//     for (const auto description : vDescription) {
//       vIterEIn.emplace_back(funIterEIn(description));
//     }
//     for (TypeSize iDescription = 0; iDescription < nDescriptions;
//          iDescription++) {
//       vvE[iDescription].clear();
//       vvE[iDescription].reserve(n);
//     }
//   }
//   virtual void PushOrSkip(Bool_t isAcceptable) {
//     for (TypeSize iDescription = 0; iDescription < nDescriptions;
//          iDescription++) {
//       if (isAcceptable)
//         vvE[iDescription].emplace_back(static_cast<E>(*vIterEIn[iDescription]));
//     }
//   }
//   // VirtualBranchMounterVector(std::vector<std::vector<E>> &vvE): vvE(vvE) {
//   //   vDescription.clear();
//   // }
//   VirtualBranchMounterVector(
//       std::function<TypeIterEIn(BranchDescription description)> funIterEIn)
//       : VirtualBranchMounterVector() {
//     this->SetFunIterEIn(funIterEIn);
//   }
//   std::vector<std::vector<E>> &vvE;

//  protected:
//   std::vector<BranchDescription> vDescription;
//   TypeSize nDescriptions;
//   std::function<TypeIterEIn(BranchDescription description)> funIterEIn;
//   std::vector<TypeIterEIn> vIterEIn;
// };

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

class BranchMounterSetterLongDouble
    : public BranchMounterSetterVector<LongDouble_t> {
  typedef LongDouble_t E;
 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
};

class BranchMounterLongDoubleScalar
    : public VirtualBranchMounterScalar<LongDouble_t> {
  typedef LongDouble_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

// #define IMPLEMENT_CLASS_VECTOR_SCALAR(NAME_TYPE_E, NAME_FUNCTION_CA_ACCESSIBLE, NAME_CLASS_VECTOR, NAME_CLASS_SCALAR) \
// template<typename TypeIterEIn=NAME_TYPE_E*> \
// class NAME_CLASS_VECTOR : public VirtualBranchMounterVector<NAME_CLASS_E, TypeIterEIn> { \
//   typedef NAME_CLASS_E E; \
//  public: \
//   Bool_t GetIsAcceptable(BranchDescription description) { \
//     return NAME_FUNCTION_CA_ACCESSIBLE(description); \
//   } \
//   using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector; \
// }; \
//  \
// class NAME_CLASS_SCALAR: public VirtualBranchMounterScalar<NAME_TYPE_E> { \
//  public: \
//   Bool_t GetIsAcceptable(BranchDescription description) { \
//     return NAME_FUNCTION_CA_ACCESSIBLE(description); \
//   } \
// };

// // ROOT preprocessor REFUSE TO UNWRAP #define-macro with parameters!!! QwQ

// IMPLEMENT_CLASS_VECTOR_SCALAR(Double_t,
// CommonAcceptability::GetIsAcceptableDouble, BranchMounterDouble,
// BranchMounterDoubleScalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Float_t,
// CommonAcceptability::GetIsAcceptableFloat, BranchMounterFloat,
// BranchMounterFloatScalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Float16_t,
// CommonAcceptability::GetIsAcceptableFloat16, BranchMounterFloat16,
// BranchMounterFloat16Scalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Float_t,
// CommonAcceptability::GetIsAcceptableFloat, BranchMounterFloat,
// BranchMounterFloatScalar); IMPLEMENT_CLASS_VECTOR_SCALAR(UInt_t,
// CommonAcceptability::GetIsAcceptableUInt, BranchMounterUInt,
// BranchMounterUIntScalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Long64_t,
// CommonAcceptability::GetIsAcceptableLong64, BranchMounterLong64,
// BranchMounterLong64Scalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Long_t,
// CommonAcceptability::GetIsAcceptableLong, BranchMounterLong,
// BranchMounterLongScalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Short_t,
// CommonAcceptability::GetIsAcceptableShort, BranchMounterShort,
// BranchMounterShortScalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Char_t,
// CommonAcceptability::GetIsAcceptableChar, BranchMounterChar,
// BranchMounterCharScalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Int_t,
// CommonAcceptability::GetIsAcceptableInt, BranchMounterInt,
// BranchMounterIntScalar); IMPLEMENT_CLASS_VECTOR_SCALAR(Bool_t,
// CommonAcceptability::GetIsAcceptableBool, BranchMounterBool,
// BranchMounterBoolScalar);

template <typename TypeIterEIn = Double_t *>
class BranchMounterDouble
    : public VirtualBranchMounterVector<Double_t, TypeIterEIn> {
  typedef Double_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterDoubleScalar : public VirtualBranchMounterScalar<Double_t> {
  typedef Double_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = Float_t *>
class BranchMounterFloat
    : public VirtualBranchMounterVector<Float_t, TypeIterEIn> {
  typedef Float_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterFloatScalar : public VirtualBranchMounterScalar<Float_t> {
  typedef Float_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = Float16_t *>
class BranchMounterFloat16
    : public VirtualBranchMounterVector<Float16_t, TypeIterEIn> {
  typedef Float16_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterFloat16Scalar
    : public VirtualBranchMounterScalar<Float16_t> {
  typedef Float16_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = Long64_t *>
class BranchMounterLong64
    : public VirtualBranchMounterVector<Long64_t, TypeIterEIn> {
  typedef Long64_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterLong64Scalar : public VirtualBranchMounterScalar<Long64_t> {
  typedef Long64_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = Long_t *>
class BranchMounterLong
    : public VirtualBranchMounterVector<Long_t, TypeIterEIn> {
  typedef Long_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterLongScalar : public VirtualBranchMounterScalar<Long_t> {
  typedef Long_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = UInt_t *>
class BranchMounterUInt
    : public VirtualBranchMounterVector<UInt_t, TypeIterEIn> {
  typedef UInt_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterUIntScalar : public VirtualBranchMounterScalar<UInt_t> {
  typedef UInt_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = Int_t *>
class BranchMounterInt : public VirtualBranchMounterVector<Int_t, TypeIterEIn> {
  typedef Int_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterIntScalar : public VirtualBranchMounterScalar<Int_t> {
  typedef Int_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = Short_t *>
class BranchMounterShort
    : public VirtualBranchMounterVector<Short_t, TypeIterEIn> {
  typedef Short_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterShortScalar : public VirtualBranchMounterScalar<Short_t> {
  typedef Short_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = Char_t *>
class BranchMounterChar
    : public VirtualBranchMounterVector<Char_t, TypeIterEIn> {
  typedef Char_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterCharScalar : public VirtualBranchMounterScalar<Char_t> {
  typedef Char_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeIterEIn = Bool_t *>
class BranchMounterBool
    : public VirtualBranchMounterVector<Bool_t, TypeIterEIn> {
  typedef Bool_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterVector<E, TypeIterEIn>::VirtualBranchMounterVector;
};

class BranchMounterBoolScalar : public VirtualBranchMounterScalar<Bool_t> {
  typedef Bool_t E;

 public:
  Bool_t GetIsAcceptable(BranchDescription description) {
    return CommonAcceptability::GetIsAcceptableLongDouble(description);
  }
  using VirtualBranchMounterScalar<E>::VirtualBranchMounterScalar;
};

template <typename TypeSize = size_t>
class BranchMounterIterableConcatenated
    : public VirtualBranchMounterIterable<TypeSize> {
  typedef VirtualBranchMounterIterable<TypeSize> TypeMounter;

 public:
  std::vector<TypeMounter *> vMounters;
  BranchMounterIterableConcatenated() { vMounters.clear(); }
  void BranchFor(BranchDescription description) {
    for (TypeMounter *&mounter : vMounters) {
      mounter->BranchFor(description);
    }
  }
  void PrepareForBranching() {
    for (TypeMounter *&mounter : vMounters) {
      mounter->PrepareForBranching();
    }
  }
  void BranchOn(TTree *treeOut) {
    for (TypeMounter *&mounter : vMounters) {
      mounter->BranchOn(treeOut);
    }
  }
  void PrepareForPushing(TypeSize n) {
    for (TypeMounter *&mounter : vMounters) {
      mounter->PrepareForPushing(n);
    }
  }
  void PushOrSkip(Bool_t isAcceptable) {
    for (TypeMounter &mounter : vMounters) {
      mounter.PushOrSkip(isAcceptable);
    }
  }
};
#endif