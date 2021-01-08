#ifndef BRANCHMOUNTERFORUNTUPLIZER
#define BRANCHMOUNTERFORUNTUPLIZER

#include "untuplizer.h"

#include "BranchMounter.h"

class BranchMounterIterableForUntuplizer: public VirtualBranchMounterIterable<size_t> {
  typedef size_t TypeSize;
 public:
  virtual Bool_t BranchFor(BranchDescription description) {
    return
    mounterLongDouble->BranchFor(description) ||
    mounterDouble->BranchFor(description) ||
    mounterFloat->BranchFor(description) ||
    mounterFloat16->BranchFor(description) ||
    mounterLong64->BranchFor(description) ||
    mounterLong->BranchFor(description) ||
    mounterUInt->BranchFor(description) ||
    mounterInt->BranchFor(description) ||
    mounterShort->BranchFor(description) ||
    mounterChar->BranchFor(description) ||
    mounterBool->BranchFor(description);
  }
  virtual void PrepareForBranching() {
    mounterLongDouble->PrepareForBranching();
    mounterDouble->PrepareForBranching();
    mounterFloat->PrepareForBranching();
    mounterFloat16->PrepareForBranching();
    mounterLong64->PrepareForBranching();
    mounterLong->PrepareForBranching();
    mounterUInt->PrepareForBranching();
    mounterInt->PrepareForBranching();
    mounterShort->PrepareForBranching();
    mounterChar->PrepareForBranching();
    mounterBool->PrepareForBranching();
  }
  virtual void BranchOn(TTree *treeOut) {
    mounterLongDouble->BranchOn(treeOut);
    mounterDouble->BranchOn(treeOut);
    mounterFloat->BranchOn(treeOut);
    mounterFloat16->BranchOn(treeOut);
    mounterLong64->BranchOn(treeOut);
    mounterLong->BranchOn(treeOut);
    mounterUInt->BranchOn(treeOut);
    mounterInt->BranchOn(treeOut);
    mounterShort->BranchOn(treeOut);
    mounterChar->BranchOn(treeOut);
    mounterBool->BranchOn(treeOut);
  }
  virtual void SetTreeReader(TreeReader &data) {
    mounterLongDouble->SetFunIterEIn([&data](BranchDescription description)->LongDouble_t*{
      return static_cast<LongDouble_t*>(data.GetPtr(description.name.Data()));
    });
    mounterDouble->SetFunIterEIn([&data](BranchDescription description)->Double_t*{
      return static_cast<Double_t*>(data.GetPtr(description.name.Data()));
    });
    mounterFloat->SetFunIterEIn([this, &data](BranchDescription description)->Float_t*{
      return data.GetPtrFloat(description.name.Data());
    });
    mounterFloat16->SetFunIterEIn([&data](BranchDescription description)->Float16_t*{
      return static_cast<Float16_t*>(data.GetPtr(description.name.Data()));
    });
    mounterLong64->SetFunIterEIn([&data](BranchDescription description)->Long64_t*{
      return data.GetPtrLong64(description.name.Data());
    });
    mounterLong->SetFunIterEIn([&data](BranchDescription description)->Long_t*{
      return static_cast<Long_t*>(data.GetPtr(description.name.Data()));
    });
    mounterUInt->SetFunIterEIn([&data](BranchDescription description)->UInt_t*{
      return static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kArrInt));
    });
    mounterInt->SetFunIterEIn([&data](BranchDescription description)->Int_t*{
      return data.GetPtrInt(description.name);
    });
    mounterShort->SetFunIterEIn([&data](BranchDescription description)->Short_t*{
      return data.GetPtrShort(description.name);
    });
    mounterChar->SetFunIterEIn([&data](BranchDescription description)->Char_t*{
      return data.GetPtrChar(description.name);
    });
    mounterBool->SetFunIterEIn([&data](BranchDescription description)->Bool_t*{
      return data.GetPtrBool(description.name);
    });
  }
  virtual void PrepareForPushing(size_t n) {
    mounterLongDouble->PrepareForPushing(n);
    mounterDouble->PrepareForPushing(n);
    mounterFloat->PrepareForPushing(n);
    mounterFloat16->PrepareForPushing(n);
    mounterLong64->PrepareForPushing(n);
    mounterLong->PrepareForPushing(n);
    mounterUInt->PrepareForPushing(n);
    mounterInt->PrepareForPushing(n);
    mounterShort->PrepareForPushing(n);
    mounterChar->PrepareForPushing(n);
    mounterBool->PrepareForPushing(n);
  }
  virtual void PushOrSkip(Bool_t isAcceptable) {
    mounterLongDouble->PushOrSkip(isAcceptable);
    mounterDouble->PushOrSkip(isAcceptable);
    mounterFloat->PushOrSkip(isAcceptable);
    mounterFloat16->PushOrSkip(isAcceptable);
    mounterLong64->PushOrSkip(isAcceptable);
    mounterLong->PushOrSkip(isAcceptable);
    mounterUInt->PushOrSkip(isAcceptable);
    mounterInt->PushOrSkip(isAcceptable);
    mounterShort->PushOrSkip(isAcceptable);
    mounterChar->PushOrSkip(isAcceptable);
    mounterBool->PushOrSkip(isAcceptable);
  }
  // struct RefVecs {
  //   std::vector<std::vector<LongDouble_t>> refVVLongDouble;
  //   std::vector<std::vector<Double_t>> refVVDouble;
  //   std::vector<std::vector<Float_t>> refVVFloat;
  //   std::vector<std::vector<Float16_t>> refVVFloat16;
  //   std::vector<std::vector<Long64_t>> refVVLong64;
  //   std::vector<std::vector<Long_t>> refVVLong;
  //   std::vector<std::vector<UInt_t>> refVVUInt;
  //   std::vector<std::vector<Int_t>> refVVInt;
  //   std::vector<std::vector<Short_t>> refVVShort;
  //   std::vector<std::vector<Char_t>> refVVChar;
  //   std::vector<std::vector<Bool_t>> refVVBool;
  //   // RefVecs(refVVLongDouble, refVVDouble, refVVFloat, refVVFloat16, refVVLong64, refVVLong, refVVUInt, refVVInt, refVVShort, refVVChar, refVVBool) {
  //   //   this->refVVLongDoule = refVVLongDouble;
  //   //   this->refVVDouble = refVVDouble;
  //   //   this->refVVFloat = refVVFloat;
  //   //   this->refVVFloat16 = refVVFloat16;
  //   //   this->refVVLong64 = refVVLong64;
  //   //   this->refVVLong = refVVLong;
  //   //   this->refVVUInt = refVVUInt;
  //   //   this->refVVShort = refVVShort;
  //   //   this->refVVChar = refVVChar;
  //   //   this->refVVBool = refVVBool;
  //   // }
  // };

//   BranchMounterIterableForUntuplizer(
//     std::vector<std::vector<LongDouble_t>> &refVVLongDouble,
//     std::vector<std::vector<Double_t>> &refVVDouble,
//     std::vector<std::vector<Float_t>> &refVVFloat,
//     std::vector<std::vector<Float16_t>> &refVVFloat16,
//     std::vector<std::vector<Long64_t>> &refVVLong64,
//     std::vector<std::vector<Long_t>> &refVVLong,
//     std::vector<std::vector<UInt_t>> &refVVUInt,
//     std::vector<std::vector<Int_t>> &refVVInt,
//     std::vector<std::vector<Short_t>> &refVVShort,
//     std::vector<std::vector<Char_t>> &refVVChar,
//     std::vector<std::vector<Bool_t>> &refVVBool) {
//     mounterLongDouble = new BranchMounterLongDouble<LongDouble_t *>(refVVLongDouble);
//     mounterDouble = new BranchMounterDouble<Double_t *>(refVVDouble);
//     mounterFloat = new BranchMounterFloat<Float_t *>(refVVFloat);
//     mounterFloat16 = new BranchMounterFloat16<Float16_t *>(refVVFloat16);
//     mounterLong64 = new BranchMounterLong64<Long64_t *>(refVVLong64);
//     mounterLong = new BranchMounterLong<Long_t *>(refVVLong);
//     mounterUInt = new BranchMounterUInt<UInt_t *>(refVVUInt);
//     mounterInt = new BranchMounterInt<Int_t *>(refVVInt);
//     mounterShort = new BranchMounterShort<Short_t *>(refVVShort);
//     mounterChar = new BranchMounterChar<Char_t *>(refVVChar);
//     mounterBool = new BranchMounterBool<Bool_t *>(refVVBool);
//   }
//  protected:
//   BranchMounterLongDouble<LongDouble_t *> *mounterLongDouble;
//   BranchMounterDouble<Double_t *> *mounterDouble;
//   BranchMounterFloat<Float_t *> *mounterFloat;
//   BranchMounterFloat16<Float16_t *> *mounterFloat16;
//   BranchMounterLong64<Long64_t *> *mounterLong64;
//   BranchMounterLong<Long_t *> *mounterLong;
//   BranchMounterUInt<UInt_t *> *mounterUInt;
//   BranchMounterInt<Int_t *> *mounterInt;
//   BranchMounterShort<Short_t *> *mounterShort;
//   BranchMounterChar<Char_t *> *mounterChar;
//   BranchMounterBool<Bool_t *> *mounterBool;
};

class BranchMounterScalarForUntuplizer: public VirtualBranchMounterScalarBase<size_t> {
 public:
  Bool_t BranchFor(BranchDescription description) {
    return
    mounterLongDoubleScalar->BranchFor(description) ||
    mounterDoubleScalar->BranchFor(description) ||
    mounterFloatScalar->BranchFor(description) ||
    mounterFloat16Scalar->BranchFor(description) ||
    mounterLong64Scalar->BranchFor(description) ||
    mounterLongScalar->BranchFor(description) ||
    mounterUIntScalar->BranchFor(description) ||
    mounterIntScalar->BranchFor(description) ||
    mounterShortScalar->BranchFor(description) ||
    mounterCharScalar->BranchFor(description) ||
    mounterBoolScalar->BranchFor(description);
  }
  void PrepareForBranching() {
    mounterLongDoubleScalar->PrepareForBranching();
    mounterDoubleScalar->PrepareForBranching();
    mounterFloatScalar->PrepareForBranching();
    mounterFloat16Scalar->PrepareForBranching();
    mounterLong64Scalar->PrepareForBranching();
    mounterLongScalar->PrepareForBranching();
    mounterUIntScalar->PrepareForBranching();
    mounterIntScalar->PrepareForBranching();
    mounterShortScalar->PrepareForBranching();
    mounterCharScalar->PrepareForBranching();
    mounterBoolScalar->PrepareForBranching();
  }
  void BranchOn(TTree *treeOut) {
    mounterLongDoubleScalar->BranchOn(treeOut);
    mounterDoubleScalar->BranchOn(treeOut);
    mounterFloatScalar->BranchOn(treeOut);
    mounterFloat16Scalar->BranchOn(treeOut);
    mounterLong64Scalar->BranchOn(treeOut);
    mounterLongScalar->BranchOn(treeOut);
    mounterUIntScalar->BranchOn(treeOut);
    mounterIntScalar->BranchOn(treeOut);
    mounterShortScalar->BranchOn(treeOut);
    mounterCharScalar->BranchOn(treeOut);
    mounterBoolScalar->BranchOn(treeOut);
  }
  void SetTreeReader(TreeReader &data) {
    mounterLongDoubleScalar->SetFunPush([&data](BranchDescription description)->LongDouble_t{
      return *static_cast<LongDouble_t*>(data.GetPtr(description.name.Data()));
    });
    mounterDoubleScalar->SetFunPush([&data](BranchDescription description)->Double_t{
      return *static_cast<Double_t*>(data.GetPtr(description.name.Data()));
    });
    mounterFloatScalar->SetFunPush([&data](BranchDescription description)->Float_t{
      return data.GetFloat(description.name.Data());
    });
    mounterFloat16Scalar->SetFunPush([&data](BranchDescription description)->Float16_t{
      return *static_cast<Float16_t*>(data.GetPtr(description.name.Data()));
    });
    mounterLong64Scalar->SetFunPush([&data](BranchDescription description)->Long64_t{
      return data.GetLong64(description.name.Data());
    });
    mounterLongScalar->SetFunPush([&data](BranchDescription description)->Long_t{
      return *static_cast<Long_t*>(data.GetPtr(description.name.Data()));
    });
    mounterUIntScalar->SetFunPush([&data](BranchDescription description)->UInt_t{
      return *static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kInt));
    });
    mounterIntScalar->SetFunPush([&data](BranchDescription description)->Int_t{
      return data.GetInt(description.name);
    });
    mounterShortScalar->SetFunPush([&data](BranchDescription description)->Short_t{
      return data.GetShort(description.name);
    });
    mounterCharScalar->SetFunPush([&data](BranchDescription description)->Char_t{
      return data.GetChar(description.name);
    });
    mounterBoolScalar->SetFunPush([&data](BranchDescription description)->Bool_t{
      return data.GetBool(description.name);
    });
  }
  void Push() {
    mounterLongDoubleScalar->Push();
    mounterDoubleScalar->Push();
    mounterFloatScalar->Push();
    mounterFloat16Scalar->Push();
    mounterLong64Scalar->Push();
    mounterLongScalar->Push();
    mounterUIntScalar->Push();
    mounterIntScalar->Push();
    mounterShortScalar->Push();
    mounterCharScalar->Push();
    mounterBoolScalar->Push();
  }
  // struct RefVecs {
  //   std::vector<LongDouble_t>& refVLongDouble;
  //   std::vector<Double_t>& refVDouble;
  //   std::vector<Float_t>& refVFloat;
  //   std::vector<Float16_t>& refVFloat16;
  //   std::vector<Long64_t>& refVLong64;
  //   std::vector<Long_t>& refVLong;
  //   std::vector<UInt_t>& refVUInt;
  //   std::vector<Int_t>& refVInt;
  //   std::vector<Short_t>& refVShort;
  //   std::vector<Char_t>& refVChar;
  //   std::vector<Bool_t>& refVBool;
  // };

  // BranchMounterScalarForUntuplizer(
  //   std::vector<LongDouble_t>& refVLongDouble,
  //   std::vector<Double_t>& refVDouble,
  //   std::vector<Float_t>& refVFloat,
  //   std::vector<Float16_t>& refVFloat16,
  //   std::vector<Long64_t>& refVLong64,
  //   std::vector<Long_t>& refVLong,
  //   std::vector<UInt_t>& refVUInt,
  //   std::vector<Int_t>& refVInt,
  //   std::vector<Short_t>& refVShort,
  //   std::vector<Char_t>& refVChar,
  //   std::vector<Bool_t>& refVBool) {
  //   mounterLongDoubleScalar = new BranchMounterLongDoubleScalar(refVLongDouble);
  //   mounterDoubleScalar = new BranchMounterDoubleScalar(refVDouble);
  //   mounterFloatScalar = new BranchMounterFloatScalar(refVFloat);
  //   mounterFloat16Scalar = new BranchMounterFloat16Scalar(refVFloat16);
  //   mounterLong64Scalar = new BranchMounterLong64Scalar(refVLong64);
  //   mounterLongScalar = new BranchMounterLongScalar(refVLong);
  //   mounterUIntScalar = new BranchMounterUIntScalar(refVUInt);
  //   mounterIntScalar = new BranchMounterIntScalar(refVInt);
  //   mounterShortScalar = new BranchMounterShortScalar(refVShort);
  //   mounterCharScalar = new BranchMounterCharScalar(refVChar);
  //   mounterBoolScalar = new BranchMounterBoolScalar(refVBool);
  // }

  BranchMounterScalarForUntuplizer(TreeReader *data) {
    SetTreeReader(data);
  }

  BranchMounterLongDoubleScalar *mounterLongDoubleScalar;
  BranchMounterDoubleScalar *mounterDoubleScalar;
  BranchMounterFloatScalar *mounterFloatScalar;
  BranchMounterFloat16Scalar *mounterFloat16Scalar;
  BranchMounterLong64Scalar *mounterLong64Scalar;
  BranchMounterLongScalar *mounterLongScalar;
  BranchMounterUIntScalar *mounterUIntScalar;
  BranchMounterIntScalar *mounterIntScalar;
  BranchMounterShortScalar *mounterShortScalar;
  BranchMounterCharScalar *mounterCharScalar;
  BranchMounterBoolScalar *mounterBoolScalar;

};
#endif