#ifndef BRANCHMOUNTERFORUNTUPLIZER
#define BRANCHMOUNTERFORUNTUPLIZER

#include "untuplizer.h"

#include "BranchMounter.h"

class BranchMounterForUntuplizer: VirtualBranchMounter<size_t> {
  typedef size_t TypeSize;
 public:
  virtual Bool_t BranchFor(BranchDescription description) {
    return
    mounterLongDouble.BranchFor(description) ||
    mounterDouble.BranchFor(description) ||
    mounterFloat.BranchFor(description) ||
    mounterFloat16.BranchFor(description) ||
    mounterLong64.BranchFor(description) ||
    mounterLong.BranchFor(description) ||
    mounterUInt.BranchFor(description) ||
    mounterInt.BranchFor(description) ||
    mounterShort.BranchFor(description) ||
    mounterChar.BranchFor(description) ||
    mounterBool.BranchFor(description);
  }
  virtual void PrepareForBranching() {
    mounterLongDouble.PrepareForBranching();
    mounterDouble.PrepareForBranching();
    mounterFloat.PrepareForBranching();
    mounterFloat16.PrepareForBranching();
    mounterLong64.PrepareForBranching();
    mounterLong.PrepareForBranching();
    mounterUInt.PrepareForBranching();
    mounterInt.PrepareForBranching();
    mounterShort.PrepareForBranching();
    mounterChar.PrepareForBranching();
    mounterBool.PrepareForBranching();
  }
  virtual void BranchOn(TTree *treeOut) {
    mounterLongDouble.BranchOn(treeOut);
    mounterDouble.BranchOn(treeOut);
    mounterFloat.BranchOn(treeOut);
    mounterFloat16.BranchOn(treeOut);
    mounterLong64.BranchOn(treeOut);
    mounterLong.BranchOn(treeOut);
    mounterUInt.BranchOn(treeOut);
    mounterInt.BranchOn(treeOut);
    mounterShort.BranchOn(treeOut);
    mounterChar.BranchOn(treeOut);
    mounterBool.BranchOn(treeOut);
  }
  virtual void SetTreeReader(TreeReader &data) {
    mounterLongDouble.SetFunIterE([&data](BranchDescription description)->LongDouble_t*{
      return static_cast<LongDouble_t*>(data.GetPtr(description.name.Data()));
    });
    mounterDouble.SetFunIterE([&data](BranchDescription description)->Double_t*{
      return static_cast<Double_t*>(data.GetPtr(description.name.Data()));
    });
    mounterFloat.SetFunIterE([&data](BranchDescription description)->Float_t*{
      return data.GetPtrFloat(description.name.Data());
    });
    mounterFloat16.SetFunIterE([&data](BranchDescription description)->Float16_t*{
      return static_cast<Float16_t*>(data.GetPtr(description.name.Data()));
    });
    mounterLong64.SetFunIterE([&data](BranchDescription description)->Long64_t*{
      return data.GetPtrLong64(description.name.Data());
    });
    mounterLong.SetFunIterE([&data](BranchDescription description)->Long_t*{
      return static_cast<Long_t*>(data.GetPtr(description.name.Data()));
    });
    mounterUInt.SetFunIterE([&data](BranchDescription description)->UInt_t*{
      return static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kArrInt));
    });
    mounterInt.SetFunIterE([&data](BranchDescription description)->Int_t*{
      return data.GetPtrInt(description.name);
    });
    mounterShort.SetFunIterE([&data](BranchDescription description)->Short_t*{
      return data.GetPtrShort(description.name);
    });
    mounterChar.SetFunIterE([&data](BranchDescription description)->Char_t*{
      return data.GetPtrChar(description.name);
    });
    mounterBool.SetFunIterE([&data](BranchDescription description)->Bool_t*{
      return data.GetPtrBool(description.name);
    });
  }
  virtual void PrepareForPushing(size_t n) {
    mounterLongDouble.PrepareForPushing(n);
    mounterDouble.PrepareForPushing(n);
    mounterFloat.PrepareForPushing(n);
    mounterFloat16.PrepareForPushing(n);
    mounterLong64.PrepareForPushing(n);
    mounterLong.PrepareForPushing(n);
    mounterUInt.PrepareForPushing(n);
    mounterInt.PrepareForPushing(n);
    mounterShort.PrepareForPushing(n);
    mounterChar.PrepareForPushing(n);
    mounterBool.PrepareForPushing(n);
  }
  virtual void PushOrSkip(Bool_t isAcceptable) {
    mounterLongDouble.PushOrSkip(isAcceptable);
    mounterDouble.PushOrSkip(isAcceptable);
    mounterFloat.PushOrSkip(isAcceptable);
    mounterFloat16.PushOrSkip(isAcceptable);
    mounterLong64.PushOrSkip(isAcceptable);
    mounterLong.PushOrSkip(isAcceptable);
    mounterUInt.PushOrSkip(isAcceptable);
    mounterInt.PushOrSkip(isAcceptable);
    mounterShort.PushOrSkip(isAcceptable);
    mounterChar.PushOrSkip(isAcceptable);
    mounterBool.PushOrSkip(isAcceptable);
  }
  BranchMounterForUntuplizer(TreeReader &data) {
    SetTreeReader(data);
  }
 protected:
  BranchMounterLongDouble<LongDouble_t *> mounterLongDouble;
  BranchMounterDouble<Double_t *> mounterDouble;
  BranchMounterFloat<Float_t *> mounterFloat;
  BranchMounterFloat16<Float16_t *> mounterFloat16;
  BranchMounterLong64<Long64_t *> mounterLong64;
  BranchMounterLong<Long_t *> mounterLong;
  BranchMounterUInt<UInt_t *> mounterUInt;
  BranchMounterInt<Int_t *> mounterInt;
  BranchMounterShort<Short_t *> mounterShort;
  BranchMounterChar<Char_t *> mounterChar;
  BranchMounterBool<Bool_t *> mounterBool;
};
#endif