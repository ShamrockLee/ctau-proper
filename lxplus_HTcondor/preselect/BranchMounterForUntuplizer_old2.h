#ifndef BRANCH_MOUNTER_FOR_UNTUPLIZER_H
#define BRANCH_MOUNTER_FOR_UNTUPLIZER_H
#include <vector>

#include "BranchMounter.h"

#include "untuplizer.h"

namespace MounterForUntuplizerMeta {
size_t NTypes = 11;
enum class Order {
  kLongDouble,
  kDouble,
  kFloat,
  kFloat16,
  kLong64,
  kLong,
  kUInt,
  kInt,
  kShort,
  kChar,
  kBool
};
} // namespace MounterForUntuplizerMeta

class DescriptionCollectorForUntuplizer: public DescriptionCollectorChained<VirtualDescriptionCollectorSingle> {
 public:
  using DescriptionCollectorChained<VirtualDescriptionCollectorSingle>::vCollectors;
  DescriptionCollectorForUntuplizer(): DescriptionCollectorChained<VirtualDescriptionCollectorSingle>({
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableLongDouble),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableDouble    ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableFloat     ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableFloat16   ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableLong64    ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableLong      ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableUInt      ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableInt       ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableShort     ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableChar      ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableBool      ),
  }) {}
};

class BranchMounterScalarForUntuplizer: public BranchMounterScalarChained<VirtualBranchMounterScalar> {
 public:
  typedef struct {
    BranchMounterScalarSingle<LongDouble_t>::Binding vLongDouble;
    BranchMounterScalarSingle<Double_t>::Binding vDouble;
    BranchMounterScalarSingle<Float_t>::Binding vFloat;
    BranchMounterScalarSingle<Float16_t>::Binding vFloat16;
    BranchMounterScalarSingle<Long64_t>::Binding vLong64;
    BranchMounterScalarSingle<Long_t>::Binding vLong;
    BranchMounterScalarSingle<UInt_t>::Binding vUInt;
    BranchMounterScalarSingle<Int_t>::Binding vInt;
    BranchMounterScalarSingle<Short_t>::Binding vShort;
    BranchMounterScalarSingle<Char_t>::Binding vChar;
    BranchMounterScalarSingle<Bool_t>::Binding vBool;
  } Binding;
  BranchMounterScalarSingle<LongDouble_t> *const mounterLongDouble;
  BranchMounterScalarSingle<Double_t    > *const mounterDouble;
  BranchMounterScalarSingle<Float_t     > *const mounterFloat;
  BranchMounterScalarSingle<Float16_t   > *const mounterFloat16;
  BranchMounterScalarSingle<Long64_t    > *const mounterLong64;
  BranchMounterScalarSingle<Long_t      > *const mounterLong;
  BranchMounterScalarSingle<UInt_t      > *const mounterUInt;
  BranchMounterScalarSingle<Int_t       > *const mounterInt;
  BranchMounterScalarSingle<Short_t     > *const mounterShort;
  BranchMounterScalarSingle<Char_t      > *const mounterChar;
  BranchMounterScalarSingle<Bool_t      > *const mounterBool;
  BranchMounterScalarForUntuplizer(Binding &binding, DescriptionCollectorForUntuplizer *setting, TreeReader &data)
      : 
        mounterLongDouble(new BranchMounterScalarSingle<LongDouble_t>(binding.vLongDouble, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLongDouble],
          [&data](BranchDescription description)->LongDouble_t{
            return *static_cast<LongDouble_t*>(data.GetPtr(description.name));
        })),
        mounterDouble(    new BranchMounterScalarSingle<Double_t    >(binding.vDouble,     setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kDouble    ],
          [&data](BranchDescription description)->Double_t    {
            return *static_cast<Double_t*>(data.GetPtr(description.name));
        })),
        mounterFloat(     new BranchMounterScalarSingle<Float_t     >(binding.vFloat,      setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat     ],
          [&data](BranchDescription description)->Float_t     {
            return data.GetFloat(description.name);
        })),
        mounterFloat16(   new BranchMounterScalarSingle<Float16_t   >(binding.vFloat16,    setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat16   ],
          [&data](BranchDescription description)->Float16_t   {
            return *static_cast<Float16_t*>(data.GetPtr(description.name));
        })),
        mounterLong64(    new BranchMounterScalarSingle<Long64_t    >(binding.vLong64,     setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong64    ],
          [&data](BranchDescription description)->Long64_t    {
            return data.GetLong64(description.name);
        })),
        mounterLong(      new BranchMounterScalarSingle<Long_t      >(binding.vLong,       setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong      ],
          [&data](BranchDescription description)->Long_t      {
            return *static_cast<Long_t*>(data.GetPtr(description.name));
        })),
        mounterUInt(      new BranchMounterScalarSingle<UInt_t      >(binding.vUInt,       setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUInt      ],
          [&data](BranchDescription description)->UInt_t      {
            return *static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kInt));
        })),
        mounterInt(       new BranchMounterScalarSingle<Int_t       >(binding.vInt, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kInt              ],
          [&data](BranchDescription description)->Int_t       {
            return data.GetInt(description.name);
        })),
        mounterShort(     new BranchMounterScalarSingle<Short_t     >(binding.vShort, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kShort          ],
          [&data](BranchDescription description)->Short_t     {
            return data.GetShort(description.name);
        })),
        mounterChar(      new BranchMounterScalarSingle<Char_t      >(binding.vChar, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kChar            ],
          [&data](BranchDescription description)->Char_t      {
            return data.GetChar(description.name);
        })),
        mounterBool(      new BranchMounterScalarSingle<Bool_t      >(binding.vBool, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kBool            ],
          [&data](BranchDescription description)->Bool_t      {
            return data.GetBool(description.name);
        })),
        BranchMounterScalarChained<VirtualBranchMounterScalar>({
          mounterLongDouble,
          mounterDouble,
          mounterFloat,
          mounterFloat16,
          mounterLong64,
          mounterLong,
          mounterUInt,
          mounterInt,
          mounterShort,
          mounterChar,
          mounterBool
      }) {}
  void BranchOn(TTree *treeOut, Binding &binding) {
        mounterLongDouble->BranchOn(treeOut, binding.vLongDouble);
        mounterDouble    ->BranchOn(treeOut, binding.vDouble    );
        mounterFloat     ->BranchOn(treeOut, binding.vFloat     );
        mounterFloat16   ->BranchOn(treeOut, binding.vFloat16   );
        mounterLong64    ->BranchOn(treeOut, binding.vLong64    );
        mounterLong      ->BranchOn(treeOut, binding.vLong      );
        mounterUInt      ->BranchOn(treeOut, binding.vUInt      );
        mounterInt       ->BranchOn(treeOut, binding.vInt       );
        mounterShort     ->BranchOn(treeOut, binding.vShort     );
        mounterChar      ->BranchOn(treeOut, binding.vChar      );
        mounterBool      ->BranchOn(treeOut, binding.vBool      );
  }
};

class BranchMounterIterableForUntuplizer: public BranchMounterIterableChained<VirtualBranchMounterIterable<size_t>, size_t> {
 public:
  typedef struct {
    BranchMounterVectorSingle<LongDouble_t, LongDouble_t*>::Binding vvLongDouble;
    BranchMounterVectorSingle<Double_t, Double_t*>::Binding vvDouble;
    BranchMounterVectorSingle<Float_t, Float_t*>::Binding vvFloat;
    BranchMounterVectorSingle<Float16_t, Float16_t*>::Binding vvFloat16;
    BranchMounterVectorSingle<Long64_t, Long64_t*>::Binding vvLong64;
    BranchMounterVectorSingle<Long_t, Long_t*>::Binding vvLong;
    BranchMounterVectorSingle<UInt_t, UInt_t*>::Binding vvUInt;
    BranchMounterVectorSingle<Int_t, Int_t*>::Binding vvInt;
    BranchMounterVectorSingle<Short_t, Short_t*>::Binding vvShort;
    BranchMounterVectorSingle<Char_t, Char_t*>::Binding vvChar;
    BranchMounterVectorSingle<Bool_t, Bool_t*>::Binding vvBool;
  } Binding;
  BranchMounterIterableForUntuplizer(Binding &binding, DescriptionCollectorForUntuplizer *setting, TreeReader &data)
    : BranchMounterIterableChained<VirtualBranchMounterIterable<size_t>, size_t>({
      new BranchMounterVectorSingle<LongDouble_t, LongDouble_t*>(binding.vvLongDouble, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLongDouble],
        [&data](BranchDescription description)->LongDouble_t*{
          return static_cast<LongDouble_t*>(data.GetPtr(description.name));
      }),
      new BranchMounterVectorSingle<Double_t, Double_t*>(binding.vvDouble, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kDouble],
        [&data](BranchDescription description)->Double_t*{
          return static_cast<Double_t*>(data.GetPtr(description.name));
      }),
      new BranchMounterVectorSingle<Float_t, Float_t*>(binding.vvFloat, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat],
        [&data](BranchDescription description)->Float_t*{
          return data.GetPtrFloat(description.name);
      }),
      new BranchMounterVectorSingle<Float16_t, Float16_t*>(binding.vvFloat16, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat16],
        [&data](BranchDescription description)->Float16_t*{
          return static_cast<Float16_t*>(data.GetPtr(description.name));
      }),
      new BranchMounterVectorSingle<Long64_t, Long64_t*>(binding.vvLong64, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong64],
        [&data](BranchDescription description)->Long64_t*{
          return data.GetPtrLong64(description.name);
      }),
      new BranchMounterVectorSingle<Long_t, Long_t*>(binding.vvLong, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong],
        [&data](BranchDescription description)->Long_t*{
          return static_cast<Long_t*>(data.GetPtr(description.name));
      }),
      new BranchMounterVectorSingle<UInt_t, UInt_t*>(binding.vvUInt, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUInt],
        [&data](BranchDescription description)->UInt_t*{
          return static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kArrInt));
      }),
      new BranchMounterVectorSingle<Int_t, Int_t*>(binding.vvInt, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kInt],
        [&data](BranchDescription description)->Int_t*{
          return data.GetPtrInt(description.name);
      }),
      new BranchMounterVectorSingle<Short_t, Short_t*>(binding.vvShort, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kShort],
        [&data](BranchDescription description)->Short_t*{
          return data.GetPtrShort(description.name);
      }),
      new BranchMounterVectorSingle<Char_t, Char_t*>(binding.vvChar, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kChar],
        [&data](BranchDescription description)->Char_t*{
          return data.GetPtrChar(description.name);
      }),
      new BranchMounterVectorSingle<Bool_t, Bool_t*>(binding.vvBool, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kBool],
        [&data](BranchDescription description)->Bool_t*{
          return data.GetPtrBool(description.name);
      })
    }) {}
};

#endif
