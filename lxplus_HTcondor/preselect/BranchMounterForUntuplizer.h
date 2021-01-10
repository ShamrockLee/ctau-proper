#ifndef BRANCH_MOUNTER_FOR_UNTUPLIZER_H
#define BRANCH_MOUNTER_FOR_UNTUPLIZER_H
#include <vector>

#include "BranchMounter.h"

#include "untuplizer.h"

namespace MounterForUntuplizerMeta {
enum class Order {
  // kLongDouble,
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
    // (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableLongDouble),
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
  ~DescriptionCollectorForUntuplizer() {
    for (auto collector: vCollectors) {
      delete collector;
    }
    DescriptionCollectorChained::~DescriptionCollectorChained();
  }
};

class BranchMounterScalarForUntuplizer: public BranchMounterScalarChained<VirtualBranchMounterScalar> {
 public:
  BranchMounterScalarForUntuplizer(const DescriptionCollectorForUntuplizer &setting, TreeReader &data)
    : BranchMounterScalarChained<VirtualBranchMounterScalar>({
        // new BranchMounterScalarSingle<LongDouble_t>(setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLongDouble],
        //   [&data](BranchDescription description)->LongDouble_t{
        //     return *static_cast<LongDouble_t*>(data.GetPtr(description.name));
        // }),
        new BranchMounterScalarSingle<Double_t    >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kDouble    ],
          [&data](BranchDescription description)->Double_t    {
            return *static_cast<Double_t*>(data.GetPtr(description.name));
        }),
        new BranchMounterScalarSingle<Float_t     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat     ],
          [&data](BranchDescription description)->Float_t     {
            return data.GetFloat(description.name);
        }),
        new BranchMounterScalarSingle<Float16_t   >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat16   ],
          [&data](BranchDescription description)->Float16_t   {
            return *static_cast<Float16_t*>(data.GetPtr(description.name));
        }),
        new BranchMounterScalarSingle<Long64_t    >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong64    ],
          [&data](BranchDescription description)->Long64_t    {
            return data.GetLong64(description.name);
        }),
        new BranchMounterScalarSingle<Long_t      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong      ],
          [&data](BranchDescription description)->Long_t      {
            return *static_cast<Long_t*>(data.GetPtr(description.name));
        }),
        new BranchMounterScalarSingle<UInt_t      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUInt      ],
          [&data](BranchDescription description)->UInt_t      {
            return *static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kInt));
        }),
        new BranchMounterScalarSingle<Int_t       >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kInt       ],
          [&data](BranchDescription description)->Int_t       {
            return data.GetInt(description.name);
        }),
        new BranchMounterScalarSingle<Short_t     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kShort     ],
          [&data](BranchDescription description)->Short_t     {
            return data.GetShort(description.name);
        }),
        new BranchMounterScalarSingle<Char_t      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kChar      ],
          [&data](BranchDescription description)->Char_t      {
            return data.GetChar(description.name);
        }),
        new BranchMounterScalarSingle<Bool_t      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kBool      ],
          [&data](BranchDescription description)->Bool_t      {
            return data.GetBool(description.name);
        })
  }) {}
  ~BranchMounterScalarForUntuplizer() {
    for (auto mounter: vMounters) {
      delete mounter;
    }
    BranchMounterScalarChained::~BranchMounterScalarChained();
  }
};

class BranchMounterIterableForUntuplizer: public BranchMounterIterableChained<VirtualBranchMounterIterable<size_t>, size_t> {
 public:
  BranchMounterIterableForUntuplizer(const DescriptionCollectorForUntuplizer &setting, TreeReader &data)
    : BranchMounterIterableChained<VirtualBranchMounterIterable<size_t>, size_t>({
        // new BranchMounterVectorSingle<LongDouble_t, LongDouble_t*>(binding.vvLongDouble, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLongDouble],
        //   [&data](BranchDescription description)->LongDouble_t*{
        //     return static_cast<LongDouble_t*>(data.GetPtr(description.name));
        // }),
        new BranchMounterVectorSingle<Double_t, Double_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kDouble],
          [&data](BranchDescription description)->Double_t*{
            return static_cast<Double_t*>(data.GetPtr(description.name));
        }),
        new BranchMounterVectorSingle<Float_t, Float_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat],
          [&data](BranchDescription description)->Float_t*{
            return data.GetPtrFloat(description.name);
        }),
        new BranchMounterVectorSingle<Float16_t, Float16_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat16],
          [&data](BranchDescription description)->Float16_t*{
            return static_cast<Float16_t*>(data.GetPtr(description.name));
        }),
        new BranchMounterVectorSingle<Long64_t, Long64_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong64],
          [&data](BranchDescription description)->Long64_t*{
            return data.GetPtrLong64(description.name);
        }),
        new BranchMounterVectorSingle<Long_t, Long_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong],
          [&data](BranchDescription description)->Long_t*{
            return static_cast<Long_t*>(data.GetPtr(description.name));
        }),
        new BranchMounterVectorSingle<UInt_t, UInt_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUInt],
          [&data](BranchDescription description)->UInt_t*{
            return static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kArrInt));
        }),
        new BranchMounterVectorSingle<Int_t, Int_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kInt],
          [&data](BranchDescription description)->Int_t*{
            return data.GetPtrInt(description.name);
        }),
        new BranchMounterVectorSingle<Short_t, Short_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kShort],
          [&data](BranchDescription description)->Short_t*{
            return data.GetPtrShort(description.name);
        }),
        new BranchMounterVectorSingle<Char_t, Char_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kChar],
          [&data](BranchDescription description)->Char_t*{
            return data.GetPtrChar(description.name);
        }),
        new BranchMounterVectorSingle<Bool_t, Bool_t*>(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kBool],
          [&data](BranchDescription description)->Bool_t*{
            return data.GetPtrBool(description.name);
        })
      }) {}
  ~BranchMounterIterableForUntuplizer() {
    for (auto mounter: vMounters) {
      delete mounter;
    }
    BranchMounterIterableChained::~BranchMounterIterableChained();
  }
};

#endif
