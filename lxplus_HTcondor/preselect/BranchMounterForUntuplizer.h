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
  kULong64,
  kLong64,
  kULong,
  kLong,
  kUInt,
  kInt,
  kUShort,
  kShort,
  kUChar,
  kChar,
  kBool
};
} // namespace MounterForUntuplizerMeta

class DescriptionCollectorForUntuplizer: public DescriptionCollectorChained<VirtualDescriptionCollectorSingle, true> {
 public:
  using DescriptionCollectorChained<VirtualDescriptionCollectorSingle, true>::vCollectors;
  DescriptionCollectorForUntuplizer(): DescriptionCollectorChained<VirtualDescriptionCollectorSingle, true>({
    // (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(CommonAcceptability::GetIsAcceptableLongDouble),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableDouble    ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableFloat     ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableFloat16   ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableULong64   ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableLong64    ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableULong     ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableLong      ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableUInt      ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableInt       ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableUShort    ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableShort     ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableUChar     ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableChar      ),
    (VirtualDescriptionCollectorSingle *) new DescriptionCollectorFromFun(MounterCommonAcceptability::GetIsAcceptableBool      )
  }) {}
};

class BranchMounterScalarForUntuplizer: public BranchMounterScalarChained<size_t, VirtualBranchMounterScalar<size_t>, true> {
 public:
  BranchMounterScalarForUntuplizer(const DescriptionCollectorForUntuplizer &setting, TreeReader &data,
    const std::function<BranchDescription(BranchDescription)> funDescriptionModified = nullptr)
    : BranchMounterScalarChained<size_t, VirtualBranchMounterScalar<size_t>, true>({
        // new BranchMounterScalarSingle<LongDouble_t>(setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLongDouble],
        //   [&data](BranchDescription description)->LongDouble_t{
        //     return *static_cast<LongDouble_t*>(data.GetPtr(description.name));
        // }),
        new BranchMounterScalarSingle<Double_t    >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kDouble    ],
          [&data](BranchDescription description)->Double_t    {
            return *static_cast<Double_t*>(data.GetPtr(description.name));
        }, funDescriptionModified),
        new BranchMounterScalarSingle<Float_t     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat     ],
          [&data](BranchDescription description)->Float_t     {
            return data.GetFloat(description.name);
        }, funDescriptionModified),
        new BranchMounterScalarSingle<Float16_t   >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat16   ],
          [&data](BranchDescription description)->Float16_t   {
            return *static_cast<Float16_t*>(data.GetPtr(description.name));
        }, funDescriptionModified),
        new BranchMounterScalarSingle<ULong64_t   >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kULong64  ],
          [&data](BranchDescription description)->ULong64_t      {
            return *static_cast<ULong64_t*>(data.GetPtr(description.name, TreeReader::kLong64));
        }, funDescriptionModified),
        new BranchMounterScalarSingle<Long64_t    >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong64    ],
          [&data](BranchDescription description)->Long64_t    {
            return data.GetLong64(description.name);
        }, funDescriptionModified),
        new BranchMounterScalarSingle<ULong_t     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kULong    ],
          [&data](BranchDescription description)->ULong_t      {
            return *static_cast<ULong_t*>(data.GetPtr(description.name));
        }, funDescriptionModified),
        new BranchMounterScalarSingle<Long_t      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong      ],
          [&data](BranchDescription description)->Long_t      {
            return *static_cast<Long_t*>(data.GetPtr(description.name));
        }, funDescriptionModified),
        new BranchMounterScalarSingle<UInt_t      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUInt      ],
          [&data](BranchDescription description)->UInt_t      {
            return *static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kInt));
        }, funDescriptionModified),
        new BranchMounterScalarSingle<Int_t       >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kInt       ],
          [&data](BranchDescription description)->Int_t       {
            return data.GetInt(description.name);
        }, funDescriptionModified),
        new BranchMounterScalarSingle<UShort_t    >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUShort   ],
          [&data](BranchDescription description)->UShort_t      {
            return *static_cast<UShort_t*>(data.GetPtr(description.name, TreeReader::kShort));
        }, funDescriptionModified),
        new BranchMounterScalarSingle<Short_t     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kShort     ],
          [&data](BranchDescription description)->Short_t     {
            return data.GetShort(description.name);
        }, funDescriptionModified),
        new BranchMounterScalarSingle<UChar_t     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUChar      ],
          [&data](BranchDescription description)->UChar_t      {
            return *static_cast<UChar_t*>(data.GetPtr(description.name, TreeReader::kChar));
        }, funDescriptionModified),
        new BranchMounterScalarSingle<Char_t      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kChar      ],
          [&data](BranchDescription description)->Char_t      {
            return data.GetChar(description.name);
        }, funDescriptionModified),
        new BranchMounterScalarSingle<Bool_t      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kBool      ],
          [&data](BranchDescription description)->Bool_t      {
            return data.GetBool(description.name);
        }, funDescriptionModified)
  }) {}
};

class BranchMounterIterableForUntuplizer: public BranchMounterIterableChained<size_t, VirtualBranchMounterIterable<size_t>, true> {
 public:
  BranchMounterIterableForUntuplizer(const DescriptionCollectorForUntuplizer &setting, TreeReader &data,
    const std::function<BranchDescription(BranchDescription)> funDescriptionModified = nullptr)
    : BranchMounterIterableChained<size_t, VirtualBranchMounterIterable<size_t>, true>({
        // new BranchMounterVectorSingle<LongDouble_t, LongDouble_t*>(binding.vvLongDouble, setting->vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLongDouble],
        //   [&data](BranchDescription description)->LongDouble_t*{
        //     return static_cast<LongDouble_t*>(data.GetPtr(description.name));
        // }),
        new BranchMounterVectorSingle<Double_t    , Double_t*    >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kDouble    ],
          [&data](BranchDescription description)->Double_t*    {
            return static_cast<Double_t*>(data.GetPtr(description.name));
        }, funDescriptionModified),
        new BranchMounterVectorSingle<Float_t     , Float_t*     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat     ],
          [&data](BranchDescription description)->Float_t*     {
            return data.GetPtrFloat(description.name);
        }, funDescriptionModified),
        new BranchMounterVectorSingle<Float16_t   , Float16_t*   >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kFloat16   ],
          [&data](BranchDescription description)->Float16_t*   {
            return static_cast<Float16_t*>(data.GetPtr(description.name));
        }, funDescriptionModified),
        new BranchMounterVectorSingle<ULong64_t   , ULong64_t*   >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kULong64   ],
          [&data](BranchDescription description)->ULong64_t*   {
            return static_cast<ULong64_t*>(data.GetPtr(description.name, TreeReader::kArrLong64));
        }, funDescriptionModified),
        new BranchMounterVectorSingle<Long64_t    , Long64_t*    >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong64    ],
          [&data](BranchDescription description)->Long64_t*    {
            return data.GetPtrLong64(description.name);
        }, funDescriptionModified),
        new BranchMounterVectorSingle<ULong_t     , ULong_t*     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kULong     ],
          [&data](BranchDescription description)->ULong_t*     {
            return static_cast<ULong_t*>(data.GetPtr(description.name));
        }, funDescriptionModified),
        new BranchMounterVectorSingle<Long_t      , Long_t*      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kLong      ],
          [&data](BranchDescription description)->Long_t*      {
            return static_cast<Long_t*>(data.GetPtr(description.name));
        }, funDescriptionModified),
        new BranchMounterVectorSingle<UInt_t      , UInt_t*      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUInt      ],
          [&data](BranchDescription description)->UInt_t*      {
            return static_cast<UInt_t*>(data.GetPtr(description.name, TreeReader::kArrInt));
        }, funDescriptionModified),
        new BranchMounterVectorSingle<Int_t       , Int_t*       >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kInt       ],
          [&data](BranchDescription description)->Int_t*       {
            return data.GetPtrInt(description.name);
        }, funDescriptionModified),
        new BranchMounterVectorSingle<UShort_t    , UShort_t*    >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUShort    ],
          [&data](BranchDescription description)->UShort_t*    {
            return static_cast<UShort_t *>(data.GetPtr(description.name, TreeReader::kArrShort));
        }, funDescriptionModified),
        new BranchMounterVectorSingle<Short_t     , Short_t*     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kShort     ],
          [&data](BranchDescription description)->Short_t*     {
            return data.GetPtrShort(description.name);
        }, funDescriptionModified),
        new BranchMounterVectorSingle<UChar_t     , UChar_t*     >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kUChar     ],
          [&data](BranchDescription description)->UChar_t*     {
            return static_cast<UChar_t *>(data.GetPtr(description.name, TreeReader::kArrChar));
        }, funDescriptionModified),
        new BranchMounterVectorSingle<Char_t      , Char_t*      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kChar      ],
          [&data](BranchDescription description)->Char_t*      {
            return data.GetPtrChar(description.name);
        }, funDescriptionModified),
        new BranchMounterVectorSingle<Bool_t      , Bool_t*      >(*setting.vCollectors[(size_t)MounterForUntuplizerMeta::Order::kBool      ],
          [&data](BranchDescription description)->Bool_t*      {
            return data.GetPtrBool(description.name);
        }, funDescriptionModified)
      }) {}
};

#endif // BRANCH_MOUNTER_FOR_UNTUPLIZER_H
