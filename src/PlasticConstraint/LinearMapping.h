#pragma once
#include <sofa/core/Mapping.h>
#include <sofa/core/MultiMapping.h>
#include <sofa/core/Multi2Mapping.h>

namespace sofa::component::mapping::linear
{

namespace crtp
{

template<class TMapping>
class CRTPLinearMapping : public TMapping
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(CRTPLinearMapping, TMapping), TMapping);
    using TMapping::TMapping;

    virtual bool isLinear() const override { return true; }
};

}

template <class TIn, class TOut>
using LinearMapping = crtp::CRTPLinearMapping<core::Mapping<TIn, TOut>>;

template <class TIn, class TOut>
using LinearMultiMapping = crtp::CRTPLinearMapping<core::MultiMapping<TIn, TOut>>;

template <class TIn1, class TIn2, class TOut>
using LinearMulti2Mapping = crtp::CRTPLinearMapping<core::Multi2Mapping<TIn1, TIn2, TOut>>;

using LinearBaseMapping = crtp::CRTPLinearMapping<sofa::core::BaseMapping>;

}
