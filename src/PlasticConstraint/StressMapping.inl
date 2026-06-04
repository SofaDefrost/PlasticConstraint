#pragma once
#include <PlasticConstraint/StressMapping.h>
#include <sofa/core/MappingHelper.h>
#include <sofa/core/visual/VisualParams.h>

#include <algorithm>  // std::sort

namespace sofa::component::mapping::linear
{

using sofa::core::objectmodel::ComponentState;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::core::visual::VisualParams;
using sofa::type::vector;

template<class TIn, class TOut>
void StressMapping<TIn, TOut>::init()
{
    m_topology = l_inputTopology.get();
    if (!m_topology)
        this->getContext()->get(m_topology);

    if (!m_topology)
    {
        msg_error() << "No input topology found. StressMapping requires a TetrahedronSetTopology.";
        this->d_componentState.setValue(ComponentState::Invalid);
        return;
    }

    const auto& tetras  = m_topology->getTetrahedra();
    const auto  nTetras = tetras.size();
    const auto  nNodes  = this->fromModel->getSize();

    this->toModel->resize(nTetras);

    Inherit::init();

    // Build J: (nTetras * NOut) x (nNodes * NIn)
    // Each output row i has 4 entries of 0.25 at columns t[0..3]
    J.compressedMatrix.resize(nTetras * NOut, nNodes * NIn);
    J.compressedMatrix.reserve(nTetras * 4 * NOut);

    for (size_t i = 0; i < nTetras; ++i)
    {
        const auto& t = tetras[i];

        sofa::type::fixed_array<sofa::Index, 4> nodes = { t[0], t[1], t[2], t[3] };
        std::sort(nodes.begin(), nodes.end());

        for (unsigned r = 0; r < NOut; ++r)
        {
            const auto row = NOut * i + r;
            J.compressedMatrix.startVec(row);
            for (auto nodeIdx : nodes)
            {
                J.compressedMatrix.insertBack(row, NIn * nodeIdx + r) = (OutReal)0.25;
            }
        }
    }
    J.compressedMatrix.finalize();
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::apply(const core::MechanicalParams* /*mparams*/, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    if (!m_topology)
    {
        msg_error() << "Missing topology in apply().";
        return;
    }

    helper::WriteOnlyAccessor<Data<VecCoord>> out = dOut;
    helper::ReadAccessor<Data<InVecCoord>>    in  = dIn;

    const auto& tetras = m_topology->getTetrahedra();
    out.resize(tetras.size());

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        const auto& t = tetras[i];
        for (unsigned int c = 0; c < Out::spatial_dimensions; ++c)
        {
            out[i][c] = (OutReal)0.25 * (in[t[0]][c] + in[t[1]][c] + in[t[2]][c] + in[t[3]][c]);
        }    
    }
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    helper::WriteOnlyAccessor<Data<VecDeriv>> out = dOut;
    helper::ReadAccessor<Data<InVecDeriv>>    in  = dIn;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    out.resize(tetras.size());

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        const auto& t = tetras[i];
        for (unsigned int c = 0; c < Out::spatial_dimensions; ++c)
        {
            out[i][c] = (OutReal)0.25 * (in[t[0]][c] + in[t[1]][c] + in[t[2]][c] + in[t[3]][c]);
        }    
    }
}

template<class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJT(const core::MechanicalParams* /*mparams*/, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{

}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJT(const core::ConstraintParams* /*cparams*/, Data<InMatrixDeriv>& dOut, const Data<MatrixDeriv>& dIn)
{
 
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::handleTopologyChange()
{
    if (this->toModel && this->fromModel &&
        this->toModel->getSize() != this->fromModel->getSize())
        this->init();
}

template <class TIn, class TOut>
const sofa::linearalgebra::BaseMatrix* StressMapping<TIn, TOut>::getJ()
{
    return &J;
}

template <class TIn, class TOut>
const typename StressMapping<TIn, TOut>::js_type* StressMapping<TIn, TOut>::getJs()
{
    return &Js;
}

} // namespace sofa::component::mapping::linear