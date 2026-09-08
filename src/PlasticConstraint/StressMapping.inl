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
    Inherit::init();

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

    const InVecCoord& p = this->fromModel->read(core::vec_id::read_access::restPosition)->getValue();
    d_initialPoints.setValue(p);

    helper::ReadAccessor<Data<InVecCoord>> X0 = d_initialPoints;
    elemShapeFun.resize(nTetras);
    B.resize(nTetras);

    for (size_t i = 0; i < nTetras; ++i)
    {
        Mat44 matVert;
        for (Index k = 0; k < 4; k++) {
            Index ix = tetras[i][k];
            matVert(k, 0) = 1.0;
            for (Index l = 1; l < 4; l++)
                matVert(k, l) = X0[ix][l-1];
        }
        const bool canInvert = type::invertMatrix(elemShapeFun[i], matVert);
        assert(canInvert);
        SOFA_UNUSED(canInvert);

        Mat44& shf = elemShapeFun[i];

        for (Index node = 0; node < 4; node++){
            B[i][0][node*3] = shf(1,node);
            B[i][1][node*3+1] = shf(2,node);
            B[i][2][node*3+2] = shf(3,node);

            B[i][3][node*3+2] = 0.5*shf(2,node) ;
            B[i][3][node*3+1] = 0.5*shf(3,node) ;

            B[i][4][node*3+2] = 0.5*shf(1,node) ;
            B[i][4][node*3] = 0.5*shf(3,node) ;

            B[i][5][node*3+1] = 0.5*shf(1,node) ;
            B[i][5][node*3] = 0.5*shf(2,node) ;
        }

    }
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::apply(const core::MechanicalParams* /*mparams*/, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    helper::WriteOnlyAccessor<Data<VecDeriv>> out = dOut;

    const InVecCoord& X = this->fromModel->read(core::vec_id::read_access::position)->getValue();
    helper::ReadAccessor<Data<InVecCoord> > X0 =  d_initialPoints;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    out.resize(tetras.size());

    type::Vec<12,Real> uFlatten;
    InVecCoord U;
    U.resize(X.size());
    for (Index i = 0; i < X0.size(); i++)
        U[i] = X[i] - X0[i];

    type::Vec<6,Real> vStrain;
    Real y = d_youngModulus.getValue()[0];
    Real p = d_poissonRatio.getValue()[0];

    for (size_t i = 0; i < tetras.size(); ++i)
    {

        const auto& it = tetras[i];

        for (Index node = 0; node < 4; node++){
            for (Index k = 0; k < 3; k++) {
                uFlatten[node*3+k] = U[it[node]][k];
            }
        }
        vStrain = B[i]*uFlatten;

        /// stress
        Real traceStrain = vStrain[0] + vStrain[1] + vStrain[2];
        VoigtTensor s;

        for (Index k = 0; k < 3; k++)
            s[k] = (vStrain[k] + (p/(1-2*p))*traceStrain) * (y/(1+p));
        for (Index k = 3; k < 6; k++)
            s[k] = vStrain[k] * (y/(1+p));

        out[i] = s;
    }

}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    helper::WriteOnlyAccessor<Data<VecDeriv>> out = dOut;
    helper::ReadAccessor<Data<InVecDeriv>> in  = dIn;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    out.resize(tetras.size());

    type::Vec<12,Real> vFlatten;
    type::Vec<6,Real> vStrain;
    Real y = d_youngModulus.getValue()[0];
    Real p = d_poissonRatio.getValue()[0];

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        const auto& it = tetras[i];

        for (Index node = 0; node < 4; node++){
            for (Index k = 0; k < 3; k++) {
                vFlatten[node*3+k] = in[it[node]][k];
            }
        }
        vStrain = B[i]*vFlatten;

        /// stress
        Real traceStrain = vStrain[0] + vStrain[1] + vStrain[2];
        VoigtTensor s;
        for (Index k = 0; k < 3; k++)
            s[k] = (vStrain[k] + (p/(1-2*p))*traceStrain) * (y/(1+p));
        for (Index k = 3; k < 6; k++)
            s[k] = vStrain[k] * (y/(1+p));
        out[i] = s;
    }
}

template<class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJT(const core::MechanicalParams* /*mparams*/, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{

    helper::WriteAccessor<Data<InVecDeriv>> out = dOut;
    helper::ReadAccessor<Data<VecDeriv>> in  = dIn;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    Real y = d_youngModulus.getValue()[0];
    Real p = d_poissonRatio.getValue()[0];
    type::Vec<12,Real> outFlatten;

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        const auto& it   = tetras[i];
        const VoigtTensor& eps = in[i];
        Real trEps = eps[0] + eps[1] + eps[2];
        VoigtTensor s; //epsilon point p

        for (Index k = 0; k < 3; k++){
            s[k] = (eps[k] + (p/(1-2*p))*trEps) * (y/(1+p));
        }
        for (Index k = 3; k < 6; k++){
            s[k] = eps[k] * (y/(1+p));
        }

        outFlatten = B[i].transposed()*s;
        for (Index node = 0; node < 4; node++){
            for (Index k = 0; k < 3; k++) {
                out[it[node]][k] += outFlatten[node*3+k];
            }
        }

    }
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJT(const core::ConstraintParams* /*cparams*/, Data<InMatrixDeriv>& dOut, const Data<MatrixDeriv>& dIn)
{
    InMatrixDeriv& out = *dOut.beginEdit();
    const MatrixDeriv& in = dIn.getValue();

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    Real y = d_youngModulus.getValue()[0];
    Real p = d_poissonRatio.getValue()[0];

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != in.end(); ++rowIt)
    {
        auto o = out.writeLine(rowIt.index());
        for (typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin(); colIt != rowIt.end(); ++colIt)
        {
            const Index tetraIdx = colIt.index();
            const VoigtTensor& eps = colIt.val();
            const auto& it  = tetras[tetraIdx];
            type::Vec<12,Real> outFlatten;

            Real trEps = eps[0] + eps[1] + eps[2];
            VoigtTensor s;

            for (Index k = 0; k < 3; k++)
                s[k] = (eps[k] + (p/(1-2*p))*trEps) * (y/(1+p));
            for (Index k = 3; k < 6; k++)
                s[k] = eps[k] * (y/(1+p));
            outFlatten = B[tetraIdx].transposed()*s;

            for (Index node = 0; node < 4; node++){
                for (Index k = 0; k < 3; k++) {
                    InDeriv d;
                    d[k] = outFlatten[node*3+k];
                    o.addCol(it[node], d);
                }
            }
        }
    }
    dOut.endEdit();
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::handleTopologyChange()
{
    if (this->toModel && this->fromModel &&
        this->toModel->getSize() != this->fromModel->getSize())
        this->init();
}

// template <class TIn, class TOut>
// const sofa::linearalgebra::BaseMatrix* StressMapping<TIn, TOut>::getJ()
// {
//     return &J;
// }

template <class TIn, class TOut>
const typename StressMapping<TIn, TOut>::js_type* StressMapping<TIn, TOut>::getJs()
{
    return &Js;
}

} // namespace sofa::component::mapping::linear