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

    //using InVecCoord = typename TIn::VecCoord;
    const InVecCoord& p = this->fromModel->read(core::vec_id::read_access::restPosition)->getValue();
    d_initialPoints.setValue(p);

    //copie & adapte de la fonction reinit dans TetrahedronFEMForceField
    helper::ReadAccessor<Data<InVecCoord>> X0 = d_initialPoints;
    elemShapeFun.resize(nTetras);

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
    }
    // elemLambda.resize(nTetras);
    // elemMu.resize(nTetras);

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
    helper::WriteOnlyAccessor<Data<VecDeriv>> out = dOut;
    helper::ReadAccessor<Data<InVecDeriv>>    in  = dIn;

    const InVecCoord& X = this->fromModel->read(core::vec_id::read_access::position)->getValue();
    helper::ReadAccessor<Data<InVecCoord> > X0 =  d_initialPoints;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra(); //tetras correspond au indexed elements du code FEMforcefield
    out.resize(tetras.size());

    InVecCoord U;
    U.resize(X.size());
    for (Index i = 0; i < X0.size(); i++)
        U[i] = X[i] - X0[i];

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        type::Vec<6,Real> vStrain;
        Mat33 gradU;
        const auto& it = tetras[i];
        Mat44& shf = elemShapeFun[i];

        /// compute gradU
        for (Index k = 0; k < 3; k++) {
            for (Index l = 0; l < 3; l++)  {
                gradU(k,l) = 0.0;
                for (Index m = 0; m < 4; m++)
                    gradU(k,l) += shf(l+1,m) * U[it[m]][k];
            }
        }

        Mat33 strain = ((Real)0.5)*(gradU + gradU.transposed() + gradU.transposed()*gradU);

        for (Index i = 0; i < 3; i++)
            vStrain[i] = strain(i,i);
        vStrain[3] = strain(1,2);
        vStrain[4] = strain(0,2);
        vStrain[5] = strain(0,1);

        Real y = d_youngModulus.getValue()[0];
        Real p = d_poissonRatio.getValue()[0];

        Real lambda = (y * p) / ((1 + p) * (1 - 2*p));
        Real mu = y / (2 * (1 + p));

        /// stress
        VoigtTensor s;

        Real traceStrain = 0.0;
        for (Index k = 0; k < 3; k++) {
            traceStrain += vStrain[k];
            s[k] = vStrain[k]*2*mu;
        }

        for (Index k = 3; k < 6; k++)
            s[k] = vStrain[k]*2*mu;

        for (Index k = 0; k < 3; k++)
            s[k] += lambda*traceStrain;
        
        out[i] = s;
    }
}

template <class TIn, class TOut>
void StressMapping<TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    helper::WriteOnlyAccessor<Data<VecDeriv>> out = dOut;
    helper::ReadAccessor<Data<InVecDeriv>>    in  = dIn;

    const InVecCoord& X = this->fromModel->read(core::vec_id::read_access::position)->getValue();
    helper::ReadAccessor<Data<InVecCoord> > X0 =  d_initialPoints;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra(); //tetras correspond au indexed elements du code FEMforcefield
    out.resize(tetras.size());

    InVecCoord U;
    U.resize(X.size());
    for (Index i = 0; i < X0.size(); i++)
        U[i] = X[i] - X0[i];

    for (size_t i = 0; i < tetras.size(); ++i)
    {
        type::Vec<6,Real> vStrain;
        Mat33 gradU;
        const auto& it = tetras[i];
        Mat44& shf = elemShapeFun[i];

        /// compute gradU
        for (Index k = 0; k < 3; k++) {
            for (Index l = 0; l < 3; l++)  {
                gradU(k,l) = 0.0;
                for (Index m = 0; m < 4; m++)
                    gradU(k,l) += shf(l+1,m) * U[it[m]][k];
            }
        }

        Mat33 strain = ((Real)0.5)*(gradU + gradU.transposed() + gradU.transposed()*gradU);

        for (Index i = 0; i < 3; i++)
            vStrain[i] = strain(i,i);
        vStrain[3] = strain(1,2);
        vStrain[4] = strain(0,2);
        vStrain[5] = strain(0,1);

        Real y = d_youngModulus.getValue()[0];
        Real p = d_poissonRatio.getValue()[0];

        Real lambda = (y * p) / ((1 + p) * (1 - 2*p));
        Real mu = y / (2 * (1 + p));

        /// stress
        VoigtTensor s;

        Real traceStrain = 0.0;
        for (Index k = 0; k < 3; k++) {
            traceStrain += vStrain[k];
            s[k] = vStrain[k]*2*mu;
        }

        for (Index k = 3; k < 6; k++)
            s[k] = vStrain[k]*2*mu;

        for (Index k = 0; k < 3; k++)
            s[k] += lambda*traceStrain;
        
        out[i] = s;
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