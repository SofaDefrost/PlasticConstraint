#pragma once
#include <PlasticConstraint/DeviatoricMapping.h>
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
void DeviatoricMapping<TIn, TOut>::init()
{
    Inherit::init();

    m_topology = l_inputTopology.get();
    if (!m_topology)
        this->getContext()->get(m_topology);

    if (!m_topology)
    {
        msg_error() << "No input topology found. DeviatoricMapping requires a TetrahedronSetTopology.";
        this->d_componentState.setValue(ComponentState::Invalid);
        return;
    }
    const auto& tetras  = m_topology->getTetrahedra();
    const auto  nTetras = tetras.size();
    const auto  nNodes  = this->fromModel->getSize();

    this->toModel->resize(nTetras);
    helper::WriteAccessor<Data<type::vector<type::vector<Real> > > > plasticF =  d_deviator;
    plasticF.resize(nTetras);
    for (int i=0; i< nTetras; i++)
        plasticF[i].resize(6);

    const InVecCoord& p = this->fromModel->read(core::vec_id::read_access::restPosition)->getValue();
    d_initialPoints.setValue(p);



}

template <class TIn, class TOut>
void DeviatoricMapping<TIn, TOut>::apply(const core::MechanicalParams* /*mparams*/, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    helper::WriteOnlyAccessor<Data<VecDeriv>> j2 = dOut;
    helper::ReadAccessor<Data<InVecDeriv>> sigma  = dIn;
    auto d_currentSigma = sigma;
    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    j2.resize(tetras.size());
    d_sigma.resize(tetras.size());
    sofa::type::vector<type::vector<Real>> tmp;
    tmp.resize(tetras.size());
    for (size_t elem = 0; elem < tetras.size(); ++elem)
    {
        Real sxx = sigma[elem][0];
        Real syy = sigma[elem][1];
        Real szz = sigma[elem][2];
        Real syz = sigma[elem][3];
        Real sxz = sigma[elem][4];
        Real sxy = sigma[elem][5];
        tmp[elem].resize(6);
        tmp[elem][0] = (2*sxx - syy - szz);
        tmp[elem][1] = (2*syy - sxx - szz);
        tmp[elem][2] = (2*szz - sxx - syy);
        tmp[elem][3] = 6*syz;
        tmp[elem][4] = 6*sxz;
        tmp[elem][5] = 6*sxy;


        for (int i=0; i<6; i++)
            d_sigma[elem][i] = sigma[elem][i];
        j2[elem] = 0.5*(pow(sxx-syy,2) + pow(syy-szz,2) + pow(szz-sxx,2) + 6*(pow(syz,2)+pow(sxz,2)+pow(sxy,2)));
    }
    d_deviator.setValue(tmp);
    // std::cout << "In apply Deviatoric !!!!!!!!!" <<  std::endl;
}

template <class TIn, class TOut>
void DeviatoricMapping<TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    // std::cout << "In applyJ Deviatoric !!!!!!!!!" <<  std::endl;
    const InVecCoord& X = this->fromModel->read(core::vec_id::read_access::position)->getValue();

    helper::WriteOnlyAccessor<Data<VecDeriv>> J2Dot = dOut;
    helper::ReadAccessor<Data<InVecDeriv>> sigmaDot  = dIn;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    J2Dot.resize(tetras.size());

    for (size_t elem = 0; elem < tetras.size(); ++elem)
    {
        Real sDotxx = sigmaDot[elem][0];
        Real sDotyy = sigmaDot[elem][1];
        Real sDotzz = sigmaDot[elem][2];
        Real sDotyz = sigmaDot[elem][3];
        Real sDotxz = sigmaDot[elem][4];
        Real sDotxy = sigmaDot[elem][5];
        Real sxx = X[elem][0];
        Real syy = X[elem][1];
        Real szz = X[elem][2];
        Real syz = X[elem][3];
        Real sxz = X[elem][4];
        Real sxy = X[elem][5];

        J2Dot[elem] = (2*sxx - syy - szz)*sDotxx
                      + (2*syy - sxx - szz)*sDotyy
                      + (2*szz - sxx - syy)*sDotzz
                      + 6*syz*sDotyz
                      + 6*sxz*sDotxz
                      + 6*sxy*sDotxy;
    }
}

template<class TIn, class TOut>
void DeviatoricMapping<TIn, TOut>::applyJT(const core::MechanicalParams* /*mparams*/, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{
    // std::cout << "In applyJT Deviatoric !!!!!!!!!" <<  std::endl;
    const InVecCoord& X = this->fromModel->read(core::vec_id::read_access::position)->getValue();

    helper::WriteAccessor<Data<InVecDeriv>> sigma = dOut;
    helper::ReadAccessor<Data<VecDeriv>> lambda  = dIn;

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();
    sigma.resize(tetras.size());

    for (size_t elem = 0; elem < tetras.size(); ++elem)
    {
        Real sxx = X[elem][0];
        Real syy = X[elem][1];
        Real szz = X[elem][2];
        Real syz = X[elem][3];
        Real sxz = X[elem][4];
        Real sxy = X[elem][5];

        sigma[elem][0] = lambda[elem][0]*(2*sxx - syy - szz);
        sigma[elem][1] = lambda[elem][0]*(2*syy - sxx - szz);
        sigma[elem][2] = lambda[elem][0]*(2*szz - sxx - syy);
        sigma[elem][3] = lambda[elem][0]*6*syz;
        sigma[elem][4] = lambda[elem][0]*6*sxz;
        sigma[elem][5] = lambda[elem][0]*6*sxy;
    }
}

template <class TIn, class TOut>
void DeviatoricMapping<TIn, TOut>::applyJT(const core::ConstraintParams* /*cparams*/, Data<InMatrixDeriv>& dOut, const Data<MatrixDeriv>& dIn)
{
    // std::cout << "In applyJT Constraint Deviatoric !!!!!!!!!" <<  std::endl;
    const InVecCoord& X = this->fromModel->read(core::vec_id::read_access::position)->getValue();

    InMatrixDeriv& out = *dOut.beginEdit();
    const MatrixDeriv& in = dIn.getValue();

    if (!m_topology) return;

    const auto& tetras = m_topology->getTetrahedra();

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != in.end(); ++rowIt)
    {
        auto o = out.writeLine(rowIt.index());
        for (typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin(); colIt != rowIt.end(); ++colIt)
        {
            const Index tetraIdx = colIt.index();
            type::Vec1 lambda = colIt.val();


            InDeriv d;
            Real sxx = X[tetraIdx][0];
            Real syy = X[tetraIdx][1];
            Real szz = X[tetraIdx][2];
            Real syz = X[tetraIdx][3];
            Real sxz = X[tetraIdx][4];
            Real sxy = X[tetraIdx][5];

            d[0] = lambda[0]*(2*sxx - syy - szz);
            d[1] = lambda[0]*(2*syy - sxx - szz);
            d[2] = lambda[0]*(2*szz - sxx - syy);
            d[3] = lambda[0]*6*syz;
            d[4] = lambda[0]*6*sxz;
            d[5] = lambda[0]*6*sxy;

            o.addCol(tetraIdx, d);

        }

    }
    dOut.endEdit();
}

template <class TIn, class TOut>
void DeviatoricMapping<TIn, TOut>::handleTopologyChange()
{
    if (this->toModel && this->fromModel &&
        this->toModel->getSize() != this->fromModel->getSize())
        this->init();
}

// template <class TIn, class TOut>
// const sofa::linearalgebra::BaseMatrix* DeviatoricMapping<TIn, TOut>::getJ()
// {
//     return &J;
// }

template <class TIn, class TOut>
const typename DeviatoricMapping<TIn, TOut>::js_type* DeviatoricMapping<TIn, TOut>::getJs()
{
    return &Js;
}

} // namespace sofa::component::mapping::linear