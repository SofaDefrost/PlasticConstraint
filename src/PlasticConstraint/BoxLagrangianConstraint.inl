#pragma once
#include <PlasticConstraint/BoxLagrangianConstraint.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/type/Vec.h>


namespace plasticconstraint::constraint
{

using sofa::core::objectmodel::ComponentState;
using sofa::core::VecCoordId;


template<class DataTypes>
BoxLagrangianConstraint<DataTypes>::BoxLagrangianConstraint(MechanicalState* object)
    : Inherit(object)
    , d_indices(initData(&d_indices, "indices", "indices of the stop constraint"))
    , d_limit(initData(&d_limit, "limit", "limit value accepted"))
    , d_min(initData(&d_min, "min", "is the limit a minumum? If false, it is a maximum"))
    , d_force(initData(&d_force, "force", "plastic strain force"))
    , d_plasticStrain(initData(&d_plasticStrain, "plasticStrain", "plastic strain"))
    , d_dofPlastic(initData(&d_dofPlastic,0, "dof", "which dof of the stress is plastifying (from 0 to 5)"))
    , d_deviator(initData(&d_deviator, "deviator", "Vector of deviator derivative at the current step (Useful only for von Mises plasticity)"))

{
}

template<>
void BoxLagrangianConstraint<sofa::defaulttype::Vec1Types>::buildConstraintMatrix(const sofa::core::ConstraintParams* /*cParams*/, DataMatrixDeriv &c_d, unsigned int &cIndex, const DataVecCoord &/*x*/)
{
    auto c = sofa::helper::getWriteAccessor(c_d);

    if (d_indices.getValue().empty())
    {
        sofa::type::vector<unsigned int> all;

        for (unsigned int i = 0; i < this->getMState()->getSize(); ++i)
            all.push_back(i);
        d_indices.setValue(all);

    }

    const auto& indices = d_indices.getValue();

    for (const auto idx : indices)
    {
        MatrixDerivRowIterator c_it = c->writeLine(cIndex++);
        c_it.setCol(idx, Coord(1));
    }

}
template<>
void BoxLagrangianConstraint<sofa::defaulttype::Vec6Types>::buildConstraintMatrix(const sofa::core::ConstraintParams* /*cParams*/, DataMatrixDeriv &c_d, unsigned int &cIndex, const DataVecCoord &/*x*/)
{
    auto c = sofa::helper::getWriteAccessor(c_d);
    int dof = d_dofPlastic.getValue();
    if (d_indices.getValue().empty())
    {
        sofa::type::vector<unsigned int> all;

        for (unsigned int i = 0; i < this->getMState()->getSize(); ++i)
            all.push_back(i);
        d_indices.setValue(all);

    }

    const auto& indices = d_indices.getValue();

    for (const auto idx : indices)
    {
        MatrixDerivRowIterator c_it = c->writeLine(cIndex++);
        Coord col = Coord(0,0,0,0,0,0);
        col(dof) = 1;
        c_it.setCol(idx, col);
    }

}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* /*cParams*/, sofa::linearalgebra::BaseVector *resV, const DataVecCoord &x, const DataVecDeriv &/*v*/)
{
    const auto constraintIndex = this->d_constraintIndex.getValue();
    const auto& indices = d_indices.getValue();
    const auto& xValue = x.getValue();
    int dof = d_dofPlastic.getValue();

    for (size_t i = 0; i < indices.size(); ++i)
    {
        resV->set(constraintIndex + i, xValue[indices[i]][dof]);
    }
}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::getConstraintResolution(const sofa::core::ConstraintParams *, std::vector<sofa::core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
    const auto& indices = d_indices.getValue();
    if (d_force.getValue().empty()){
        sofa::type::vector<double> tmp;
        sofa::type::vector<VoigtTensor> voigtTmp;
        tmp.resize(indices.size());
        tmp.fill(0.0);
        voigtTmp.resize(indices.size());
        for (size_t i = 0; i < indices.size(); ++i){
            voigtTmp[i].fill(0.0);
        }
        d_force.setValue(tmp);
        d_plasticStrain.setValue(voigtTmp);
    }
    for (size_t i = 0; i < indices.size(); ++i){
        resTab[offset++] = new BoxLagrangianConstraintResolution1Dof(d_limit.getValue(), d_min.getValue(), d_force.getValue()[i]);


    }
}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::storeLambda(const sofa::core::ConstraintParams* cParams,
                                                     sofa::core::MultiVecDerivId res,
                                                     const sofa::linearalgebra::BaseVector* lambda)
{

    const auto& indices = d_indices.getValue();
    sofa::type::vector<double> tmp;
    sofa::type::vector<VoigtTensor> voigtTmp;
    tmp.resize(indices.size());
    voigtTmp.resize(indices.size());

    for (size_t i = 0; i < indices.size(); ++i){
        tmp[i] = d_force.getValue()[i] + lambda->element(this->d_constraintIndex.getValue() + i);
        voigtTmp[i].fill(0);
        voigtTmp[i][d_dofPlastic.getValue()] = tmp[i];
    }
    d_force.setValue(tmp);
    d_plasticStrain.setValue(voigtTmp);
}

} //namespace sofa::component::constraint::lagrangian::model
