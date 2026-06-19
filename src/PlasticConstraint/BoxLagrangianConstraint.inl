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
    , d_index(initData(&d_index, 0, "index", "index of the stop constraint"))
    , d_min(initData(&d_min, -100.0_sreal, "min", "minimum value accepted"))
    , d_max(initData(&d_max, 100.0_sreal, "max", "maximum value accepted"))
{
}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::buildConstraintMatrix(const sofa::core::ConstraintParams* /*cParams*/, DataMatrixDeriv &c_d, unsigned int &cIndex, const DataVecCoord &/*x*/)
{
    auto c = sofa::helper::getWriteAccessor(c_d);

    MatrixDerivRowIterator c_it = c->writeLine(cIndex++);
    c_it.setCol(d_index.getValue(), Coord(1,0,0));
}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* /*cParams*/, sofa::linearalgebra::BaseVector *resV, const DataVecCoord &x, const DataVecDeriv &/*v*/)
{
    const auto constraintIndex = this->d_constraintIndex.getValue();
    resV->set(constraintIndex, x.getValue()[d_index.getValue()][0]);
}

template<class DataTypes>
void BoxLagrangianConstraint<DataTypes>::getConstraintResolution(const sofa::core::ConstraintParams *, std::vector<sofa::core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
    for(int i=0; i<1; i++)
        resTab[offset++] = new BoxLagrangianConstraintResolution1Dof(d_min.getValue(), d_max.getValue());
}

} //namespace sofa::component::constraint::lagrangian::model
