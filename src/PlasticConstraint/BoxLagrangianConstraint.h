#pragma once
#include <PlasticConstraint/config.h>

#include <sofa/core/behavior/Constraint.h>
#include <sofa/core/behavior/ConstraintResolution.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/OdeSolver.h>


namespace plasticconstraint::constraint
{

class BoxLagrangianConstraintResolution1Dof : public sofa::core::behavior::ConstraintResolution
{
protected:
    double _invW, _w, _limit, _force ;
    bool _min;


public:

    BoxLagrangianConstraintResolution1Dof(const double &limit, const bool &min, const double &force)
        : sofa::core::behavior::ConstraintResolution(1)
        , _limit(limit)
        , _min(min)
        , _force(force)

    {
    }

    void init(int line, SReal** w, SReal*force) override
    {
        _w = w[line][line];
        _invW = 1.0/_w;
        force[line  ] = 0.0;
    }

    void resolution(int line, SReal** /*w*/, SReal* d, SReal* force, SReal*) override
    {
        const double dfree = d[line] - _w * force[line] + _force * _w;

        if (_min){
            if (dfree < _limit)
                force[line] = (_limit - dfree) * _invW ;
            else
                force[line] = 0.0;
        }
        else
        {
            if (dfree > _limit)
                force[line] = (_limit - dfree) * _invW ;
            else
                force[line] = 0.0;
        }

    }

};

template< class DataTypes >
class BoxLagrangianConstraint : public sofa::core::behavior::Constraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BoxLagrangianConstraint,DataTypes), SOFA_TEMPLATE(sofa::core::behavior::Constraint,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename Coord::value_type Real;
    typedef typename sofa::core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename sofa::core::behavior::Constraint<DataTypes> Inherit;

    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef sofa::core::objectmodel::Data<VecCoord>		DataVecCoord;
    typedef sofa::core::objectmodel::Data<VecDeriv>		DataVecDeriv;
    typedef sofa::core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;
    /// Symmetrical tensor written as a vector following the Voigt notation
    typedef sofa::type::VecNoInit<6,Real> VoigtTensor;


protected:


    sofa::Data<sofa::type::vector<unsigned int>> d_indices; ///< indices of the stop constraints
    sofa::Data<SReal> d_limit; ///< limit value accepted
    sofa::Data<bool> d_min; ///< is the limit a minimum ? (if not it is a maximum)
    sofa::Data<sofa::type::vector<double>> d_force;
    sofa::Data<sofa::type::vector<VoigtTensor>> d_plasticStrain;
    sofa::Data<sofa::type::vector<sofa::type::vector<double>>> d_deviator;
    sofa::Data<int> d_dofPlastic;


    BoxLagrangianConstraint(MechanicalState* object = nullptr);

    virtual ~BoxLagrangianConstraint() {}

    // virtual void init() override;

    virtual sofa::type::vector<std::string> getConstraintIdentifiers() //override final
    {
        sofa::type::vector<std::string> ids = getStopperIdentifiers();
        ids.push_back("Stopper");
        ids.push_back("Unilateral");
        return ids;
    }

    virtual sofa::type::vector<std::string> getStopperIdentifiers(){ return {}; }



public:
    // void init() override;

    void buildConstraintMatrix(const sofa::core::ConstraintParams* cParams, DataMatrixDeriv &c_d, unsigned int &cIndex, const DataVecCoord &x) override;
    void getConstraintViolation(const sofa::core::ConstraintParams* cParams, sofa::linearalgebra::BaseVector *resV, const DataVecCoord &x, const DataVecDeriv &v) override;

    void getConstraintResolution(const sofa::core::ConstraintParams *, std::vector<sofa::core::behavior::ConstraintResolution*>& resTab, unsigned int& offset) override;
    void storeLambda(const sofa::core::ConstraintParams* cParams, sofa::core::MultiVecDerivId res, const sofa::linearalgebra::BaseVector* lambda) override;
};

#if !defined(PLASTICCONSTRAINT_BOXLAGRANGIANCONSTRAINT_CPP)
extern template class SOFA_COMPONENT_CONSTRAINT_LAGRANGIAN_MODEL_API BoxLagrangianConstraint<sofa::defaulttype::Vec1Types>;
extern template class SOFA_COMPONENT_CONSTRAINT_LAGRANGIAN_MODEL_API BoxLagrangianConstraint<sofa::defaulttype::Vec6Types>;

#endif

} //namespace sofa::component::constraint::lagrangian::model
