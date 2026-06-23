#define PLASTICCONSTRAINT_BOXLAGRANGIANCONSTRAINT_CPP
#include <PlasticConstraint/BoxLagrangianConstraint.inl>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace plasticconstraint::constraint
{

using namespace sofa::defaulttype;

void registerBoxLagrangianConstraint(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Lagrangian-based constraint forcing a 1D DoF to be  inside a given range.")
        //.add< BoxLagrangianConstraint<Vec3Types> >()); //pour test avec vec3
        .add< BoxLagrangianConstraint<Vec6Types> >());

}

//template class PLASTICCONSTRAINT_API BoxLagrangianConstraint<Vec3Types>; //test avec vec3
template class PLASTICCONSTRAINT_API BoxLagrangianConstraint<Vec6Types>;

} //namespace sofa::component::constraint::lagrangian::model
