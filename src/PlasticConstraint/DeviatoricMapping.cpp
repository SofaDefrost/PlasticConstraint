#define PLASTICCONSTRAINT_DEVIATORICMAPPING_CPP
#include <PlasticConstraint/DeviatoricMapping.inl>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::component::mapping::linear
{

using namespace sofa::defaulttype;

void registerDeviatoricMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData(
        "PlasticConstraint identity mapping.")
        .add< DeviatoricMapping< Vec6Types,   Vec1Types   > >());
}

template class DeviatoricMapping< Vec6Types,   Vec1Types   >;

} 