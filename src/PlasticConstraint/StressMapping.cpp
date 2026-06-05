#define PLASTICCONSTRAINT_STRESSMAPPING_CPP
#include <PlasticConstraint/StressMapping.inl>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::component::mapping::linear
{

using namespace sofa::defaulttype;

void registerStressMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData(
        "PlasticConstraint identity mapping.")
        .add< StressMapping< Vec3Types,   Vec6Types   > >());
}

template class StressMapping< Vec3Types,   Vec6Types   >;

} 