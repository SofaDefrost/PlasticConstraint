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
        .add< StressMapping< Vec3Types,   Vec3Types   > >()
        .add< StressMapping< Vec2Types,   Vec2Types   > >()
        .add< StressMapping< Vec1Types,   Vec1Types   > >()
        .add< StressMapping< Vec6Types,   Vec3Types   > >()
        .add< StressMapping< Vec6Types,   Vec6Types   > >()
        .add< StressMapping< Rigid3Types, Rigid3Types > >()
        .add< StressMapping< Rigid2Types, Rigid2Types > >()
        .add< StressMapping< Rigid3Types, Vec3Types   > >()
        .add< StressMapping< Rigid2Types, Vec2Types   > >());
}

template class StressMapping< Vec3Types,   Vec3Types   >;
template class StressMapping< Vec2Types,   Vec2Types   >;
template class StressMapping< Vec1Types,   Vec1Types   >;
template class StressMapping< Vec6Types,   Vec3Types   >;
template class StressMapping< Vec6Types,   Vec6Types   >;
template class StressMapping< Rigid3Types, Rigid3Types >;
template class StressMapping< Rigid2Types, Rigid2Types >;
template class StressMapping< Rigid3Types, Vec3Types   >;
template class StressMapping< Rigid2Types, Vec2Types   >;

} 