#include <PlasticConstraint/init.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::component::mapping::linear {
    extern void registerStressMapping(sofa::core::ObjectFactory* factory);
}

namespace plasticconstraint
{

void initializePlugin() 
{
    static bool first = true;
    if (first) {
        first = false;
        sofa::core::ObjectFactory* factory = sofa::core::ObjectFactory::getInstance();
        sofa::component::mapping::linear::registerStressMapping(factory);    }
}

}

extern "C" 
{
    PLASTICCONSTRAINT_API void initExternalModule() 
    {
        plasticconstraint::initializePlugin();
    }

    PLASTICCONSTRAINT_API const char* getModuleName() 
    {
        return plasticconstraint::MODULE_NAME;
    }

    PLASTICCONSTRAINT_API const char* getModuleVersion() 
    {
        return plasticconstraint::MODULE_VERSION;
    }

    PLASTICCONSTRAINT_API const char* getModuleLicense() 
    {
        return "LGPL";
    }

    PLASTICCONSTRAINT_API const char* getModuleDescription() 
    {
        return "SOFA plugin for PlasticConstraint";
    }
}
