import os
path = os.path.dirname(os.path.abspath(__file__)) + '/'

def createScene(rootNode):
    rootNode.addObject("RequiredPlugin", pluginName='PlasticConstraint')
    rootNode.addObject("RequiredPlugin", pluginName='SofaPython3')
    rootNode.addObject('RequiredPlugin', pluginName=[
        "Sofa.Component.Collision.Detection.Algorithm",
        "Sofa.Component.Collision.Detection.Intersection",
        "Sofa.Component.Collision.Geometry",
        "Sofa.Component.Collision.Response.Contact",
        "Sofa.Component.Constraint.Projective",
        "Sofa.Component.IO.Mesh",
        "Sofa.Component.LinearSolver.Iterative",
        "Sofa.Component.Mapping.Linear",
        "Sofa.Component.Mass",
        "Sofa.Component.ODESolver.Backward",
        "Sofa.Component.SolidMechanics.FEM.Elastic",
        "Sofa.Component.StateContainer",
        "Sofa.Component.Topology.Container.Constant",
        "Sofa.Component.Topology.Container.Dynamic",
        "Sofa.Component.Visual",
        "Sofa.GL.Component.Rendering3D",
    ])

    rootNode.addObject('VisualStyle', displayFlags='showBehaviorModels showMappings')
    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('CollisionPipeline', verbose=0)
    rootNode.addObject('BruteForceBroadPhase')
    rootNode.addObject('BVHNarrowPhase')
    rootNode.addObject('CollisionResponse', response='PenalityContactForceField')
    rootNode.addObject('MinProximityIntersection', name='Proximity', alarmDistance=0.8, contactDistance=0.5)
    
    tshirt = rootNode.addChild('tshirt')
    tshirt.addObject('EulerImplicitSolver', rayleighStiffness=0.1, rayleighMass=0.1)
    tshirt.addObject('CGLinearSolver', iterations=25, tolerance=1e-5, threshold=1e-5)
    tshirt.addObject('MeshGmshLoader', name='loader', filename=path + 'tshirt_0.msh')
    tshirt.addObject('MeshTopology', src='@loader')
    tshirt.addObject('MechanicalObject', src='@loader', scale=10)
    tshirt.addObject('UniformMass', vertexMass=1)
    tshirt.addObject('FixedProjectiveConstraint', indices=[38, 39, 40, 41, 42, 43, 123, 124, 137])
    tshirt.addObject('TriangleFEMForceField', name='FEM', youngModulus=50000, poissonRatio=0.3, method='large')
    tshirt.addObject('TriangleCollisionModel')

    visu = tshirt.addChild('Visu')
    visu.addObject('OglModel', name='Visual', color='red')
    visu.addObject('StressMapping', input='@..', output='@Visual')

    return rootNode