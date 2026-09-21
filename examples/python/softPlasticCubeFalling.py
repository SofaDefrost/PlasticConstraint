import Sofa
import Sofa.Core
import os
pathSceneFile = os.path.dirname(os.path.abspath(__file__))
pathMesh = os.path.dirname(os.path.abspath(__file__))+'/meshes/'
# Units: mm, kg, s.     Pressure in kPa = k (kg/(m.s^2)) = k (g/(mm.s^2) =  kg/(mm.s^2)

plugins=['SofaPython3',
         'Sofa.Component.Collision.Detection.Algorithm',
         'Sofa.Component.Visual',
         'Sofa.Component.AnimationLoop',
         'Sofa.Component.Collision.Detection.Intersection',
         'Sofa.Component.Collision.Geometry',
         'Sofa.Component.Collision.Response.Contact',
         'Sofa.Component.Constraint.Lagrangian.Correction',
         'Sofa.Component.Constraint.Lagrangian.Solver',
         'Sofa.Component.IO.Mesh',
         'Sofa.Component.LinearSolver.Direct',
         'Sofa.Component.Mapping.Linear',
         'Sofa.Component.Mass',
         'Sofa.Component.ODESolver.Backward',
         'Sofa.Component.SolidMechanics.FEM.Elastic',
         'Sofa.Component.StateContainer',
         'Sofa.Component.Topology.Container.Constant',
         'Sofa.Component.Topology.Container.Dynamic',
         'Sofa.Component.Visual',
         'Sofa.GL.Component.Rendering3D'
        ]

def createScene(rootNode):
                rootNode.addObject('RequiredPlugin', pluginName=plugins, printLog=False)
                rootNode.addObject("RequiredPlugin", pluginName='PlasticConstraint')

                rootNode.addObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields showInteractionForceFields hideWireframe')
                rootNode.gravity=[0,-9810, 0]
                rootNode.dt = 0.001
                rootNode.addObject('FreeMotionAnimationLoop')
                rootNode.addObject('BlockGaussSeidelConstraintSolver', name='GSSolver', maxIterations='10000', tolerance='1e-15')
                rootNode.addObject('CollisionPipeline', verbose='0')
                rootNode.addObject('BruteForceBroadPhase', name='N2')
                rootNode.addObject('BVHNarrowPhase')
                rootNode.addObject('CollisionResponse', response="FrictionContactConstraint", responseParams="mu=0.0")
                rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance="8.0", contactDistance="0.5", angleCone="0.01")
                cubeElasticTranslation = [-15,10,0]
                cubePlasticTranslation = [0,10,0]

                solverNodeElastic = rootNode.addChild('solverNodeElastic')
                solverNodeElastic.addObject('EulerImplicitSolver', name='odesolver',firstOrder="false", rayleighStiffness='0.0', rayleighMass='0.0')
                solverNodeElastic.addObject('SparseLDLSolver', name="preconditioner", template="CompressedRowSparseMatrixd")
                solverNodeElastic.addObject('GenericConstraintCorrection', linearSolver='@preconditioner',printLog=True, name='ResReso')

                ##########################################
                # FEM Model                              #
                ##########################################
                
                modelElastic = solverNodeElastic.addChild('model')
                modelElastic.addObject('MeshVTKLoader', name='loader', filename=pathMesh+'softCubeSmall.vtk', translation=cubeElasticTranslation)
                modelElastic.addObject('MeshTopology',src = '@loader')
                modelElastic.addObject('MechanicalObject', name='tetras', template='Vec3d', showIndices='false', showIndicesScale='4e-5', rx='0',printLog="0")
                modelElastic.addObject('UniformMass', totalMass='0.65', printLog='0')
                modelElastic.addObject('TetrahedronFEMForceField', template='Vec3d',youngModulus=1.0e3,poissonRatio=0.45, method='small')


                modelCollisElastic = modelElastic.addChild('modelCollisElastic')
                modelCollisElastic.addObject('MeshSTLLoader', name='loader', filename=pathMesh+'softCubeSmall.stl', rotation="0 0 0", translation=cubeElasticTranslation)
                modelCollisElastic.addObject('TriangleSetTopologyContainer', src='@loader', name='container')
                modelCollisElastic.addObject('MechanicalObject', name='collisMO', template='Vec3d')
                modelCollisElastic.addObject('TriangleCollisionModel',group="0")
                modelCollisElastic.addObject('LineCollisionModel',group="0")
                modelCollisElastic.addObject('PointCollisionModel',group="0")
                modelCollisElastic.addObject('BarycentricMapping')

                visuElastic = modelElastic.addChild('visu')
                visuElastic.addObject(  'MeshSTLLoader', name= 'loader', filename=pathMesh+'softCubeSmall.stl')
                visuElastic.addObject('OglModel', src="@loader", color='0.3 0.7 0.5 0.9', translation=cubeElasticTranslation)
                visuElastic.addObject('BarycentricMapping')


                solverNode = rootNode.addChild('solverNode')
                solverNode.addObject('EulerImplicitSolver', name='odesolver',firstOrder="false", rayleighStiffness='0.0', rayleighMass='0.0')
                solverNode.addObject('SparseLDLSolver', name="preconditioner", template="CompressedRowSparseMatrixd")
                solverNode.addObject('GenericConstraintCorrection', linearSolver='@preconditioner',printLog=True, name='ResReso')

                model = solverNode.addChild('model')
                model.addObject('MeshVTKLoader', name='loader', filename=pathMesh+'softCubeSmall.vtk', translation=cubePlasticTranslation)
                model.addObject('MeshTopology', name="topo", src = '@loader')
                model.addObject('MechanicalObject', name='tetras', template='Vec3d', showIndices='false', showIndicesScale='4e-5', rx='0',printLog="0")
                model.addObject('UniformMass', totalMass='0.65', printLog='0')


                ##########################################
                # Plastic Constraint                     #
                ##########################################

                plasticList = list(range(0,24));
                Plastic = model.addChild('Plastic')
                Plastic.addObject('MechanicalObject', template='Vec6', name="Stress")
                Plastic.addObject('BoxLagrangianConstraint', name="boxConstraint", limit=-1, min=True, dof=1, indices=plasticList)

                Plastic.addObject('StressMapping', template='Vec3,Vec6', input='@../tetras', output='@Stress', inputTopology='@../topo', youngModulus=1.0e3, poissonRatio=0.45)
                model.addObject('TetrahedronFEMForceField', template='Vec3d',plasticForcePerElement="@/solverNode/model/Plastic/boxConstraint.plasticStrain",youngModulus=1.0e3, poissonRatio=0.45, method='small')

                # j2 = Plastic.addChild('j2')
                # j2.addObject('MechanicalObject', template='Vec1', name="deviator")
                # j2.addObject('DeviatoricMapping', template='Vec6,Vec1', input='@../Stress', output='@deviator', inputTopology='@../../loader', youngModulus=1.0e3, poissonRatio=0.45)
                # j2.addObject('BoxLagrangianConstraint', name="boxConstraint", min=-100, max=100, indices=plasticList)

                modelCollis = model.addChild('modelCollis')
                modelCollis.addObject('MeshSTLLoader', name='loader', filename=pathMesh+'softCubeSmall.stl', rotation="0 0 0", translation=cubePlasticTranslation)
                #modelCollis.addObject('MeshSTLLoader', name='loader', filename=pathMesh+'sphereFine.stl', rotation="0 0 0", translation=sphereTranslation)
                modelCollis.addObject('TriangleSetTopologyContainer', src='@loader', name='container')
                modelCollis.addObject('MechanicalObject', name='collisMO', template='Vec3d')
                modelCollis.addObject('TriangleCollisionModel',group="0")
                modelCollis.addObject('LineCollisionModel',group="0")
                modelCollis.addObject('PointCollisionModel',group="0")
                modelCollis.addObject('BarycentricMapping')

                rotation=[0,0,0]
                #rotation=[30,0,0]
                planeNode = rootNode.addChild('Plane')
                planeNode.addObject('MeshOBJLoader', name='loader', filename="mesh/floorFlat.obj", triangulate="true", rotation=rotation)
                planeNode.addObject('MeshTopology', src="@loader")
                planeNode.addObject('MechanicalObject', src="@loader", rotation="0 0 0", translation="0 0 0", scale="1")
                planeNode.addObject('TriangleCollisionModel',simulated="0", moving="0",group="1",name='TriPlane')
                planeNode.addObject('LineCollisionModel',simulated="0", moving="0",group="1")
                planeNode.addObject('PointCollisionModel',simulated="0", moving="0",group="1")
                planeNode.addObject('OglModel',name="Visual", src="@loader", color="1 0 0 1", rotation=rotation, translation="0 0 0", scale="1")

                visu = model.addChild('visu')
                visu.addObject(  'MeshSTLLoader', name= 'loader', filename=pathMesh+'softCubeSmall.stl')
                visu.addObject('OglModel', src="@loader", color='0.7 0.3 0.5 0.9', translation=cubePlasticTranslation)
                visu.addObject('BarycentricMapping')


                return rootNode

  