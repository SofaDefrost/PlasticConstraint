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
                cubeElasticTranslation = [-15,0,0]
                cubePlasticTranslation = [0,0,0]


                solverNode = rootNode.addChild('solverNode')
                solverNode.addObject('EulerImplicitSolver', name='odesolver',firstOrder="false", rayleighStiffness='0.0', rayleighMass='0.0')
                solverNode.addObject('SparseLDLSolver', name="preconditioner", template="CompressedRowSparseMatrixd")
                solverNode.addObject('GenericConstraintCorrection', linearSolver='@preconditioner',printLog=True, name='ResReso')

                model = solverNode.addChild('model')
                model.addObject('MeshVTKLoader', name='loader', filename=pathMesh+'softCubeSmall.vtk', translation=cubePlasticTranslation)
                model.addObject('MeshTopology',src = '@loader',name="topo")
                model.addObject('MechanicalObject', name='tetras', template='Vec3d', showIndices='false', showIndicesScale='4e-5', rx='0',printLog="0")
                model.addObject('BoxROI', name='ROI1', box='-0.5 -0.5 -0.5 10.5 0.5 10.5', drawBoxes='true')
                model.addObject('BoxROI', name='ROI1force', box='-0.5 9.5 -0.5 10.5 10.5 10.5', drawBoxes='true')

                model.addObject('UniformMass',name="mass", totalMass='0.2', printLog='0')

                model.addObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness = '1e8')
                model.addObject('ConstantForceField', indices='@ROI1force.indices', totalForce = [0,-30000,0])

                ##########################################
                # Plastic Constraint                     #
                ##########################################

                plasticList = list(range(0,24));
                Plastic = model.addChild('Plastic')
                Plastic.addObject('MechanicalObject', template='Vec6', name="Stress")
                Plastic.addObject('BoxLagrangianConstraint', name="boxConstraint", limit=-100, min=True, dof=1, indices=plasticList)
                Plastic.addObject('StressMapping', template='Vec3,Vec6', input='@../tetras', output='@Stress', inputTopology='@../topo', youngModulus=1.0e3, poissonRatio=0.45)

                # j2 = Plastic.addChild('j2')
                # j2.addObject('MechanicalObject', template='Vec1', name="deviator")
                # j2.addObject('DeviatoricMapping',name="deviatorMapping", template='Vec6,Vec1', input='@../Stress', output='@deviator', inputTopology='@../../topo', youngModulus=1.0e3, poissonRatio=0.45)
                # j2.addObject('BoxLagrangianConstraint', name="boxConstraint", min=-10, max=10000, indices=plasticList, deviator="@deviatorMapping.deviateur")

                model.addObject('TetrahedronFEMForceField', template='Vec3d',plasticForcePerElement="@/solverNode/model/Plastic/boxConstraint.plasticStrain",youngModulus=1.0e3, poissonRatio=0.45, method='small')

                visu = model.addChild('visu')
                visu.addObject(  'MeshSTLLoader', name= 'loader', filename=pathMesh+'softCubeSmall.stl')
                visu.addObject('OglModel', src="@loader", color='0.7 0.3 0.5 0.9')
                visu.addObject('BarycentricMapping')


                solverNodeElastic = rootNode.addChild('solverNodeElastic')
                solverNodeElastic.addObject('EulerImplicitSolver', name='odesolver',firstOrder="false", rayleighStiffness='0.0', rayleighMass='0.0')
                solverNodeElastic.addObject('SparseLDLSolver', name="preconditioner", template="CompressedRowSparseMatrixd")
                solverNodeElastic.addObject('GenericConstraintCorrection', linearSolver='@preconditioner',printLog=True, name='ResReso')

                modelElastic = solverNodeElastic.addChild('modelElastic')
                modelElastic.addObject('MeshVTKLoader', name='loader', filename=pathMesh+'softCubeSmall.vtk', translation=cubeElasticTranslation)
                modelElastic.addObject('MeshTopology',src = '@loader')
                modelElastic.addObject('MechanicalObject', name='tetras', template='Vec3d', showIndices='false', showIndicesScale='4e-5', rx='0',printLog="0")
                modelElastic.addObject('BoxROI', name='ROI1', box='-15.5 -0.5 -0.5 -4.5 0.5 10.5', drawBoxes='true')
                modelElastic.addObject('BoxROI', name='ROI1force', box='-15.5 9.5 -0.5 -4.5 10.5 10.5', drawBoxes='true')

                modelElastic.addObject('UniformMass', totalMass='0.2', printLog='0')
                modelElastic.addObject('TetrahedronFEMForceField', template='Vec3d',youngModulus=1.0e3,poissonRatio=0.45, method='small')
                modelElastic.addObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness = '1e8')
                modelElastic.addObject('ConstantForceField', indices='@ROI1force.indices', totalForce = [0,-30000,0])

                visuElastic = modelElastic.addChild('visu')
                visuElastic.addObject(  'MeshSTLLoader', name= 'loader', filename=pathMesh+'softCubeSmall.stl')
                visuElastic.addObject('OglModel', src="@loader", color='0.3 0.7 0.5 0.9', translation=cubeElasticTranslation)
                visuElastic.addObject('BarycentricMapping')

                return rootNode

  