import os

import Sofa
import numpy as np


class validateCylinder(Sofa.PythonScriptController):

    def __init__(self, node):
        self.node = node
        self.bodyID = "cylinderPDMS"
        self.numVolumeElems = 7675
        # self.numVolumeElems = 12724
        # self.numVolumeElems = 23665
        # self.numVolumeElems = 52979
        # self.numVolumeElems = 169436
        # self.numVolumeElems = 288074
        # self.numVolumeElems = 491200
        # self.numVolumeElems = 814336
        self.counter = 0
        self.convergence_threshold = 1e-6
        self.material = "MjedMR"
        self.dt = 0.1
        self.exportVTK = 1
        self.exportStep = 1
        self.exportBegin = 0
        self.exportEnd = 0
        self.steps = 10  # number of gravity loading

        self.simuTag = "mesh_convergence"

        self.color = "1 0 0 1"

        self.simuID = self.bodyID + "_" + str(self.numVolumeElems) + "_" + self.material + "_" + self.simuTag

        self.outDir = "../outData/" + self.simuID
        if not os.path.isdir(self.outDir):
            print("Creating output directory " + self.outDir)
            os.mkdir(self.outDir)

        self.createGraph()

    def createGraph(self):
        self.node.findData('dt').value = self.dt
        self.node.findData('gravity').value = [0, 0, 0]
        self.node.createObject('RequiredPlugin', name='SofaMJEDFEM')
        self.node.createObject('RequiredPlugin', name='SofaExporter')
        self.node.createObject('RequiredPlugin', name='SofaOpenglVisual')
        self.node.createObject('RequiredPlugin', name='SofaPython')
        self.node.createObject('RequiredPlugin', name='SofaSparseSolver')
        self.node.createObject('VisualStyle', displayFlags='showBehaviorModels showForceFields showCollisionModels')

        simuNode1 = self.node.createChild('simulationNode1')
        self.createSimulationNode(simuNode1)

        # self.node/expTip
        expTip = self.node.createChild('expTip')
        expTip.createObject('MechanicalObject', position='0 -0.11625 0.128')
        expTip.createObject('SphereCollisionModel', radius='0.001')
        return 0

    def createSimulationNode(self, simuNode):
        volumeFileName = '../inData/' + self.bodyID + str(self.numVolumeElems) + 'Vol.vtk'
        scale3D = '1 1 1'
        fixedBox = '-0.05 -0.05 -0.002   0.05 0.05 0.0'

        # # simuNode.createObject('EulerImplicitSolver', firstOrder="1")
        self.newton = simuNode.createObject('NewtonStaticSolver', maxIt='30', name='NewtonStatic',
                                            correctionTolerance='1e-6', convergeOnResidual='0',
                                            residualTolerance='1e-6', printLog='1')
        # # simuNode.createObject('StepPCGLinearSolver', name="StepPCG", iterations="10000", tolerance="1e-12", preconditioners="precond", verbose="1", precondOnTimeStep="1")
        # simuNode.createObject('SparsePARDISOSolver', name='precond', exportDataToFolder=self.matDir, symmetric="1", iterativeSolverNumbering="1")
        # # simuNode.createObject('SparseLDLSolver', name='precond', saveMatrixToFile="1")
        # # simuNode.createObject('PCGLinearSolver')
        # # simuNode.createObject('PCGLinearSolver', preconditioners="precond")
        simuNode.createObject('SparseLDLSolver', name='precond', exportDataToDir='inMM')

        # simuNode.createObject("EulerImplicitSolver", name="cg_odesolver", printLog=False)
        # simuNode.createObject('StaticSolver', newton_iterations=100, residual_tolerance_threshold=1e-12)

        # simuNode.createObject("CGLinearSolver", name="linear_solver", tolerance=1e-9, threshold=1e-9)
        simuNode.createObject('MeshVTKLoader', name='loader', filename=volumeFileName, scale3d=scale3D)
        self.dofs = simuNode.createObject('MechanicalObject', src='@loader', name='Volume')
        simuNode.createObject('TetrahedronSetTopologyContainer', src='@loader', name='Container')
        simuNode.createObject('TetrahedronSetTopologyModifier', name='Modifier')
        simuNode.createObject('TetrahedronSetTopologyAlgorithms', name='TopoAlgs')
        simuNode.createObject('TetrahedronSetGeometryAlgorithms', name='GeomAlgs')
        simuNode.createObject('BoxROI', box=fixedBox, drawBoxes='1', name='fixedBox1')
        simuNode.createObject('FixedConstraint', indices='@fixedBox1.indices')

        self.createFEM(simuNode)

        if self.exportVTK:
            fileName = self.outDir + "/" + self.simuID + ".vtk"
            simuNode.createObject('VTKExporter', listening='1', filename=fileName, XMLformat='0',
                                  exportEveryNumberOfSteps=self.exportStep,
                                  exportAtBegin=self.exportBegin, exportAtEnd=self.exportEnd, tetras='1', edges='0')
        return 0

    def createFEM(self, simuNode):

        C01 = 151065.460
        C10 = 101709.668
        K = 1E+07
        objectDensity = '963'
        materialParams = '{} {} {}'.format(C01, C10, K)

        simuNode.createObject('MeshMatrixMass', printMass='0', lumping='1', massDensity=objectDensity, name='mass')
        simuNode.createObject('MJEDTetrahedralForceField', name='FEM', materialName='MooneyRivlin',
                              ParameterSet=materialParams)
        return 0

    def incremental_gravity(self, start, steps):
        self.g = np.linspace(start, -9.81, steps)
        return [0, self.g[self.counter], 0]

    def onBeginAnimationStep(self, deltaTime):
        if self.counter >= self.steps:
            self.position_begin = np.linalg.norm(np.array(self.dofs.position))

    def onEndAnimationStep(self, deltaTime):
        if self.counter < self.steps:
            gravity = self.incremental_gravity(0, self.steps)
            print("gravity = ", gravity)
            self.node.findData('gravity').value = gravity
            self.counter += 1
        else:
            print(self.node.findData('gravity').value)
            self.node.animate = False
            print("Convergence threshold reached")


def createScene(rootNode):
    validateCylinder(rootNode)
    return 0
