#> This is an example program which solves a weakly coupled 
#> Finite Elasticity-ALE NavierStokes equation using OpenCMISS
#> calls.
#>
#> By Chris Bradley, Andreas Hessenthaler, Soroush Safaei, Zohreh Ekhlasi
#>

#================================================================================================================================
#  Symbol Definitions
#================================================================================================================================

FLUID = 1
SOLID = 2
FSI = 3

LINEAR_LAGRANGE = 1
QUADRATIC_LAGRANGE = 2
CUBIC_LAGRANGE = 3
CUBIC_HERMITE = 4
LINEAR_SIMPLEX = 5
QUADRATIC_SIMPLEX = 6
CUBIC_SIMPLEX = 7

NOTHING = 1
VELOCITY = 2
PRESSURE = 3
REFPRESSURE = 4

#================================================================================================================================
#  User changeable example parameters
#================================================================================================================================

problemType = FSI
#problemType = FLUID

width = 3.0
height = 1.5

numberOfSolidXElements = 1
numberOfSolidYElements = 2
numberOfFluidX1Elements = 2
numberOfFluidX2Elements = 4
numberOfFluidYElements = 2

#uInterpolation = QUADRATIC_LAGRANGE
#pInterpolation = LINEAR_LAGRANGE
uInterpolation = QUADRATIC_SIMPLEX
pInterpolation = LINEAR_SIMPLEX

#RBS = False
RBS = True

outputFrequency = 1 # Result output frequency

setupOutput = True
progressDiagnostics = True
debugLevel = 3

# Temporal information
startTime = 0.0
stopTime  = 5.01
timeStep  = 0.1

# Inlet velocity parameters
A = 0.5
B = 2.0
C = -0.5

# Material properties
fluidDynamicViscosity = 0.05  # kg / (m s)
fluidDensity  = 100           # kg m^-3
solidDensity  = 300           # kg m^-3
mooneyRivlin1 = 2.0           # N / m^2
mooneyRivlin2 = 4.0
# Moving mesh
movingMeshKParameter   = 1.0       #default

solidPRef = 0.0
fluidPRef = 0.0

solidPInit = -mooneyRivlin1
fluidPInit = fluidPRef

# Set solver parameters
fsiDynamicSolverTheta    = [1.0]
nonlinearMaximumIterations      = 100000000 #default: 100000
nonlinearRelativeTolerance      = 1.0E-4    #default: 1.0E-05
nonlinearAbsoluteTolerance      = 1.0E-4    #default: 1.0E-10
nonlinearMaxFunctionEvaluations = 100000
nonlinearLinesearchAlpha        = 1.0
linearMaximumIterations      = 100000000 #default: 100000
linearRelativeTolerance      = 1.0E-4    #default: 1.0E-05
linearAbsoluteTolerance      = 1.0E-4    #default: 1.0E-10
linearDivergenceTolerance    = 1.0E5     #default: 1.0E5
linearRestartValue           = 30        #default: 30

#================================================================================================================================
#  Should not need to change anything below here.
#================================================================================================================================

numberOfDimensions = 2
numberOfInterfaceDimensions = 1

if (uInterpolation == LINEAR_LAGRANGE):
    numberOfNodesXi = 2
    numberOfGaussXi = 2
    simplex = False
elif (uInterpolation == QUADRATIC_LAGRANGE):
    numberOfNodesXi = 3
    numberOfGaussXi = 3
    simplex = False
elif (uInterpolation == CUBIC_LAGRANGE):
    numberOfNodesXi = 4
    numberOfGaussXi = 3
    simplex = False
elif (uInterpolation == CUBIC_HERMITE):
    numberOfNodesXi = 2
    numberOfGaussXi = 3
    simplex = False
elif (uInterpolation == LINEAR_SIMPLEX):
    numberOfNodesXi = 2
    gaussOrder = 2
    simplex = True
    simplexOrder = 1
elif (uInterpolation == QUADRATIC_SIMPLEX):
    numberOfNodesXi = 3
    gaussOrder = 4
    simplex = True
    simplexOrder = 2
elif (uInterpolation == CUBIC_SIMPLEX):
    numberOfNodesXi = 4
    gaussOrder = 5
    simplex = True
    simplexOrder = 3
else:
    print('Invalid u interpolation')
    exit()
    
if (simplex):
    numberOfSubElements = 2
    if (not ((pInterpolation == LINEAR_SIMPLEX) or (pInterpolation == QUADRATIC_SIMPLEX) or (pInterpolation == CUBIC_SIMPLEX))):
        print('If the u interpolation is of simplex type the p interpolation must be of simplex type.')
        exit()
else:
    numberOfSubElements = 1

xElementSize = width/(numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements)
yElementSize = height/(numberOfSolidYElements+numberOfFluidYElements)
solidXSize = numberOfSolidXElements*xElementSize
solidYSize = numberOfSolidYElements*yElementSize
fluidX1Size = numberOfFluidX1Elements*xElementSize
fluidX2Size = numberOfFluidX2Elements*xElementSize
fluidYSize = numberOfFluidYElements*yElementSize

numberOfSolidXNodes = numberOfSolidXElements*(numberOfNodesXi-1)+1
numberOfSolidYNodes = numberOfSolidYElements*(numberOfNodesXi-1)+1
numberOfSolidNodes = (numberOfSolidXElements*(numberOfNodesXi-1)+1)*(numberOfSolidYElements*(numberOfNodesXi-1)+1)
numberOfFluidX1Nodes = numberOfFluidX1Elements*(numberOfNodesXi-1)+1
numberOfFluidX2Nodes = numberOfFluidX2Elements*(numberOfNodesXi-1)+1
numberOfFluidXNodes1 = numberOfFluidX1Nodes+numberOfFluidX2Nodes+numberOfSolidXNodes-2
numberOfFluidXNodes2 = numberOfFluidX1Nodes+numberOfFluidX2Nodes
numberOfFluidYNodes = numberOfFluidYElements*(numberOfNodesXi-1)+1
numberOfFluidNodes = (numberOfFluidX1Elements*(numberOfNodesXi-1)+1)*(numberOfSolidYElements*(numberOfNodesXi-1))+ \
                    (numberOfFluidX2Elements*(numberOfNodesXi-1)+1)*(numberOfSolidYElements*(numberOfNodesXi-1))+ \
                    ((numberOfFluidX1Elements+numberOfFluidX2Elements+numberOfSolidXElements)*(numberOfNodesXi-1)+1)* \
                    (numberOfFluidYElements*(numberOfNodesXi-1)+1)
numberOfInterfaceNodes = (numberOfSolidXElements*(numberOfNodesXi-1)+1)+2*numberOfSolidYElements*(numberOfNodesXi-1)
numberOfSolidElements = numberOfSolidXElements*numberOfSolidYElements*numberOfSubElements
numberOfFluidXElements1 = (numberOfFluidX1Elements+numberOfFluidX2Elements+numberOfSolidXElements)*numberOfSubElements
numberOfFluidXElements2 = (numberOfFluidX1Elements+numberOfFluidX2Elements)*numberOfSubElements
numberOfFluidElements = ((numberOfFluidX1Elements+numberOfFluidX2Elements+numberOfSolidXElements)*numberOfSubElements* \
                        (numberOfSolidYElements+numberOfFluidYElements) - numberOfSolidElements)
numberOfInterfaceElements = numberOfSolidXElements+2*numberOfSolidYElements

numberOfLocalInterfaceNodes = numberOfNodesXi
localNodeIdx0=0
localNodeIdx1=numberOfNodesXi-1
if (simplex):
    if (uInterpolation == LINEAR_SIMPLEX):
        numberOfLocalNodes = 3
    elif (uInterpolation == QUADRATIC_SIMPLEX):
        numberOfLocalNodes = 6
    elif (uInterpolation == CUBIC_SIMPLEX):
        numberOfLocalNodes = 10
    localNodeIdx00 = 0
    localNodeIdx100 = 0
    localNodeIdx010 = 1
    localNodeIdx001 = 2
else:
    numberOfLocalNodes = numberOfNodesXi*numberOfNodesXi
    localNodeIdx00 = 0
    localNodeIdx10 = numberOfNodesXi-1
    localNodeIdx01 = numberOfNodesXi*(numberOfNodesXi-1)
    localNodeIdx11 = numberOfNodesXi*numberOfNodesXi-1
    
solidCoordinateSystemUserNumber     = 1
fluidCoordinateSystemUserNumber     = 2
interfaceCoordinateSystemUserNumber = 3

solidRegionUserNumber = 1
fluidRegionUserNumber = 2
interfaceUserNumber   = 3

uBasisUserNumber = 1
pBasisUserNumber = 2
interfaceBasisUserNumber = 3

solidMeshUserNumber     = 1
fluidMeshUserNumber     = 2
interfaceMeshUserNumber = 3
movingMeshUserNumber    = 4
          
solidDecompositionUserNumber     = 1
fluidDecompositionUserNumber     = 2
interfaceDecompositionUserNumber = 3
          
solidGeometricFieldUserNumber     = 11
solidFibreFieldUserNumber     = 12
solidEquationsSetFieldUserNumber = 13
solidDependentFieldUserNumber = 14
solidMaterialsFieldUserNumber = 15
solidSourceFieldUserNumber = 16

fluidGeometricFieldUserNumber     = 21
fluidEquationsSetFieldUserNumber = 22
fluidDependentFieldUserNumber = 23
fluidMaterialsFieldUserNumber = 24
fluidIndependentFieldUserNumber = 25
bcCellMLModelsFieldUserNumber = 26
bcCellMLStateFieldUserNumber = 27
bcCellMLParametersFieldUserNumber = 28
bcCellMLIntermediateFieldUserNumber = 29

movingMeshEquationsSetFieldUserNumber = 31
movingMeshDependentFieldUserNumber    = 32
movingMeshMaterialsFieldUserNumber    = 33
movingMeshIndependentFieldUserNumber  = 34

interfaceGeometricFieldUserNumber = 41
interfaceLagrangeFieldUserNumber  = 42

solidEquationsSetUserNumber  = 1
fluidEquationsSetUserNumber  = 2
movingMeshEquationsSetUserNumber = 3

bcCellMLUserNumber = 1

interfaceConditionUserNumber = 1
          
fsiProblemUserNumber = 1

#================================================================================================================================
#  Define functions
#================================================================================================================================

def GetElementNodes2D(elementNumber,subElementNumber,localNodes2D,numberOfXNodes1,numberOfXNodes2):
    if (uInterpolation == LINEAR_LAGRANGE or uInterpolation == CUBIC_HERMITE):
        localNodes2D[localNodeIdx10] = localNodes2D[localNodeIdx00]+(numberOfNodesXi-1)
        localNodes2D[localNodeIdx01] = localNodes2D[localNodeIdx00]+numberOfXNodes2
        localNodes2D[localNodeIdx11] = localNodes2D[localNodeIdx01]+(numberOfNodesXi-1)
        uNodes2D = [localNodes2D[localNodeIdx00],localNodes2D[localNodeIdx10],localNodes2D[localNodeIdx01],localNodes2D[localNodeIdx11]]
    elif (uInterpolation == QUADRATIC_LAGRANGE):
        localNodes2D[localNodeIdx10] = localNodes2D[localNodeIdx00]+(numberOfNodesXi-1)
        localNodes2D[localNodeIdx01] = localNodes2D[localNodeIdx00]+numberOfXNodes1+numberOfXNodes2
        localNodes2D[localNodeIdx11] = localNodes2D[localNodeIdx01]+(numberOfNodesXi-1)
        localNodes2D[1] = localNodes2D[localNodeIdx00] + 1
        localNodes2D[3] = localNodes2D[localNodeIdx00] + numberOfXNodes1
        localNodes2D[4] = localNodes2D[3] + 1
        localNodes2D[5] = localNodes2D[4] + 1
        localNodes2D[7] = localNodes2D[localNodeIdx01] + 1
        uNodes2D = localNodes2D
    elif (uInterpolation == CUBIC_LAGRANGE):
        localNodes2D[localNodeIdx10] = localNodes2D[localNodeIdx00]+(numberOfNodesXi-1)
        localNodes2D[localNodeIdx01] = localNodes2D[localNodeIdx00]+2*numberOfXNodes1+numberOfXNodes2
        localNodes2D[localNodeIdx11] = localNodes2D[localNodeIdx01]+(numberOfNodesXi-1)
        localNodes2D[1] = localNodes2D[localNodeIdx00] + 1
        localNodes2D[2] = localNodes2D[1] + 1
        localNodes2D[4] = localNodes2D[localNodeIdx00] + numberOfXNodes1
        localNodes2D[5] = localNodes2D[4] + 1
        localNodes2D[6] = localNodes2D[5] + 1
        localNodes2D[7] = localNodes2D[6] + 1
        localNodes2D[8] = localNodes2D[4] + numberOfXNodes1
        localNodes2D[9] = localNodes2D[8] + 1
        localNodes2D[10] = localNodes2D[9] + 1
        localNodes2D[11] = localNodes2D[10] + 1
        localNodes2D[13] = localNodes2D[localNodeIdx01] + 1
        localNodes2D[14] = localNodes2D[13] + 1
        uNodes2D = localNodes2D
    elif (uInterpolation == LINEAR_SIMPLEX):
        if (subElementNumber == 1):
            localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx100] + numberOfXNodes2 + 1
            localNodes2D[localNodeIdx001] = localNodes2D[localNodeIdx100] + numberOfXNodes2
        else:
            localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx100] + 1
            localNodes2D[localNodeIdx001] = localNodes2D[localNodeIdx100] + numberOfXNodes2 + 1
        uNodes2D = [localNodes2D[localNodeIdx100],localNodes2D[localNodeIdx010],localNodes2D[localNodeIdx001]]
    elif (uInterpolation == QUADRATIC_SIMPLEX):
        if (subElementNumber == 1):
            localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx100] + numberOfXNodes1 + numberOfXNodes2 + 2
            localNodes2D[localNodeIdx001] = localNodes2D[localNodeIdx100] + numberOfXNodes1 + numberOfXNodes2
            localNodes2D[3] = localNodes2D[localNodeIdx100] + numberOfXNodes1 + 1
            localNodes2D[4] = localNodes2D[localNodeIdx100] + numberOfXNodes1 + numberOfXNodes2 + 1
            localNodes2D[5] = localNodes2D[localNodeIdx100] + numberOfXNodes1
        else:
	    localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx100] + 2
            localNodes2D[localNodeIdx001] = localNodes2D[localNodeIdx100] + numberOfXNodes1 + numberOfXNodes2 + 2
            localNodes2D[3] = localNodes2D[localNodeIdx100] + 1
            localNodes2D[4] = localNodes2D[localNodeIdx100] + numberOfXNodes1 + 2
            localNodes2D[5] = localNodes2D[localNodeIdx00] + numberOfXNodes1 + 1           
        uNodes2D = [localNodes2D[localNodeIdx100],localNodes2D[localNodeIdx010],localNodes2D[localNodeIdx001],localNodes2D[3],localNodes2D[4],localNodes2D[5]]
    elif (uInterpolation == CUBIC_SIMPLEX):
        if (subElementNumber == 1):
            localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx100] + 2*numberOfXNodes1 + numberOfXNodes2 + 3
            localNodes2D[localNodeIdx001] = localNodes2D[localNodeIdx100] + 2*numberOfXNodes1 + numberOfXNodes2
            localNodes2D[3] = localNodes2D[localNodeIdx100] + numberOfXNodes1 + 1
            localNodes2D[4] = localNodes2D[localNodeIdx100] + 2*numberOfXNodes1 + 2
            localNodes2D[5] = localNodes2D[localNodeIdx100] + 2*numberOfXNodes1 + numberOfXNodes2 + 2
            localNodes2D[6] = localNodes2D[localNodeIdx100] + 2*numberOfXNodes1 + numberOfXNodes2 + 1
            localNodes2D[7] = localNodes2D[localNodeIdx100] + 2*numberOfXNodes1
            localNodes2D[8] = localNodes2D[localNodeIdx100] + numberOfXNodes1
            localNodes2D[9] = localNodes2D[localNodeIdx100] + 2*numberOfXNodes1 + 1
        else:
            localNodes2D[localNodeIdx010] = localNodes2D[localNodeIdx100] + 3
            localNodes2D[localNodeIdx001] = localNodes2D[localNodeIdx100] + 2*numberOfXNodes1+numberOfXNodes2 + 3
            localNodes2D[3] = localNodes2D[localNodeIdx100] + 1
            localNodes2D[4] = localNodes2D[localNodeIdx100] + 2
            localNodes2D[5] = localNodes2D[localNodeIdx00] + numberOfXNodes1 + 3
            localNodes2D[6] = localNodes2D[localNodeIdx00] + 2*numberOfXNodes1 + 3
            localNodes2D[7] = localNodes2D[localNodeIdx00] + 2*numberOfXNodes1 + 2
            localNodes2D[8] = localNodes2D[localNodeIdx00] + numberOfXNodes1 + 1
            localNodes2D[9] = localNodes2D[localNodeIdx00] + numberOfXNodes1 + 2
        uNodes2D = [localNodes2D[localNodeIdx100],localNodes2D[localNodeIdx010],localNodes2D[localNodeIdx001], \
                    localNodes2D[3],localNodes2D[4],localNodes2D[5],localNodes2D[6],localNodes2D[7],localNodes2D[8], \
                    localNodes2D[9]]
    else:
        print('Invalid u interpolation')
        exit()
    if (simplex):
        pNodes2D = [localNodes2D[localNodeIdx100],localNodes2D[localNodeIdx010],localNodes2D[localNodeIdx001]]
    else:
        pNodes2D = [localNodes2D[localNodeIdx00],localNodes2D[localNodeIdx10],localNodes2D[localNodeIdx01],localNodes2D[localNodeIdx11]]          
    if (debugLevel > 2):
        print('    Element %8d' %(elementNumber))
        print('      Sub-element %8d' %(subElementNumber))
        print('        U Nodes: '+str(uNodes2D))
        print('        P Nodes: '+str(pNodes2D))
    return uNodes2D,pNodes2D

def GetElementNodes1D(elementNumber,localNodes1D):
    localNodes1D[localNodeIdx1] = localNodes1D[localNodeIdx0]+(numberOfNodesXi-1)
    if (uInterpolation == LINEAR_LAGRANGE or uInterpolation == CUBIC_HERMITE or uInterpolation == LINEAR_SIMPLEX):
        uNodes1D = [localNodes1D[localNodeIdx0],localNodes1D[localNodeIdx1]]
    elif (uInterpolation == QUADRATIC_LAGRANGE or uInterpolation == QUADRATIC_SIMPLEX):
        localNodes1D[1] = localNodes1D[localNodeIdx0]+1
        uNodes1D = [localNodes1D[localNodeIdx0],localNodes1D[1],localNodes1D[localNodeIdx1]]
    elif (uInterpolation == CUBIC_LAGRANGE or uInterpolation == CUBIC_SIMPLEX):
        localNodes1D[1] = localNodes1D[localNodeIdx0]+1
        localNodes1D[2] = localNodes1D[1]+1
        uNodes1D = [localNodes1D[localNodeIdx0],localNodes1D[1],localNodes1D[2],localNodes1D[localNodeIdx1]]
    else:
        print('Invalid u interpolation')
        exit()                    
    pNodes1D = [localNodes1D[localNodeIdx0],localNodes1D[localNodeIdx1]]
    if (debugLevel > 2):
        print('    Element %8d' %(elementNumber))
        print('      U Nodes: '+str(uNodes1D))
        print('      P Nodes: '+str(pNodes1D))
    return uNodes1D,pNodes1D

def SetNodeParameters2D(nodeNumber,field,xPosition,yPosition):
    field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
    field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
    if (debugLevel > 2):
        print('      Node        %d:' % (nodeNumber))
        print('         Position         = [ %.2f, %.2f ]' % (xPosition,yPosition))                 
    if (uInterpolation == CUBIC_HERMITE):
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,1.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,0.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,0.0)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,1.0)
        if (debugLevel > 2):
            print('        S1 derivative    = [ %.2f, %.2f ]' % (1.0,0.0))                 
            print('        S2 derivative    = [ %.2f, %.2f ]' % (0.0,1.0))                 
            
def SetNodeParameters1D(nodeNumber,field,xPosition,yPosition,xTangent,yTangent):
    field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
    field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                   1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
    if (debugLevel > 2):
        print('      Node        %d:' % (nodeNumber))
        print('         Position         = [ %.2f, %.2f ]' % (xPosition,yPosition))                 
    if (uInterpolation == CUBIC_HERMITE):
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xTangent)
        field.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                       1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,yTangent)
        if (debugLevel > 2):
            print('        S1 derivative    = [ %.2f, %.2f ]' % (xTangent,yTangent))                 
            
#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import numpy,csv,time,sys,os,pdb
from opencmiss.iron import iron

# Ensure output directories exist
if not os.path.exists('./output'):
    os.makedirs('./output')
if not os.path.exists('./output/Fluid'):
    os.makedirs('./output/Fluid')
if not os.path.exists('./output/Solid'):
    os.makedirs('./output/Solid')
if not os.path.exists('./output/Interface'):
    os.makedirs('./output/Interface')

# Diagnostics
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#iron.ErrorHandlingModeSet(iron.ErrorHandlingModes.TRAP_ERROR)
iron.OutputSetOn("Testing")

# Get the computational nodes info
computationEnvironment = iron.ComputationEnvironment()
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()
          
#================================================================================================================================
#  Initial Data & Default Values
#================================================================================================================================

# (NONE/TIMING/MATRIX/ELEMENT_MATRIX/NODAL_MATRIX)
fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
#fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
fluidEquationsOutputType = iron.EquationsOutputTypes.NONE
#fluidEquationsOutputType = iron.EquationsOutputTypes.TIMING
#fluidEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#fluidEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
solidEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
#solidEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
solidEquationsOutputType = iron.EquationsOutputTypes.NONE
#solidEquationsOutputType = iron.EquationsOutputTypes.TIMING
#solidEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#solidEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
movingMeshEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
#movingMeshEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
movingMeshEquationsOutputType = iron.EquationsOutputTypes.NONE
#movingMeshEquationsOutputType = iron.EquationsOutputTypes.TIMING
#movingMeshEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#movingMeshEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
interfaceConditionOutputType = iron.InterfaceConditionOutputTypes.NONE
#interfaceConditionOutputType = iron.InterfaceConditionOutputTypes.PROGRESS
interfaceEquationsOutputType = iron.EquationsOutputTypes.NONE
#interfaceEquationsOutputType = iron.EquationsOutputTypes.TIMING
#interfaceEquationsOutputType = iron.EquationsOutputTypes.PROGRESS
#interfaceEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#interfaceEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
# (NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
movingMeshLinearSolverOutputType = iron.SolverOutputTypes.NONE
#movingMeshLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#movingMeshLinearSolverOutputType = iron.SolverOutputTypes.MATRIX
#fsiDynamicSolverOutputType = iron.SolverOutputTypes.NONE
fsiDynamicSolverOutputType = iron.SolverOutputTypes.MONITOR
#fsiDynamicSolverOutputType = iron.SolverOutputTypes.MATRIX
#fsiNonlinearSolverOutputType = iron.SolverOutputTypes.NONE
fsiNonlinearSolverOutputType = iron.SolverOutputTypes.MONITOR
#fsiNonlinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fsiNonlinearSolverOutputType = iron.SolverOutputTypes.MATRIX
#fsiLinearSolverOutputType = iron.SolverOutputTypes.NONE
fsiLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fsiLinearSolverOutputType = iron.SolverOutputTypes.MATRIX

if (setupOutput):
    print('SUMMARY')
    print('=======')
    print(' ')
    print('  Temporal parameters')
    print('  -------------------')
    print(' ')
    print('  Start time:     %.3f' % (startTime))
    print('  Stop time:      %.3f' % (stopTime))
    print('  Time increment: %.5f' % (timeStep))
    print(' ')
    print('  Material parameters')
    print('  -------------------')
    print(' ')
    if (problemType != SOLID):
        print('    Fluid:')
        print('      Dynamic viscosity: {0:.3f}'.format(fluidDynamicViscosity))
        print('      Density: {0:.3f}'.format(fluidDensity))
        print(' ')
    if (problemType != FLUID):
        print('    Solid:')
        print('      Density: {0:.3f}'.format(solidDensity))
        print('      Mooney Rivlin 1: {0:.3f}'.format(mooneyRivlin1))
        print('      Mooney Rivlin 2: {0:.3f}'.format(mooneyRivlin2))
        print(' ')
    print('  Mesh parameters')
    print('  -------------------')
    print(' ')
    if (uInterpolation == LINEAR_LAGRANGE):
        print('    U interpolation: LINEAR_LAGRANGE')
    elif (uInterpolation == QUADRATIC_LAGRANGE):
        print('    U interpolation: QUADRATIC_LAGRANGE')
    elif (uInterpolation == CUBIC_LAGRANGE):
        print('    U interpolation: CUBIC_LAGRANGE')
    elif (uInterpolation == CUBIC_HERMITE):
        print('    U interpolation: CUBIC_HERMITE')
    elif (uInterpolation == LINEAR_SIMPLEX):
        print('    U interpolation: LINEAR_SIMPLEX')
    elif (uInterpolation == QUADRATIC_SIMPLEX):
        print('    U interpolation: QUADRATIC_SIMPLEX')
    elif (uInterpolation == CUBIC_SIMPLEX):
        print('    U interpolation: CUBIC_SIMPLEX')
    else:
        print('Invalid u interpolation')
        exit()            
    if (problemType != SOLID):
        print('    Fluid:')
        print('      Number of X1 elements: {0:d}'.format(numberOfFluidX1Elements))
        print('      Number of X2 elements: {0:d}'.format(numberOfFluidX2Elements))
        print('      Number of Y  elements: {0:d}'.format(numberOfFluidYElements))
        print('      Number of nodes: {0:d}'.format(numberOfFluidNodes))
        print('      Number of elements: {0:d}'.format(numberOfFluidElements))
    if (problemType != FLUID):
        print('    Solid:')
        print('      Number of X elements: {0:d}'.format(numberOfSolidXElements))
        print('      Number of Y elements: {0:d}'.format(numberOfSolidYElements))
        print('      Number of nodes: {0:d}'.format(numberOfSolidNodes))
        print('      Number of elements: {0:d}'.format(numberOfSolidElements))
    if (problemType == FSI):
        print('    Interface:')
        print('      Number of nodes: {0:d}'.format(numberOfInterfaceNodes))
        print('      Number of elements: {0:d}'.format(numberOfInterfaceElements))
        
#================================================================================================================================
#  Coordinate Systems
#================================================================================================================================

if (progressDiagnostics):
    print(' ')
    print('Coordinate systems ...')

if (problemType != FLUID):
    # Create a RC coordinate system for the solid region
    solidCoordinateSystem = iron.CoordinateSystem()
    solidCoordinateSystem.CreateStart(solidCoordinateSystemUserNumber)
    solidCoordinateSystem.DimensionSet(numberOfDimensions)
    solidCoordinateSystem.CreateFinish()
if (problemType != SOLID):
    # Create a RC coordinate system for the fluid region
    fluidCoordinateSystem = iron.CoordinateSystem()
    fluidCoordinateSystem.CreateStart(fluidCoordinateSystemUserNumber)
    fluidCoordinateSystem.DimensionSet(numberOfDimensions)
    fluidCoordinateSystem.CreateFinish()
if (problemType == FSI):
    # Create a RC coordinate system for the interface region
    interfaceCoordinateSystem = iron.CoordinateSystem()
    interfaceCoordinateSystem.CreateStart(interfaceCoordinateSystemUserNumber)
    interfaceCoordinateSystem.DimensionSet(numberOfDimensions)
    interfaceCoordinateSystem.CreateFinish()
    
if (progressDiagnostics):
    print('Coordinate systems ... Done')
    
#================================================================================================================================
#  Regions
#================================================================================================================================

if (progressDiagnostics):
    print('Regions ...')

if (problemType != FLUID):
    # Create a solid region
    solidRegion = iron.Region()
    solidRegion.CreateStart(solidRegionUserNumber,iron.WorldRegion)
    solidRegion.label = 'SolidRegion'
    solidRegion.coordinateSystem = solidCoordinateSystem
    solidRegion.CreateFinish()
if (problemType != SOLID):
    # Create a fluid region
    fluidRegion = iron.Region()
    fluidRegion.CreateStart(fluidRegionUserNumber,iron.WorldRegion)
    fluidRegion.label = 'FluidRegion'
    fluidRegion.coordinateSystem = fluidCoordinateSystem
    fluidRegion.CreateFinish()
    
if (progressDiagnostics):
    print('Regions ... Done')

#================================================================================================================================
#  Bases
#================================================================================================================================

if (progressDiagnostics):
    print('Basis functions ...')
          
pBasis = iron.Basis()
pBasis.CreateStart(pBasisUserNumber)
pBasis.NumberOfXiSet(numberOfDimensions)
if (simplex):
    pBasis.TypeSet(iron.BasisTypes.SIMPLEX)
    pBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*numberOfDimensions)
    pBasis.QuadratureOrderSet(gaussOrder)
else:
    pBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    pBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfDimensions)
    pBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfDimensions)
pBasis.CreateFinish()

uBasis = iron.Basis()
uBasis.CreateStart(uBasisUserNumber)
uBasis.NumberOfXiSet(numberOfDimensions)
if (simplex):
    uBasis.TypeSet(iron.BasisTypes.SIMPLEX)
    if (uInterpolation == LINEAR_SIMPLEX):
        uBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*numberOfDimensions)
    elif (uInterpolation == QUADRATIC_SIMPLEX):
        uBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX]*numberOfDimensions)
    elif (uInterpolation == CUBIC_SIMPLEX):
        uBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_SIMPLEX]*numberOfDimensions)
    else:
        print('Invalid u interpolation for simplex')
        exit()
    uBasis.QuadratureOrderSet(gaussOrder)
else:
    uBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    if (uInterpolation == LINEAR_LAGRANGE):
        uBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfDimensions)
    elif (uInterpolation == QUADRATIC_LAGRANGE):
        uBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfDimensions)
    elif (uInterpolation == CUBIC_LAGRANGE):
        uBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*numberOfDimensions)
    elif (uInterpolation == CUBIC_HERMITE):
        uBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*numberOfDimensions)
    else:
        print('Invalid u interpolation for non simplex')
        exit()
    uBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfDimensions)
uBasis.CreateFinish()

if (problemType == FSI):
    interfaceBasis = iron.Basis()
    interfaceBasis.CreateStart(interfaceBasisUserNumber)
    interfaceBasis.NumberOfXiSet(numberOfInterfaceDimensions)
    if (simplex):
        interfaceBasis.TypeSet(iron.BasisTypes.SIMPLEX)
        if (uInterpolation == LINEAR_SIMPLEX):
            interfaceBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*numberOfInterfaceDimensions)
        elif (uInterpolation == QUADRATIC_SIMPLEX):
            interfaceBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX]*numberOfInterfaceDimensions)
        elif (uInterpolation == CUBIC_SIMPLEX):
            interfaceBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_SIMPLEX]*numberOfInterfaceDimensions)
        else:
            print('Invalid u interpolation for simplex')
            exit()
        interfaceBasis.QuadratureOrderSet(gaussOrder)
    else:
        interfaceBasis.Type(iron.BasisTypes.LAGRANGE_HERMITE_TP)
        if (uInterpolation == LINEAR_LAGRANGE):
            interfaceBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfInterfaceDimensions)
        elif (uInterpolation == QUADRATIC_LAGRANGE):
            interfaceBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfInterfaceDimensions)
        elif (uInterpolation == CUBIC_LAGRANGE):
            interfaceBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*numberOfInterfaceDimensions)
        elif (uInterpolation == CUBIC_HERMITE):
            interfaceBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*numberOfInterfaceDimensions)
        else:
            print('Invalid u interpolation for non simplex')
            exit()
        interfaceBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfInterfaceDimensions)
    interfaceBasis.CreateFinish()

if (progressDiagnostics):
    print('Basis functions ... Done')
  
#================================================================================================================================
#  Mesh
#================================================================================================================================

if (progressDiagnostics):
    print('Meshes ...')    
                   
pNodes2D = [0]*4
uNodes2D = [0]*numberOfLocalNodes
localNodes2D = [0]*numberOfLocalNodes

if (problemType != FLUID):
    solidNodes = iron.Nodes()
    solidNodes.CreateStart(solidRegion,numberOfSolidNodes)
    solidNodes.CreateFinish()

    solidMesh = iron.Mesh()
    solidMesh.CreateStart(solidMeshUserNumber,solidRegion,numberOfDimensions)
    solidMesh.NumberOfElementsSet(numberOfSolidElements)
    solidMesh.NumberOfComponentsSet(2)

    solidUElements = iron.MeshElements()
    solidUElements.CreateStart(solidMesh,1,uBasis)
    
    solidPElements = iron.MeshElements()
    solidPElements.CreateStart(solidMesh,2,pBasis)
                
    # Solid mesh elements
    if (debugLevel > 2):
        print('  Solid Elements:')
    for yElementIdx in range(1,numberOfSolidYElements+1):
        for xElementIdx in range(1,numberOfSolidXElements+1):
            for subElementIdx in range(1,numberOfSubElements+1):
                elementNumber = subElementIdx+((xElementIdx-1)+(yElementIdx-1)*numberOfSolidXElements)*numberOfSubElements
                localNodes2D[localNodeIdx00]=(xElementIdx-1)*(numberOfNodesXi-1)+1+ \
                                              (yElementIdx-1)*(numberOfNodesXi-1)*(numberOfSolidXNodes)
                [uNodes2D,pNodes2D] = GetElementNodes2D(elementNumber,subElementIdx,localNodes2D,numberOfSolidXNodes,numberOfSolidXNodes)
                solidUElements.NodesSet(elementNumber,uNodes2D)
                solidPElements.NodesSet(elementNumber,pNodes2D)

    solidUElements.CreateFinish()
    solidPElements.CreateFinish()

    solidMesh.CreateFinish()

if (problemType != SOLID):
    fluidNodes = iron.Nodes()
    fluidNodes.CreateStart(fluidRegion,numberOfFluidNodes)
    fluidNodes.CreateFinish()

    fluidMesh = iron.Mesh()
    fluidMesh.CreateStart(fluidMeshUserNumber,fluidRegion,numberOfDimensions)
    fluidMesh.NumberOfElementsSet(numberOfFluidElements)
    fluidMesh.NumberOfComponentsSet(2)

    fluidUElements = iron.MeshElements()
    fluidUElements.CreateStart(fluidMesh,1,uBasis)
    
    fluidPElements = iron.MeshElements()
    fluidPElements.CreateStart(fluidMesh,2,pBasis)
                        
    # Fluid mesh elements
    if (debugLevel > 2):
        print('  Fluid Elements:')
    for yElementIdx in range(1,numberOfSolidYElements+1):
        # Elements to the left of the solid
        for xElementIdx in range(1,numberOfFluidX1Elements+1):
            for subElementIdx in range(1,numberOfSubElements+1):
                elementNumber = subElementIdx+(xElementIdx-1)*numberOfSubElements+(yElementIdx-1)*numberOfFluidXElements2
                localNodes2D[localNodeIdx00] = (xElementIdx-1)*(numberOfNodesXi-1)+1+\
                                               (yElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidXNodes2
                [uNodes2D,pNodes2D] = GetElementNodes2D(elementNumber,subElementIdx,localNodes2D,numberOfFluidXNodes2,numberOfFluidXNodes2)
                fluidUElements.NodesSet(elementNumber,uNodes2D)
                fluidPElements.NodesSet(elementNumber,pNodes2D)
        # Elements to the right of the solid
        for xElementIdx in range(1,numberOfFluidX2Elements+1):
            for subElementIdx in range(1,numberOfSubElements+1):
                elementNumber = numberOfFluidX1Elements*numberOfSubElements+subElementIdx+(xElementIdx-1)*numberOfSubElements+\
                                (yElementIdx-1)*numberOfFluidXElements2
                localNodes2D[localNodeIdx00] = numberOfFluidX1Nodes+(xElementIdx-1)*(numberOfNodesXi-1)+1+\
                                               (yElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidXNodes2
                if(yElementIdx == numberOfSolidYElements):
                    [uNodes2D,pNodes2D] = GetElementNodes2D(elementNumber,subElementIdx,localNodes2D,numberOfFluidXNodes2,numberOfFluidXNodes1)
                else: 
                    [uNodes2D,pNodes2D] = GetElementNodes2D(elementNumber,subElementIdx,localNodes2D,numberOfFluidXNodes2,numberOfFluidXNodes2)
                fluidUElements.NodesSet(elementNumber,uNodes2D)
                fluidPElements.NodesSet(elementNumber,pNodes2D)
    # Elements above the solid
    for yElementIdx in range(1,numberOfFluidYElements+1):
        for xElementIdx in range(1,numberOfFluidX1Elements+numberOfSolidXElements+numberOfFluidX2Elements+1):
            for subElementIdx in range(1,numberOfSubElements+1):
                elementNumber = numberOfFluidXElements2*numberOfSolidYElements+subElementIdx+(xElementIdx-1)*numberOfSubElements+\
                                (yElementIdx-1)*numberOfFluidXElements1
                localNodes2D[localNodeIdx00] = numberOfFluidXNodes2*(numberOfSolidYNodes-1)+ \
                                            (xElementIdx-1)*(numberOfNodesXi-1)+1+\
                                            (yElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidXNodes1
                [uNodes2D,pNodes2D] = GetElementNodes2D(elementNumber,subElementIdx,localNodes2D,numberOfFluidXNodes1,numberOfFluidXNodes1)
                fluidUElements.NodesSet(elementNumber,uNodes2D)
                fluidPElements.NodesSet(elementNumber,pNodes2D)

    fluidUElements.CreateFinish()
    fluidPElements.CreateFinish()

    fluidMesh.CreateFinish()

if (progressDiagnostics):
    print('Meshes ... Done')    

#================================================================================================================================
#  Interface
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Interface ...')
    
        # Create an interface between the two meshes
        interface = iron.Interface()
        interface.CreateStart(interfaceUserNumber,iron.WorldRegion)
        interface.LabelSet('Interface')
        # Add in the two meshes
        solidMeshIndex = interface.MeshAdd(solidMesh)
        fluidMeshIndex = interface.MeshAdd(fluidMesh)
        interface.CoordinateSystemSet(interfaceCoordinateSystem)
        interface.CreateFinish()
        
    if (progressDiagnostics):
        print('Interface ... Done')
            
#================================================================================================================================
#  Interface Mesh
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Interface Mesh ...')
    
    pNodes1D = [0]*2
    uNodes1D = [0]*numberOfLocalInterfaceNodes
    localNodes1D = [0]*numberOfLocalInterfaceNodes

    # Create an interface mesh
    InterfaceNodes = iron.Nodes()
    InterfaceNodes.CreateStartInterface(interface,numberOfInterfaceNodes)
    InterfaceNodes.CreateFinish()
    
    interfaceMesh = iron.Mesh()
    interfaceMesh.CreateStartInterface(interfaceMeshUserNumber,interface,numberOfInterfaceDimensions)
    interfaceMesh.NumberOfElementsSet(numberOfInterfaceElements)
    interfaceMesh.NumberOfComponentsSet(1)
    
    interfaceElements = iron.MeshElements()
    interfaceElements.CreateStart(interfaceMesh,1,interfaceBasis)
        
    if (debugLevel > 2):
        print('  Interface Elements:')
    elementNumber = 0
    for interfaceElementIdx in range(1,numberOfSolidXElements + 2*numberOfSolidYElements + 1):
        elementNumber = elementNumber + 1
        localNodes1D[localNodeIdx00] = (interfaceElementIdx-1)*(numberOfNodesXi-1)+1
        [uNodes1D,pNodes1D] = GetElementNodes1D(elementNumber,localNodes1D)
        interfaceElements.NodesSet(elementNumber,uNodes1D)

    interfaceElements.CreateFinish()

    interfaceMesh.CreateFinish()

    if (progressDiagnostics):
        print('Interface Mesh ... Done')
    

#================================================================================================================================
#  Mesh Connectivity
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Interface Mesh Connectivity ...')

    # Couple the interface meshes
    interfaceMeshConnectivity = iron.InterfaceMeshConnectivity()
    interfaceMeshConnectivity.CreateStart(interface,interfaceMesh)
    interfaceMeshConnectivity.BasisSet(interfaceBasis)
        
    interfaceElementNumber = 0
    interfaceNodes = [0]*(numberOfInterfaceNodes)
    solidNodes = [0]*(numberOfInterfaceNodes)
    fluidNodes = [0]*(numberOfInterfaceNodes)
    localInterfaceNodes = [0]*numberOfNodesXi
    localSolidNodes = [0]*numberOfNodesXi
    localFluidNodes = [0]*numberOfNodesXi
    # Left edge of solid
    for interfaceElementIdx in range(1,numberOfSolidYElements+1):
        interfaceElementNumber = interfaceElementNumber + 1
        if (debugLevel > 2):
            print('  Interface Element %8d:' % (interfaceElementNumber))        
        solidElementNumber = (interfaceElementIdx - 1)*numberOfSolidXElements*numberOfSubElements + 1
        fluidElementNumber = numberOfFluidX1Elements*numberOfSubElements+(interfaceElementIdx - 1)*numberOfFluidXElements2
        # Map interface elements
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,solidMeshIndex,solidElementNumber)
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber)
        if (debugLevel > 2):
            print('    Solid Element %8d; Fluid Element %8d' % (solidElementNumber,fluidElementNumber))        
        localInterfaceNodes[0] = (interfaceElementIdx-1)*(numberOfNodesXi-1) + 1
        localSolidNodes[0] = (interfaceElementIdx-1)*(numberOfNodesXi-1)*numberOfSolidXNodes+1
        localFluidNodes[0] = numberOfFluidX1Nodes + (interfaceElementIdx-1)*(numberOfNodesXi-1)*numberOfFluidXNodes2
        if (uInterpolation == QUADRATIC_LAGRANGE or uInterpolation == QUADRATIC_SIMPLEX):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] + numberOfSolidXNodes
            localFluidNodes[1] = localFluidNodes[0] + numberOfFluidXNodes2
        elif (uInterpolation == CUBIC_LAGRANGE or uInterpolation == CUBIC_SIMPLEX):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] + numberOfSolidXNodes
            localFluidNodes[1] = localFluidNodes[0] + numberOfFluidXNodes2
            localInterfaceNodes[2] = localInterfaceNodes[1]+1
            localSolidNodes[2] = localSolidNodes[1] + numberOfSolidXNodes
            localFluidNodes[2] = localFluidNodes[1] + numberOfFluidXNodes2        
        localInterfaceNodes[numberOfNodesXi-1] = localInterfaceNodes[0] + (numberOfNodesXi-1)
        localSolidNodes[numberOfNodesXi-1] = localSolidNodes[0] + (numberOfNodesXi-1)*numberOfSolidXNodes
        localFluidNodes[numberOfNodesXi-1] = localFluidNodes[0] + (numberOfNodesXi-1)*numberOfFluidXNodes2
        # Map interface xi
        for localNodeIdx in range(0,numberOfNodesXi):
            xi=float(localNodeIdx)/float(numberOfNodesXi-1)
            if (simplex):
                solidXi = [xi,1.0]
                fluidXi = [1.0,xi]
            else:
                solidXi = [0.0,xi]
                fluidXi = [1.0,xi]
            interfaceNodes[localInterfaceNodes[localNodeIdx]-1]=localInterfaceNodes[localNodeIdx]
            solidNodes[localInterfaceNodes[localNodeIdx]-1]=localSolidNodes[localNodeIdx]
            fluidNodes[localInterfaceNodes[localNodeIdx]-1]=localFluidNodes[localNodeIdx]
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,solidMeshIndex,solidElementNumber,localNodeIdx+1,1,solidXi)
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber,localNodeIdx+1,1,fluidXi)
            if (debugLevel > 2):
                print('    Local node    %8d:' % (localNodeIdx+1))        
                print('      Interface node    %8d:' % (localInterfaceNodes[localNodeIdx]))        
                print('      Solid node        %8d; Solid xi = [%.2f, %.2f ]' % (localSolidNodes[localNodeIdx],solidXi[0],solidXi[1]))
                print('      Fluid node        %8d; Fluid xi = [%.2f, %.2f ]' % (localFluidNodes[localNodeIdx],fluidXi[0],fluidXi[1]))
    # Top edge of solid
    for interfaceElementIdx in range(1,numberOfSolidXElements+1):
        interfaceElementNumber = interfaceElementNumber + 1
        if (debugLevel > 2):
            print('  Interface Element %8d:' % (interfaceElementNumber))        
        solidElementNumber = 1+(interfaceElementIdx-1)*numberOfSubElements + numberOfSolidXElements*numberOfSubElements*(numberOfSolidYElements-1)
        fluidElementNumber = numberOfSubElements + (interfaceElementIdx-1)*numberOfSubElements + \
                             (numberOfFluidX1Elements*numberOfSubElements + numberOfFluidXElements2*numberOfSolidYElements)
        # Map interface elements
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,solidMeshIndex,solidElementNumber)
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber)
        if (debugLevel > 2):
            print('    Solid Element %8d; Fluid Element %8d' % (solidElementNumber,fluidElementNumber))        
        localInterfaceNodes[0] = (interfaceElementIdx-1)*(numberOfNodesXi-1) + numberOfSolidYNodes
        localSolidNodes[0] = (interfaceElementIdx-1)*(numberOfNodesXi-1) + 1 + numberOfSolidXNodes*(numberOfSolidYNodes-1)        
        localFluidNodes[0] = (interfaceElementIdx-1)*(numberOfNodesXi-1) + numberOfFluidX1Nodes + numberOfFluidXNodes2*(numberOfSolidYNodes-1)
        if (uInterpolation == QUADRATIC_LAGRANGE or uInterpolation == QUADRATIC_SIMPLEX):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] + 1
            localFluidNodes[1] = localFluidNodes[0] + 1
        elif (uInterpolation == CUBIC_LAGRANGE or uInterpolation == CUBIC_SIMPLEX):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] + 1
            localFluidNodes[1] = localFluidNodes[0] + 1
            localInterfaceNodes[2] = localInterfaceNodes[1]+1
            localSolidNodes[2] = localSolidNodes[1] + 1
            localFluidNodes[2] = localFluidNodes[1] + 1
        localInterfaceNodes[numberOfNodesXi-1] = localInterfaceNodes[0] + (numberOfNodesXi-1)
        localSolidNodes[numberOfNodesXi-1] = localSolidNodes[0] + (numberOfNodesXi-1)
        localFluidNodes[numberOfNodesXi-1] = localFluidNodes[0] + (numberOfNodesXi-1)        
        # Map interface xi
        for localNodeIdx in range(0,numberOfNodesXi):
            xi=float(localNodeIdx)/float(numberOfNodesXi-1)
            if (simplex):
                solidXi = [1.0,1.0-xi]
                fluidXi = [xi,1.0-xi]
            else:
                solidXi = [xi,1.0]
                fluidXi = [xi,0.0]
            interfaceNodes[localInterfaceNodes[localNodeIdx]-1]=localInterfaceNodes[localNodeIdx]
            solidNodes[localInterfaceNodes[localNodeIdx]-1]=localSolidNodes[localNodeIdx]
            fluidNodes[localInterfaceNodes[localNodeIdx]-1]=localFluidNodes[localNodeIdx]
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,solidMeshIndex,solidElementNumber,localNodeIdx+1,1,solidXi)
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber,localNodeIdx+1,1,fluidXi)
            if (debugLevel > 2):
                print('    Local node    %8d:' % (localNodeIdx+1))        
                print('      Interface node    %8d:' % (localInterfaceNodes[localNodeIdx]))        
                print('      Solid node        %8d; Solid xi = [%.2f, %.2f ]' % (localSolidNodes[localNodeIdx],solidXi[0],solidXi[1]))
                print('      Fluid node        %8d; Fluid xi = [%.2f, %.2f ]' % (localFluidNodes[localNodeIdx],fluidXi[0],fluidXi[1]))
    # right edge of solid
    for interfaceElementIdx in range(1,numberOfSolidYElements+1):
        if (interfaceElementIdx == 1):
            offset = numberOfSolidXNodes-1
        else:
            offset = 1
        interfaceElementNumber = interfaceElementNumber + 1
        if (debugLevel > 2):
            print('  Interface Element %8d:' % (interfaceElementNumber))
        solidElementNumber = (numberOfSolidYElements - interfaceElementIdx + 1)*numberOfSolidXElements*numberOfSubElements
        fluidElementNumber = (numberOfSolidYElements - interfaceElementIdx)*numberOfFluidXElements2 + \
                             numberOfFluidX1Elements*numberOfSubElements + 1
        # Map interface elements
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,solidMeshIndex,solidElementNumber)
        interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber)
        if (debugLevel > 2):
            print('    Solid Element %8d; Fluid Element %8d' % (solidElementNumber,fluidElementNumber))        
        localInterfaceNodes[0] = (interfaceElementIdx-1)*(numberOfNodesXi-1) + \
                                 (numberOfSolidXElements + numberOfSolidYElements)*(numberOfNodesXi - 1) + 1
        localSolidNodes[0] = (numberOfSolidXElements*(numberOfNodesXi-1) + 1)* \
                             ((numberOfSolidYElements - interfaceElementIdx + 1)*(numberOfNodesXi-1)+1)
        localFluidNodes[0] = numberOfFluidXNodes2*(numberOfSolidYElements - interfaceElementIdx + 1)*(numberOfNodesXi-1) + \
                             numberOfFluidX1Nodes + offset
        if (uInterpolation == QUADRATIC_LAGRANGE or uInterpolation == QUADRATIC_SIMPLEX):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] - numberOfSolidXNodes 
            localFluidNodes[1] = localFluidNodes[0] - numberOfFluidXNodes2 - offset + 1
        elif (uInterpolation == CUBIC_LAGRANGE or uInterpolation == CUBIC_SIMPLEX):
            localInterfaceNodes[1] = localInterfaceNodes[0]+1
            localSolidNodes[1] = localSolidNodes[0] - numberOfSolidXNodes 
            localFluidNodes[1] = localFluidNodes[0] - numberOfFluidXNodes2 - offset + 1
            localInterfaceNodes[2] = localInterfaceNodes[1]+1
            localSolidNodes[2] = localSolidNodes[1] - numberOfSolidXNodes 
            localFluidNodes[2] = localFluidNodes[1] - numberOfFluidXNodes2 - offset + 1
        localInterfaceNodes[numberOfNodesXi-1] = localInterfaceNodes[0] + numberOfNodesXi - 1
        localSolidNodes[numberOfNodesXi-1] = localSolidNodes[0] - (numberOfNodesXi-1)*numberOfSolidXNodes
        localFluidNodes[numberOfNodesXi-1] = localFluidNodes[0] - (numberOfNodesXi-1)*numberOfFluidXNodes2 - offset + 1
        # Map interface xi
        for localNodeIdx in range(0,numberOfNodesXi):
            xi=float(numberOfNodesXi-localNodeIdx-1)/float(numberOfNodesXi-1)
            if (simplex):
                solidXi = [1.0,xi]
                fluidXi = [xi,1.0]
            else:
                solidXi = [1.0,xi]
                fluidXi = [0.0,xi]
            interfaceNodes[localInterfaceNodes[localNodeIdx]-1]=localInterfaceNodes[localNodeIdx]
            solidNodes[localInterfaceNodes[localNodeIdx]-1]=localSolidNodes[localNodeIdx]
            fluidNodes[localInterfaceNodes[localNodeIdx]-1]=localFluidNodes[localNodeIdx]
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,solidMeshIndex,solidElementNumber,localNodeIdx+1,1,solidXi)
            interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber,localNodeIdx+1,1,fluidXi)
            if (debugLevel > 2):
                print('    Local node    %8d:' % (localNodeIdx+1))        
                print('      Interface node    %8d:' % (localInterfaceNodes[localNodeIdx]))        
                print('      Solid node        %8d; Solid xi = [%.2f, %.2f ]' % (localSolidNodes[localNodeIdx],solidXi[0],solidXi[1]))
                print('      Fluid node        %8d; Fluid xi = [%.2f, %.2f ]' % (localFluidNodes[localNodeIdx],fluidXi[0],fluidXi[1]))
    # Map interface nodes
    interfaceMeshConnectivity.NodeNumberSet(interfaceNodes,solidMeshIndex,solidNodes,fluidMeshIndex,fluidNodes)        

    interfaceMeshConnectivity.CreateFinish()

    if (progressDiagnostics):
        print('Interface Mesh Connectivity ... Done')

#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (progressDiagnostics):
    print('Decomposition ...')
    
if (problemType != FLUID):
    # Create a decomposition for the solid mesh
    solidDecomposition = iron.Decomposition()
    solidDecomposition.CreateStart(solidDecompositionUserNumber,solidMesh)
    solidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
    solidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    solidDecomposition.CalculateFacesSet(True)
    solidDecomposition.CreateFinish()

if (problemType != SOLID):
    # Create a decomposition for the fluid mesh
    fluidDecomposition = iron.Decomposition()
    fluidDecomposition.CreateStart(fluidDecompositionUserNumber,fluidMesh)
    fluidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
    fluidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    fluidDecomposition.CalculateFacesSet(True)
    fluidDecomposition.CreateFinish()

if (problemType == FSI):
    # Create a decomposition for the interface mesh
    interfaceDecomposition = iron.Decomposition()
    interfaceDecomposition.CreateStart(interfaceDecompositionUserNumber,interfaceMesh)
    interfaceDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
    interfaceDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
    interfaceDecomposition.CreateFinish()

if (progressDiagnostics):
    print('Decomposition ... Done')
    
#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (progressDiagnostics):
    print('Geometric Field ...')

if (problemType != FLUID):    
    # Start to create a default (geometric) field on the solid region
    solidGeometricField = iron.Field()
    solidGeometricField.CreateStart(solidGeometricFieldUserNumber,solidRegion)
    # Set the decomposition to use
    solidGeometricField.MeshDecompositionSet(solidDecomposition)
    solidGeometricField.meshDecomposition = solidDecomposition
    # Set the scaling to use
    if (uInterpolation == CUBIC_HERMITE):
        solidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    else:
        solidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
    solidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'SolidGeometry')
    # Set the domain to be used by the field components.
    solidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    solidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    # Finish creating the first field
    solidGeometricField.CreateFinish()

if (problemType != SOLID):
    # Start to create a default (geometric) field on the fluid region
    fluidGeometricField = iron.Field()
    fluidGeometricField.CreateStart(fluidGeometricFieldUserNumber,fluidRegion)
    # Set the decomposition to use
    fluidGeometricField.MeshDecompositionSet(fluidDecomposition)
    # Set the scaling to use
    if (uInterpolation == CUBIC_HERMITE):
        fluidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    else:
        fluidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
    fluidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidGeometry')
    # Set the domain to be used by the field components.
    fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    # Finish creating the second field
    fluidGeometricField.CreateFinish()

if (problemType == FSI):
    # Start to create a default (geometric) field on the Interface
    interfaceGeometricField = iron.Field()
    interfaceGeometricField.CreateStartInterface(interfaceGeometricFieldUserNumber,interface)
    # Set the decomposition to use
    interfaceGeometricField.MeshDecompositionSet(interfaceDecomposition)
    # Set the scaling to use
    if (uInterpolation == CUBIC_HERMITE):
        interfaceGeometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    else:
        interfaceGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
    interfaceGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'InterfaceGeometry')
    # Set the domain to be used by the field components.
    interfaceGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    interfaceGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    # Finish creating the first field
    interfaceGeometricField.CreateFinish()

if (progressDiagnostics):
    print('Geometric Field ... Done')
    
if (progressDiagnostics):
    print('Geometric Parameters ...')
    
if (problemType != FLUID):
    # Solid nodes
    if (debugLevel > 2):
        print('  Solid Nodes:')
    for yNodeIdx in range(1,numberOfSolidYNodes+1):
        for xNodeIdx in range(1,numberOfSolidXNodes+1):
            nodeNumber = xNodeIdx+(yNodeIdx-1)*numberOfSolidXNodes
            nodeDomain = solidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                xPosition = fluidX1Size + float(xNodeIdx-1)/float(numberOfSolidXNodes-1)*solidXSize
                yPosition = float(yNodeIdx-1)/float(numberOfSolidYNodes-1)*solidYSize
                SetNodeParameters2D(nodeNumber,solidGeometricField,xPosition,yPosition)
    # Update fields            
    solidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    solidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (problemType != SOLID):                        
    if (debugLevel > 2):
        print('  Fluid Nodes:')
    for yNodeIdx in range(1,numberOfSolidYNodes):
        # Nodes to the left of the solid
        for xNodeIdx in range(1,numberOfFluidX1Nodes+1):
            nodeNumber = xNodeIdx+(yNodeIdx-1)*numberOfFluidXNodes2
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                xPosition = float(xNodeIdx-1)/float(numberOfFluidX1Elements*(numberOfNodesXi-1))*fluidX1Size
                yPosition = float(yNodeIdx-1)/float(numberOfSolidYElements*(numberOfNodesXi-1))*solidYSize
                SetNodeParameters2D(nodeNumber,fluidGeometricField,xPosition,yPosition)
        # Nodes to the right of the solid
        for xNodeIdx in range(1,numberOfFluidX2Nodes+1):
            nodeNumber = xNodeIdx+numberOfFluidX1Nodes+(yNodeIdx-1)*numberOfFluidXNodes2
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                xPosition = fluidX1Size+solidXSize+float(xNodeIdx-1)/float(numberOfFluidX1Nodes-1)*fluidX1Size
                yPosition = float(yNodeIdx-1)/float(numberOfSolidYNodes-1)*solidYSize
                SetNodeParameters2D(nodeNumber,fluidGeometricField,xPosition,yPosition)
    # Nodes to the top of the solid
    for yNodeIdx in range(1,numberOfFluidYNodes+1):
        for xNodeIdx in range(1,numberOfFluidXNodes1+1):
            nodeNumber = numberOfFluidXNodes2*(numberOfSolidYNodes-1)+xNodeIdx+\
                         (yNodeIdx-1)*numberOfFluidXNodes1
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                xPosition = float(xNodeIdx-1)/float(numberOfFluidXNodes1-1)*(fluidX1Size+solidXSize+fluidX2Size)
                yPosition = solidYSize + float(yNodeIdx-1)/float(numberOfFluidYNodes-1)*fluidYSize
                SetNodeParameters2D(nodeNumber,fluidGeometricField,xPosition,yPosition)
    # Update fields            
    fluidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    fluidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (problemType == FSI):
    if (debugLevel > 2):
        print('  Interface Nodes:')
    # Left edge of interface nodes    
    for yNodeIdx in range(1,numberOfSolidYNodes):
        nodeNumber = yNodeIdx
        #nodeDomain = interfaceDecomposition.NodeDomainGet(nodeNumber,1)
        nodeDomain = computationalNodeNumber
        if (nodeDomain == computationalNodeNumber):
            xPosition = fluidX1Size
            yPosition = float(yNodeIdx-1)/float(numberOfSolidYNodes-1)*solidYSize
            SetNodeParameters1D(nodeNumber,interfaceGeometricField,xPosition,yPosition,0.0,1.0)
    # Top edge of interface nodes    
    for xNodeIdx in range(1,numberOfSolidXNodes+1):
        nodeNumber = xNodeIdx+numberOfSolidYNodes-1
        #nodeDomain = interfaceDecomposition.NodeDomainGet(nodeNumber,1)
        nodeDomain = computationalNodeNumber
        if (nodeDomain == computationalNodeNumber):
            xPosition = fluidX1Size+float(xNodeIdx-1)/float(numberOfSolidXNodes-1)*solidXSize
            yPosition = solidYSize
            SetNodeParameters1D(nodeNumber,interfaceGeometricField,xPosition,yPosition,1.0,0.0)
    # Right edge of interface nodes    
    for yNodeIdx in range(1,numberOfSolidYNodes):
        nodeNumber = yNodeIdx+(numberOfSolidYElements+numberOfSolidXElements)*(numberOfNodesXi-1)+1
        #nodeDomain = interfaceDecomposition.NodeDomainGet(nodeNumber,1)
        nodeDomain = computationalNodeNumber
        if (nodeDomain == computationalNodeNumber):
            xPosition = fluidX1Size+solidXSize
            yPosition = solidYSize-float(yNodeIdx)/float(numberOfSolidYNodes-1)*solidYSize
            SetNodeParameters1D(nodeNumber,interfaceGeometricField,xPosition,yPosition,0.0,-1.0)

    # Update fields            
    interfaceGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    interfaceGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Geometric Parameters ... Done')

#================================================================================================================================
#  Equations Set
#================================================================================================================================

if (progressDiagnostics):
    print('Equations Sets ...')

if (problemType != FLUID):
    # Create the equations set for the solid region 
    solidEquationsSetField = iron.Field()
    solidEquationsSet = iron.EquationsSet()
    solidEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                      iron.EquationsSetTypes.FINITE_ELASTICITY,
                                      iron.EquationsSetSubtypes.MOONEY_RIVLIN]
    solidEquationsSet.CreateStart(solidEquationsSetUserNumber,solidRegion,solidGeometricField,
                                  solidEquationsSetSpecification,solidEquationsSetFieldUserNumber,
                                  solidEquationsSetField)
    solidEquationsSet.OutputTypeSet(solidEquationsSetOutputType)
    solidEquationsSet.CreateFinish()
    
if (problemType != SOLID):
    # Create the equations set for the fluid region - ALE Navier-Stokes
    fluidEquationsSetField = iron.Field()
    fluidEquationsSet = iron.EquationsSet()
    if RBS:
        if (problemType == FSI):
            fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                              iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                              iron.EquationsSetSubtypes.ALE_RBS_NAVIER_STOKES]
        else:
            fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                              iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                              iron.EquationsSetSubtypes.TRANSIENT_RBS_NAVIER_STOKES]            
    else:
        if (problemType == FSI):
            fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                              iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                              iron.EquationsSetSubtypes.ALE_NAVIER_STOKES]
        else:
            fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                              iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                              iron.EquationsSetSubtypes.TRANSIENT_NAVIER_STOKES]
        
    fluidEquationsSet.CreateStart(fluidEquationsSetUserNumber,fluidRegion,fluidGeometricField,
                                  fluidEquationsSetSpecification,fluidEquationsSetFieldUserNumber,
                                  fluidEquationsSetField)
    fluidEquationsSet.OutputTypeSet(fluidEquationsSetOutputType)
    fluidEquationsSet.CreateFinish()

    if RBS:
        # Set boundary retrograde flow stabilisation scaling factor (default 0- do not use)
        fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                           iron.FieldParameterSetTypes.VALUES,1,1.0)
        # Set max CFL number (default 1.0)
        fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                           iron.FieldParameterSetTypes.VALUES,2,1.0E20)
        # Set time increment (default 0.0)
        fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                           iron.FieldParameterSetTypes.VALUES,3,timeStep)
        # Set stabilisation type (default 1.0 = RBS)
        fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                           iron.FieldParameterSetTypes.VALUES,4,1.0)
        
if (problemType == FSI):
    # Create the equations set for the moving mesh
    movingMeshEquationsSetField = iron.Field()
    movingMeshEquationsSet = iron.EquationsSet()
    movingMeshEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                           iron.EquationsSetTypes.LAPLACE_EQUATION,
                                           iron.EquationsSetSubtypes.MOVING_MESH_LAPLACE]
    movingMeshEquationsSet.CreateStart(movingMeshEquationsSetUserNumber,fluidRegion,fluidGeometricField,
                                       movingMeshEquationsSetSpecification,movingMeshEquationsSetFieldUserNumber,
                                       movingMeshEquationsSetField)
    movingMeshEquationsSet.OutputTypeSet(movingMeshEquationsSetOutputType)
    movingMeshEquationsSet.CreateFinish()
    
if (progressDiagnostics):
    print('Equations Sets ... Done')


#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (progressDiagnostics):
    print('Dependent Fields ...')

if (problemType != FLUID):
    # Create the equations set dependent field variables for the solid equations set
    solidDependentField = iron.Field()
    solidEquationsSet.DependentCreateStart(solidDependentFieldUserNumber,solidDependentField)
    solidDependentField.VariableLabelSet(iron.FieldVariableTypes.U,'SolidDependent')
    solidDependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'SolidTraction')
    for componentIdx in range(1,numberOfDimensions+1):
        solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
        solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,componentIdx,1)
    solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,numberOfDimensions+1,2)
    solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,2)
    solidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,numberOfDimensions+1,iron.FieldInterpolationTypes.NODE_BASED)
    solidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,iron.FieldInterpolationTypes.NODE_BASED)
    if (uInterpolation == CUBIC_HERMITE):
        solidDependentField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    else:
        solidDependentField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
    solidEquationsSet.DependentCreateFinish()

    # Initialise the solid dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
    for componentIdx in range(1,numberOfDimensions+1):
        solidGeometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,\
                                                                     componentIdx,solidDependentField,iron.FieldVariableTypes.U,
                                                                     iron.FieldParameterSetTypes.VALUES,componentIdx)
    solidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    numberOfDimensions+1,solidPInit)
    
    solidDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    solidDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (problemType != SOLID):
    # Create the equations set dependent field variables for dynamic Navier-Stokes
    fluidDependentField = iron.Field()
    fluidEquationsSet.DependentCreateStart(fluidDependentFieldUserNumber,fluidDependentField)
    fluidDependentField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidDependent')
    # Set the mesh component to be used by the field components.
    for componentIdx in range(1,numberOfDimensions+1):
        fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
        fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,componentIdx,1)
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,numberOfDimensions+1,2)
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,2)
    # fluidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,numberOfDimensions+1,iron.FieldInterpolationTypes.NODE_BASED)
    # fluidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,iron.FieldInterpolationTypes.NODE_BASED)
    # Finish the equations set dependent field variables
    fluidEquationsSet.DependentCreateFinish()

    # Initialise the fluid dependent field
    for componentIdx in range(1,numberOfDimensions+1):
        fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,componentIdx,0.0)
    # Initialise pressure component
    fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                    numberOfDimensions+1,fluidPInit)
    if RBS:
        fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.PRESSURE_VALUES,3,fluidPInit)
        
    fluidDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    fluidDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (problemType == FSI):        
    # Create the equations set dependent field variables for moving mesh
    movingMeshDependentField = iron.Field()
    movingMeshEquationsSet.DependentCreateStart(movingMeshDependentFieldUserNumber,movingMeshDependentField)
    movingMeshDependentField.VariableLabelSet(iron.FieldVariableTypes.U,'MovingMeshDependent')
    # Set the mesh component to be used by the field components.
    for componentIdx in range(1,numberOfDimensions+1):
        movingMeshDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
        movingMeshDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,componentIdx,1)
    # Finish the equations set dependent field variables
    movingMeshEquationsSet.DependentCreateFinish()

    # Initialise dependent field moving mesh
    for ComponentIdx in range(1,numberOfDimensions+1):
        movingMeshDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                             componentIdx,0.0)

    movingMeshDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    movingMeshDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Dependent Fields ... Done')
     
#================================================================================================================================
#  Materials Field
#================================================================================================================================

if (progressDiagnostics):
    print('Materials Fields ...')

if (problemType != FLUID):
    # Create the solid materials field
    solidMaterialsField = iron.Field()
    solidEquationsSet.MaterialsCreateStart(solidMaterialsFieldUserNumber,solidMaterialsField)
    solidMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'SolidMaterials')
    solidMaterialsField.VariableLabelSet(iron.FieldVariableTypes.V,'SolidDensity')
    solidEquationsSet.MaterialsCreateFinish()
    # Set Mooney-Rivlin constants c10 and c01 respectively
    solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,mooneyRivlin1)
    solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,mooneyRivlin2)
    solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,solidDensity)

if (problemType != SOLID):
    # Create the equations set materials field variables for dynamic Navier-Stokes
    fluidMaterialsField = iron.Field()
    fluidEquationsSet.MaterialsCreateStart(fluidMaterialsFieldUserNumber,fluidMaterialsField)
    # Finish the equations set materials field variables
    fluidEquationsSet.MaterialsCreateFinish()
    fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,fluidDynamicViscosity)
    fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,fluidDensity)
    
if (problemType == FSI):    
    # Create the equations set materials field variables for moving mesh
    movingMeshMaterialsField = iron.Field()
    movingMeshEquationsSet.MaterialsCreateStart(movingMeshMaterialsFieldUserNumber,movingMeshMaterialsField)
    # Finish the equations set materials field variables
    movingMeshEquationsSet.MaterialsCreateFinish()

    movingMeshMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,\
                                                         movingMeshKParameter)
   
if (progressDiagnostics):
    print('Materials Fields ... Done')
    
#================================================================================================================================
# Independent Field
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Independent Fields ...')

    # Create fluid mesh velocity independent field 
    fluidIndependentField = iron.Field()
    fluidEquationsSet.IndependentCreateStart(fluidIndependentFieldUserNumber,fluidIndependentField)
    fluidIndependentField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidIndependent')
    # Set the mesh component to be used by the field components.
    for componentIdx in range(1,numberOfDimensions+1):
        fluidIndependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
    # Finish the equations set independent field variables
    fluidEquationsSet.IndependentCreateFinish()
  
    # Create the moving mesh independent field 
    movingMeshIndependentField = iron.Field()
    movingMeshEquationsSet.IndependentCreateStart(movingMeshIndependentFieldUserNumber,movingMeshIndependentField)
    movingMeshIndependentField.VariableLabelSet(iron.FieldVariableTypes.U,'MovingMeshIndependent')
    # Set the mesh component to be used by the field components.
    for componentIdx in range(1,numberOfDimensions+1):
        movingMeshIndependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)    
    # Finish the equations set independent field variables
    movingMeshEquationsSet.IndependentCreateFinish()

    # Initialise independent field moving mesh
    movingMeshIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,movingMeshKParameter)

    if (progressDiagnostics):
        print('Independent Fields ... Done')

#================================================================================================================================
#  Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Equations ...')

if (problemType != FLUID):
    # Solid equations
    solidEquations = iron.Equations()
    solidEquationsSet.EquationsCreateStart(solidEquations)
    solidEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    solidEquations.outputType = solidEquationsOutputType
    solidEquationsSet.EquationsCreateFinish()

if (problemType != SOLID):
    # Fluid equations 
    fluidEquations = iron.Equations()
    fluidEquationsSet.EquationsCreateStart(fluidEquations)
    fluidEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    fluidEquations.outputType = fluidEquationsOutputType
    fluidEquationsSet.EquationsCreateFinish()

if (problemType == FSI):
    # Moving mesh equations
    movingMeshEquations = iron.Equations()
    movingMeshEquationsSet.EquationsCreateStart(movingMeshEquations)
    movingMeshEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    movingMeshEquations.outputType = movingMeshEquationsOutputType
    movingMeshEquationsSet.EquationsCreateFinish()

if (progressDiagnostics):
    print('Equations ... Done')

#================================================================================================================================
#  CellML
#================================================================================================================================

if (progressDiagnostics):
    print('CellML ...')

if (problemType != SOLID):
    # Create CellML equations for the temporal fluid boundary conditions
    bcCellML = iron.CellML()
    bcCellML.CreateStart(bcCellMLUserNumber,fluidRegion)
    bcCellMLIdx = bcCellML.ModelImport("input/exponentialrampupinletbc.cellml")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/A")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/B")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/C")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/x")
    bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/y")
    bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inletx")
    bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inlety")
    bcCellML.CreateFinish()

    # Create CellML <--> OpenCMISS field maps
    bcCellML.FieldMapsCreateStart()
    # Map geometric field to x0 and y0
    bcCellML.CreateFieldToCellMLMap(fluidGeometricField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/x",iron.FieldParameterSetTypes.VALUES)
    bcCellML.CreateFieldToCellMLMap(fluidGeometricField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/y",iron.FieldParameterSetTypes.VALUES)
    # Map fluid velocity to ensure dependent field isn't cleared when the velocities are copied back
    bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/inletx",iron.FieldParameterSetTypes.VALUES)
    bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/inlety",iron.FieldParameterSetTypes.VALUES)
    # Map inletx and inlety to dependent field
    bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inletx",iron.FieldParameterSetTypes.VALUES,
	                            fluidDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
    bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inlety",iron.FieldParameterSetTypes.VALUES,
	                            fluidDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES)
    bcCellML.FieldMapsCreateFinish()


    # Create the CellML models field
    bcCellMLModelsField = iron.Field()
    bcCellML.ModelsFieldCreateStart(bcCellMLModelsFieldUserNumber,bcCellMLModelsField)
    bcCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"BCModelMap")
    bcCellML.ModelsFieldCreateFinish()

    # Only evaluate BC on inlet nodes
    bcCellMLModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0)
    if (debugLevel > 2):
        print('  CellML Boundary Conditions:')
        print('    Inlet Model Set:')
    for yNodeIdx in range(2,numberOfSolidYNodes):
        nodeNumber = (yNodeIdx-1)*numberOfFluidXNodes2+1
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                           1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
    for yNodeIdx in range(1,numberOfFluidYNodes):
        nodeNumber = (numberOfSolidYNodes-1)*numberOfFluidXNodes2 + (yNodeIdx-1)*numberOfFluidXNodes1+1
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                           1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))

    # Create the CellML state field
    bcCellMLStateField = iron.Field()
    bcCellML.StateFieldCreateStart(bcCellMLStateFieldUserNumber,bcCellMLStateField)
    bcCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCState")
    bcCellML.StateFieldCreateFinish()

    # Create the CellML parameters field
    bcCellMLParametersField = iron.Field()
    bcCellML.ParametersFieldCreateStart(bcCellMLParametersFieldUserNumber,bcCellMLParametersField)
    bcCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"BCParameters")
    bcCellML.ParametersFieldCreateFinish()

    # Get the component numbers
    AComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/A")
    BComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/B")
    CComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/C")
    # Set up the parameters field
    bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,AComponentNumber,A)
    bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,BComponentNumber,B)
    bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,CComponentNumber,C)

    # Create the CELL intermediate field
    bcCellMLIntermediateField = iron.Field()
    bcCellML.IntermediateFieldCreateStart(bcCellMLIntermediateFieldUserNumber,bcCellMLIntermediateField)
    bcCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCIntermediate")
    bcCellML.IntermediateFieldCreateFinish()

if (progressDiagnostics):
    print('CellML ... Done')

#================================================================================================================================
#  Interface Condition
#================================================================================================================================

if (problemType == FSI):
    if (progressDiagnostics):
        print('Interface Conditions ...')

    # Create an interface condition between the two meshes
    interfaceCondition = iron.InterfaceCondition()
    interfaceCondition.CreateStart(interfaceConditionUserNumber,interface,interfaceGeometricField)
    # Specify the method for the interface condition
    interfaceCondition.MethodSet(iron.InterfaceConditionMethods.LAGRANGE_MULTIPLIERS)
    # Specify the type of interface condition operator
    interfaceCondition.OperatorSet(iron.InterfaceConditionOperators.SOLID_FLUID)
    # Add in the dependent variables from the equations sets
    interfaceCondition.DependentVariableAdd(solidMeshIndex,solidEquationsSet,iron.FieldVariableTypes.U)
    interfaceCondition.DependentVariableAdd(fluidMeshIndex,fluidEquationsSet,iron.FieldVariableTypes.U)
    # Set the label
    interfaceCondition.LabelSet("FSI Interface Condition")
    # Set the output type
    interfaceCondition.OutputTypeSet(interfaceConditionOutputType)
    # Finish creating the interface condition
    interfaceCondition.CreateFinish()

    if (progressDiagnostics):
        print('Interface Conditions ... Done')

    if (progressDiagnostics):
        print('Interface Lagrange Field ...')
    
    # Create the Lagrange multipliers field
    interfaceLagrangeField = iron.Field()
    interfaceCondition.LagrangeFieldCreateStart(interfaceLagrangeFieldUserNumber,interfaceLagrangeField)
    interfaceLagrangeField.VariableLabelSet(iron.FieldVariableTypes.U,'InterfaceLagrange')
    # Finish the Lagrange multipliers field
    interfaceCondition.LagrangeFieldCreateFinish()
    
    for componentIdx in range(1,numberOfDimensions+1):
        interfaceLagrangeField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,componentIdx,0.0)

        interfaceLagrangeField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        interfaceLagrangeField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

    if (progressDiagnostics):
        print('Interface Lagrange Field ... Done')

    if (progressDiagnostics):
        print('Interface Equations ...')

    # Create the interface condition equations
    interfaceEquations = iron.InterfaceEquations()
    interfaceCondition.EquationsCreateStart(interfaceEquations)
    # Set the interface equations sparsity
    interfaceEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    # Set the interface equations output
    interfaceEquations.outputType = interfaceEquationsOutputType
    # Finish creating the interface equations
    interfaceCondition.EquationsCreateFinish()

    if (progressDiagnostics):
        print('Interface Equations ... Done')

#================================================================================================================================
#  Problem
#================================================================================================================================

if (progressDiagnostics):
    print('Problems ...')

# Create a FSI problem
fsiProblem = iron.Problem()
if (problemType == SOLID):
   fsiProblemSpecification = [iron.ProblemClasses.ELASTICITY,
                              iron.ProblemTypes.FINITE_ELASTICITY,
                              iron.ProblemSubtypes.QUASISTATIC_FINITE_ELASTICITY]
elif (problemType == FLUID):
    if RBS:
        fsiProblemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                                   iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                                   iron.ProblemSubtypes.TRANSIENT_RBS_NAVIER_STOKES]
    else:
        fsiProblemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                                   iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                                   iron.ProblemSubtypes.TRANSIENT_NAVIER_STOKES]
elif (problemType == FSI):
    if RBS:
        fsiProblemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                                   iron.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                                   iron.ProblemSubtypes.FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE]
    else:
        fsiProblemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                                   iron.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                                   iron.ProblemSubtypes.FINITE_ELASTICITY_NAVIER_STOKES_ALE]
        
fsiProblem.CreateStart(fsiProblemUserNumber,fsiProblemSpecification)
fsiProblem.CreateFinish()

if (progressDiagnostics):
    print('Problems ... Done')

#================================================================================================================================
#  Control Loop
#================================================================================================================================

if (progressDiagnostics):
    print('Control Loops ...')

# Create the fsi problem control loop
fsiControlLoop = iron.ControlLoop()
fsiProblem.ControlLoopCreateStart()
fsiProblem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],fsiControlLoop)
fsiControlLoop.LabelSet('TimeLoop')
fsiControlLoop.TimesSet(startTime,stopTime,timeStep)
fsiControlLoop.TimeOutputSet(outputFrequency)
fsiProblem.ControlLoopCreateFinish()

if (progressDiagnostics):
    print('Control Loops ... Done')

#================================================================================================================================
#  Solvers
#================================================================================================================================

if (progressDiagnostics):
    print('Solvers ...')

# Create the problem solver
bcCellMLEvaluationSolver = iron.Solver()
fsiDynamicSolver = iron.Solver()
fsiNonlinearSolver = iron.Solver()
fsiLinearSolver = iron.Solver()
movingMeshLinearSolver = iron.Solver()

fsiProblem.SolversCreateStart()
if (problemType == SOLID):
    # Solvers for growth Finite Elasticity problem
    # Get the BC CellML solver
    fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
    bcCellMLEvaluationSolver.outputType = iron.SolverOutputTypes.PROGRESS
    # Get the nonlinear solver
    fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,fsiNonlinearSolver)
    fsiNonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
    fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
    #fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
    fsiNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
    fsiNonlinearSolver.OutputTypeSet(fsiNonlinearSolverOutputType)
    fsiNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
    fsiNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
    fsiNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
    fsiNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
    # Get the dynamic nonlinear linear solver
    fsiNonlinearSolver.NewtonLinearSolverGet(fsiLinearSolver)
    #fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
    #fsiLinearSolver.LinearIterativeTypeSet(iron.IterativeLinearSolverTypes.GMRES)
    #fsiLinearSolver.LinearIterativeGMRESRestartSet(linearRestartValue)
    #fsiLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
    #fsiLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
    #fsiLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
    #fsiLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
    fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
    fsiLinearSolver.OutputTypeSet(fsiLinearSolverOutputType)
elif (problemType == FLUID):
    # Solvers for coupled FiniteElasticity NavierStokes problem
    # Get the BC CellML solver
    fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
    bcCellMLEvaluationSolver.outputType = iron.SolverOutputTypes.PROGRESS
    # Get the dynamic ALE solver
    fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,fsiDynamicSolver)
    fsiDynamicSolver.OutputTypeSet(fsiDynamicSolverOutputType)
    fsiDynamicSolver.DynamicThetaSet(fsiDynamicSolverTheta)
    # Get the dynamic nonlinear solver
    fsiDynamicSolver.DynamicNonlinearSolverGet(fsiNonlinearSolver)
    fsiNonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
    fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
    #fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
    fsiNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
    fsiNonlinearSolver.OutputTypeSet(fsiNonlinearSolverOutputType)
    fsiNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
    fsiNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
    fsiNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
    fsiNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
    # Get the dynamic nonlinear linear solver
    fsiNonlinearSolver.NewtonLinearSolverGet(fsiLinearSolver)
    #fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
    #fsiLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
    #fsiLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
    #fsiLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
    #fsiLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
    fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
    fsiLinearSolver.OutputTypeSet(fsiLinearSolverOutputType)
elif (problemType == FSI):
    # Solvers for coupled FiniteElasticity NavierStokes problem
    # Get the BC CellML solver
    fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
    bcCellMLEvaluationSolver.outputType = iron.SolverOutputTypes.PROGRESS
    # Get the dynamic ALE solver
    fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,fsiDynamicSolver)
    fsiDynamicSolver.OutputTypeSet(fsiDynamicSolverOutputType)
    fsiDynamicSolver.DynamicThetaSet(fsiDynamicSolverTheta)
    # Get the dynamic nonlinear solver
    fsiDynamicSolver.DynamicNonlinearSolverGet(fsiNonlinearSolver)
    fsiNonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
    #fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
    fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
    fsiNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
    fsiNonlinearSolver.OutputTypeSet(fsiNonlinearSolverOutputType)
    fsiNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
    fsiNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
    fsiNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
    fsiNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
    # Get the dynamic nonlinear linear solver
    fsiNonlinearSolver.NewtonLinearSolverGet(fsiLinearSolver)
    #fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
    #fsiLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
    #fsiLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
    #fsiLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
    #fsiLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
    fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
    fsiLinearSolver.OutputTypeSet(fsiLinearSolverOutputType)
    # Linear solver for moving mesh
    fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,movingMeshLinearSolver)
    movingMeshLinearSolver.OutputTypeSet(movingMeshLinearSolverOutputType)
# Finish the creation of the problem solver
fsiProblem.SolversCreateFinish()

if (progressDiagnostics):
    print('Solvers ... Done')

#================================================================================================================================
#  CellML Equations
#================================================================================================================================

if (progressDiagnostics):
    print('CellML Equations ...')

if (problemType != SOLID):
    # Create CellML equations and add BC equations to the solver
    bcEquations = iron.CellMLEquations()
    fsiProblem.CellMLEquationsCreateStart()
    bcCellMLEvaluationSolver.CellMLEquationsGet(bcEquations)
    bcEquationsIndex = bcEquations.CellMLAdd(bcCellML)
    fsiProblem.CellMLEquationsCreateFinish()

if (progressDiagnostics):
    print('CellML Equations ... Done')

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Solver Equations ...')

# Start the creation of the fsi problem solver equations
fsiProblem.SolverEquationsCreateStart()
# Get the fsi dynamic solver equations
fsiSolverEquations = iron.SolverEquations()
if (problemType == SOLID):
    fsiNonlinearSolver.SolverEquationsGet(fsiSolverEquations)
else:
    fsiDynamicSolver.SolverEquationsGet(fsiSolverEquations)
fsiSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
if (problemType != FLUID):
    fsiSolidEquationsSetIndex = fsiSolverEquations.EquationsSetAdd(solidEquationsSet)
if (problemType != SOLID):
    fsiFluidEquationsSetIndex = fsiSolverEquations.EquationsSetAdd(fluidEquationsSet)
if (problemType == FSI):
    fsiInterfaceConditionIndex = fsiSolverEquations.InterfaceConditionAdd(interfaceCondition)
    # Set the time dependence of the interface matrix to determine the interface matrix coefficient in the solver matrix
    # (basiy position in big coupled matrix system)
    interfaceEquations.MatrixTimeDependenceTypeSet(fsiSolidEquationsSetIndex,True, \
                                                   [iron.InterfaceMatricesTimeDependenceTypes.STATIC,\
                                                    iron.InterfaceMatricesTimeDependenceTypes.FIRST_ORDER_DYNAMIC])
    interfaceEquations.MatrixTimeDependenceTypeSet(fsiFluidEquationsSetIndex,True, \
                                                   [iron.InterfaceMatricesTimeDependenceTypes.STATIC,\
                                                    iron.InterfaceMatricesTimeDependenceTypes.STATIC])
    
    # Create the moving mesh solver equations
    movingMeshSolverEquations = iron.SolverEquations()
    # Get the linear moving mesh solver equations
    movingMeshLinearSolver.SolverEquationsGet(movingMeshSolverEquations)
    movingMeshSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
    # Add in the equations set
    movingMeshEquationsSetIndex = movingMeshSolverEquations.EquationsSetAdd(movingMeshEquationsSet)
    
# Finish the creation of the fsi problem solver equations
fsiProblem.SolverEquationsCreateFinish()

if (progressDiagnostics):
    print('Solver Equations ...')

#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (progressDiagnostics):
    print('Boundary Conditions ...')

# Start the creation of the fsi boundary conditions
fsiBoundaryConditions = iron.BoundaryConditions()
fsiSolverEquations.BoundaryConditionsCreateStart(fsiBoundaryConditions)
if (problemType != FLUID):
    # Set no displacement boundary conditions on the bottom edge of the solid
    if (debugLevel > 2):
        print('  Solid Boundary Conditions:')
        print('    No Displacement Boundary conditions:')
    for xNodeIdx in range(1,numberOfSolidXElements*(numberOfNodesXi-1)+2):
        nodeNumber = xNodeIdx
        nodeDomain = solidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Displacement      = [ %.2f, %.2f ]' % (0.0,0.0))                 
            if (uInterpolation == CUBIC_HERMITE):
                fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debugLevel > 2):
        print('    Reference Solid Pressure Boundary Condition:')
        nodeNumber = numberOfSolidXNodes
        nodeDomain = solidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,solidPRef)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         =   %.2f' % (solidPRef))

if (problemType != SOLID):                
    # Set inlet boundary conditions on the left hand edge
    if (debugLevel > 2):
        print('  Fluid Boundary Conditions:')
        print('    Inlet Boundary conditions:')
    for yNodeIdx in range(2,numberOfSolidYNodes+1):
        nodeNumber = (yNodeIdx-1)*numberOfFluidXNodes2+1
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Velocity         = [ %.2f, %.2f ]' % (0.0,0.0))                 
            if (uInterpolation == CUBIC_HERMITE):
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
    for yNodeIdx in range(1,numberOfFluidYNodes):
        nodeNumber = (numberOfSolidYNodes-1)*numberOfFluidXNodes2 + \
                     (yNodeIdx-1)*numberOfFluidXNodes1+1
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Velocity         = [ %.2f, %.2f ]' % (0.0,0.0))                 
            if (uInterpolation == CUBIC_HERMITE):
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
    # Set outlet boundary conditions on the right hand edge to have zero pressure
    if (debugLevel > 2):
        print('    Outlet Boundary conditions:')
    # Elements to the right of the solid
    for yElementIdx in range(2,numberOfSolidYElements+1):
        nodeNumber = (yElementIdx-1)*(numberOfNodesXi-1)*(numberOfFluidXNodes2)+numberOfFluidXNodes2
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
        if (nodeDomain == computationalNodeNumber):
            if RBS:
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber,numberOfDimensions+1,iron.BoundaryConditionsTypes.PRESSURE,fluidPRef)
            else:
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.DELUDELN,1, \
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber,numberOfDimensions+1,iron.BoundaryConditionsTypes.FIXED,fluidPRef)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         =   %.2f' % (fluidPRef))                 
        if RBS:
            # Set the element normals for outlet stabilisation
            elementNumber = (numberOfFluidX1Elements+numberOfFluidX2Elements)*numberOfSubElements+\
                            (yElementIdx-2)*(numberOfFluidX1Elements+numberOfFluidX2Elements)*numberOfSubElements
            elementDomain = fluidDecomposition.ElementDomainGet(elementNumber)
            if (elementDomain == computationalNodeNumber):
                # Set the outflow normal to (0,0,+1)
                fluidEquationsSetField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.V, \
                                                                   iron.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,5,+1.0)
                fluidEquationsSetField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.V, \
                                                                   iron.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,6,0.0)
                # Set the boundary type
                fluidEquationsSetField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.V, \
                                                                   iron.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,9, \
                                                                   iron.BoundaryConditionsTypes.PRESSURE)                                                
                if (debugLevel > 2):
                    print('      Element     %d:' % (elementNumber))
                    print('         Normal          = [ %.2f, %.2f ]' % (+1.0,0.0))
    # Elements above the solid
    for yElementIdx in range(1,numberOfFluidYElements+2):
        nodeNumber = (numberOfSolidYNodes-1)*numberOfFluidXNodes2 + ((yElementIdx-1)*(numberOfNodesXi-1)+1)*numberOfFluidXNodes1
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
        if (nodeDomain == computationalNodeNumber):
            if RBS:
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber,numberOfDimensions+1,iron.BoundaryConditionsTypes.PRESSURE,fluidPRef)
            else:
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.DELUDELN,1, \
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                              nodeNumber,numberOfDimensions+1,iron.BoundaryConditionsTypes.FIXED,fluidPRef)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         =   %.2f' % (fluidPRef))
        if RBS:
            # Set the element normals for outlet stabilisation
            elementNumber = numberOfFluidXElements2*numberOfSolidYElements+\
                            numberOfFluidXElements1+(yElementIdx-2)*(numberOfFluidXElements1)
            elementDomain = fluidDecomposition.ElementDomainGet(elementNumber)
            if (elementDomain == computationalNodeNumber):
                # Set the outflow normal to (0,0,+1)
                fluidEquationsSetField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.V, \
                                                                   iron.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,5,+1.0)
                fluidEquationsSetField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.V, \
                                                                   iron.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,6,0.0)
                # Set the boundary type
                fluidEquationsSetField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.V, \
                                                                   iron.FieldParameterSetTypes.VALUES, \
                                                                   elementNumber,9, \
                                                                   iron.BoundaryConditionsTypes.PRESSURE)
                if (debugLevel > 2):
                        print('      Element     %d:' % (elementNumber))
                        print('         Normal          = [ %.2f, %.2f ]' % (+1.0,0.0))
                
    # Set no-slip boundary conditions on the bottom edge
    if (debugLevel > 2):
        print('    No-slip Boundary conditions:')
    for xNodeIdx in range(1,numberOfFluidX1Nodes+1):
        nodeNumber = xNodeIdx
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (uInterpolation == CUBIC_HERMITE):
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    for xNodeIdx in range(1,numberOfFluidX2Nodes+1):
        nodeNumber = numberOfFluidX1Nodes+xNodeIdx
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (uInterpolation == CUBIC_HERMITE):
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (problemType == FLUID):
        # Set no slip around the solid
        # Left and right solid edge nodes
        for yNodeIdx in range(2,numberOfSolidYNodes):
            nodeNumber1 = (yNodeIdx-1)*numberOfFluidXNodes2+numberOfFluidX1Nodes
            nodeNumber2 = nodeNumber1+1
            nodeDomain1 = fluidDecomposition.NodeDomainGet(nodeNumber1,1)
            nodeDomain2 = fluidDecomposition.NodeDomainGet(nodeNumber2,1)
            if (nodeDomain1 == computationalNodeNumber):
                fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,1,
                                              iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,2,
                                              iron.BoundaryConditionsTypes.FIXED,0.0)
                if (debugLevel > 2):
                    print('      Node        %d:' % (nodeNumber1))
                if (uInterpolation == CUBIC_HERMITE):    
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
            if (nodeDomain2 == computationalNodeNumber):
                fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,1,
                                              iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,2,
                                              iron.BoundaryConditionsTypes.FIXED,0.0)
                if (debugLevel > 2):
                    print('      Node        %d:' % (nodeNumber2))
                if (uInterpolation == CUBIC_HERMITE):    
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
        # Top solid edge nodes
        for xNodeIdx in range(1,numberOfSolidXNodes+1):
            nodeNumber = xNodeIdx+(numberOfSolidYNodes-1)*numberOfFluidXNodes2+numberOfFluidX1Nodes-1
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                              iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                              iron.BoundaryConditionsTypes.FIXED,0.0)
                if (debugLevel > 2):
                    print('      Node        %d:' % (nodeNumber))
                if (uInterpolation == CUBIC_HERMITE):    
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
                    fsiBoundaryConditions.AddNode(fluidDependentField,iron.FieldVariableTypes.U,1,
                                                  iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                  iron.BoundaryConditionsTypes.FIXED,0.0)
    # Set slip boundary conditions on the top edge
    if (debugLevel > 2):
        print('    Slip Boundary conditions:')
    for xNodeIdx in range(1,numberOfFluidXNodes1+1):
        nodeNumber = numberOfFluidXNodes2*(numberOfSolidYNodes-1)+ \
                     numberOfFluidXNodes1*(numberOfFluidYNodes-1)+xNodeIdx
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (uInterpolation == CUBIC_HERMITE):
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                              iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2, \
                                              nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debugLevel > 2):
        print('    Reference Fluid Pressure Boundary Condition:')
        nodeNumber = numberOfFluidXNodes2
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,fluidPRef)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         =   %.2f' % (fluidPRef))

if (problemType == FSI):
    # Remove dof's at nodes where solid displacement and zero velocity is set (first n last interface node)
    if (debugLevel > 2):
        print('  Lagrange Boundary Conditions:')
        print('    Fixed Boundary conditions:')
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,1, \
                                  iron.BoundaryConditionsTypes.FIXED,0.0)
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,1,2, \
                                  iron.BoundaryConditionsTypes.FIXED,0.0)
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,numberOfInterfaceNodes,1, \
                                  iron.BoundaryConditionsTypes.FIXED,0.0)
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,numberOfInterfaceNodes,2, \
                                  iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debugLevel > 2):
        print('      Node        %d:' % (1))
        print('      Node        %d:' % (numberOfInterfaceNodes))
    if (uInterpolation == CUBIC_HERMITE):
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,1,1, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,1,1, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,1,1, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,1,2, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,1,2, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,1,2, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,numberOfInterfaceNodes,1, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,numberOfInterfaceNodes,1, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,numberOfInterfaceNodes,1, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,numberOfInterfaceNodes,2, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,numberOfInterfaceNodes,2, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,numberOfInterfaceNodes,2, \
                                      iron.BoundaryConditionsTypes.FIXED,0.0)               
# Finish FSI boundary conditions
fsiSolverEquations.BoundaryConditionsCreateFinish()

if (problemType == FSI):
    # Start the creation of the moving mesh boundary conditions
    movingMeshBoundaryConditions = iron.BoundaryConditions()
    movingMeshSolverEquations.BoundaryConditionsCreateStart(movingMeshBoundaryConditions)
    if (debugLevel > 2):
        print('  Moving Mesh Boundary Conditions:')
        print('    Fixed Wall Boundary conditions:')
    # Bottom edge nodes
    for xNodeIdx in range(1,numberOfFluidXNodes2+1):
        nodeNumber = xNodeIdx
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
    # Side edges nodes
    for yNodeIdx in range(2,numberOfSolidYNodes):
        nodeNumber1 = (yNodeIdx-1)*numberOfFluidXNodes2+1
        nodeNumber2 = yNodeIdx*numberOfFluidXNodes2
        nodeDomain1 = fluidDecomposition.NodeDomainGet(nodeNumber1,1)
        nodeDomain2 = fluidDecomposition.NodeDomainGet(nodeNumber2,1)
        if (nodeDomain1 == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,1,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,2,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber1))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
        if (nodeDomain2 == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,1,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,2,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber2))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
    for yNodeIdx in range(1,numberOfFluidYNodes):
        nodeNumber1 = (yNodeIdx-1)*numberOfFluidXNodes1+1+\
                      numberOfFluidXNodes2*(numberOfSolidYNodes-1)
        nodeNumber2 = yNodeIdx*numberOfFluidXNodes1+\
                      numberOfFluidXNodes2*(numberOfSolidYNodes-1)
        nodeDomain1 = fluidDecomposition.NodeDomainGet(nodeNumber1,1)
        nodeDomain2 = fluidDecomposition.NodeDomainGet(nodeNumber2,1)
        if (nodeDomain1 == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,1,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,2,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber1))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
        if (nodeDomain2 == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,1,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,2,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber2))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
    # Top edge nodes
    for xNodeIdx in range(1,numberOfFluidXNodes1+1):
        nodeNumber = xNodeIdx+numberOfFluidXNodes2*(numberOfSolidYNodes-1)+\
                          (numberOfFluidYNodes-1)*numberOfFluidXNodes1
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                movingMeshBoundaryConditions.SetNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
    if (debugLevel > 2):
        print('    Moving Wall Boundary conditions:')
    # Left and right solid edge nodes
    for yNodeIdx in range(2,numberOfSolidYNodes):
        nodeNumber1 = (yNodeIdx-1)*numberOfFluidXNodes2+numberOfFluidX1Nodes
        nodeNumber2 = nodeNumber1+1
        nodeDomain1 = fluidDecomposition.NodeDomainGet(nodeNumber1,1)
        nodeDomain2 = fluidDecomposition.NodeDomainGet(nodeNumber2,1)
        if (nodeDomain1 == computationalNodeNumber):
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,1,
                                                iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber1,2,
                                                iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber1))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber1,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
        if (nodeDomain2 == computationalNodeNumber):
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,1,
                                                iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber2,2,
                                                iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber2))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber2,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
    # Top solid edge nodes
    for xNodeIdx in range(1,numberOfSolidXNodes+1):
        nodeNumber = xNodeIdx+(numberOfSolidYNodes-1)*numberOfFluidXNodes2+numberOfFluidX1Nodes-1
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,
                                                iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,
                                                iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
            if (debugLevel > 2):
                print('      Node        %d:' % (nodeNumber))
            if (uInterpolation == CUBIC_HERMITE):    
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,1,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                    iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1_S2,nodeNumber,2,
                                                    iron.BoundaryConditionsTypes.MOVED_WALL,0.0)

    # Finish moving mesh boundary conditions
    movingMeshSolverEquations.BoundaryConditionsCreateFinish()

if (progressDiagnostics):
    print('Boundary Conditions ... Done')

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

#quit()

# Solve the problem
print('Solving problem...')
start = time.time()
fsiProblem.Solve()
end = time.time()
elapsed = end - start
print('Calculation Time = %3.4f' %elapsed)
print('Problem solved!')
print('#')

#================================================================================================================================
#  Finish Program
#================================================================================================================================
