#! /usr/bin/env pvpython

from math import *
from paraview.simple import *

# Vyhodnoceni probehne v rovinach 20mm pred a 20mm za lopatkou
traversOffset = 20e-3


foamfoam = OpenFOAMReader(FileName='foam.foam')
foamfoam.CellArrays = ['p', 'U']
finalTime = foamfoam.TimestepValues[-1]

# Urci pozici nabezne a odtokove hrany a nastavi souradnice rezu 1 a 2
foamfoam.MeshRegions = ['PROF']
foamfoam.UpdatePipeline(time=finalTime)
boundingBox = foamfoam.GetDataInformation().GetBounds()
x1 = boundingBox[0] - traversOffset
x2 = boundingBox[1] + traversOffset


# Prida k datum celkovy tlak
data = Calculator(Input=foamfoam)
data.ResultArrayName = 'pTot'
data.Function = 'p+0.5*mag(U)^2'
data.AttributeType = 'Cell Data'


# Pocita vazeny prumer 
def WeightedAverage(Input, weight, variable):
    calc1 = Calculator(Input=Input)
    calc1.ResultArrayName = 'weight'
    calc1.Function = weight
    calc1.AttributeType = 'Cell Data'

    calc2 = Calculator(Input=calc1)
    calc2.ResultArrayName = 'weightedVar'
    calc2.Function = weight + '*' + variable
    calc2.AttributeType = 'Cell Data'
    
    integrals = IntegrateVariables(Input=calc2)
    results = servermanager.Fetch(integrals)
    w = results.GetCellData().GetArray('weight').GetValue(0)
    var = results.GetCellData().GetArray('weightedVar').GetValue(0)
    return var/w 
    
def FlowRate(Input, velocity):
	calc1 = Calculator(Input=Input)
	calc1.ResultArrayName = 'velocity'
	calc1.Function = velocity
	calc1.AttributeType = 'Cell Data'
	integrals=IntegrateVariables(Input=calc1)
	results = servermanager.Fetch(integrals)
	
	return results.GetCellData().GetArray('velocity').GetValue(0)

foamfoam.MeshRegions = ['internalMesh']

stage1 = Slice(Input=data)
stage1.SliceType = 'Plane'
stage1.SliceType.Origin = [x1, 0.0, 0.0]
stage1.SliceType.Normal = [1.0, 0.0, 0.0]

stage2 = Slice(Input=data)
stage2.SliceType = 'Plane'
stage2.SliceType.Origin = [x2, 0.0, 0.0]
stage2.SliceType.Normal = [1.0, 0.0, 0.0]

p1    = WeightedAverage(stage1, 'U_X', 'p')
pTot1 = WeightedAverage(stage1, 'U_X', 'pTot')

p2    = WeightedAverage(stage2, 'U_X', 'p')
flow2 = FlowRate(stage2, 'U_X')
pTot2 = WeightedAverage(stage2, 'U_X', 'pTot')

foamfoam.MeshRegions = ['PRESSURE-OUTLET']
pOut = WeightedAverage(data, 'U_X', 'p')
pTotOut = WeightedAverage(data, 'U_X', 'pTot')
flowOut = FlowRate(data, 'U_X')

foamfoam.MeshRegions = ['PRESSURE-INLET']
pIn = WeightedAverage(data, 'U_X', 'p')
pTotIn = WeightedAverage(data, 'U_X', 'pTot')

print("Pressure loss (1-2)        [%]:", (pTot1 - pTot2)/(pTot1 - p1)*100)
print("Pressure loss (In-Out)     [%]:", (pTotIn - pTotOut)/(pTotIn - pIn)*100)
print("Compression ratio (1-2)    [1]:", pTot1/pTot2)
print("Compression ratio (In-Out) [1]:", pTotIn/pTotOut)
print("Volume flow (2)         [ms-2]:", flow2)
print("Volume flow (Out)       [ms-2]:", flowOut)


