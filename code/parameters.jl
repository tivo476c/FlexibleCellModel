## General setup:

domainL = 0.5                               # 0.5 / 40.0
domain = (-domainL, domainL)                # domain where cells can move: [-5.0, 5.0]^2 
NumberOfCellWallPoints = 6                  # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 400                           # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells
D = 1                                       # diffusitivity constant 
radius = 0.005                              # cell radius 

## Force scalings: 
hardness = 0.5                                       # tells how hard the cells are hardness =1 for hard cells, hardness =0 for soft cells, or something in between 

areaForceFactor = 4e8
edgeForceFactor = 3e4
interiorAngleForceFactor = 1e-1
overlapForceFactor = 6e4
# areaForceFactor = (1 - hardness) * 1e6 + hardness * 1e10
# edgeForceFactor = (1 - hardness) * 1e1 + hardness * 1e5
# interiorAngleForceFactor = (1 - hardness) * 1e0 + hardness * 1e2
# overlapForceFactor = 2e4
overlapForceTypes = ["bachelorThesis", "radiusBilliard", "combination"]
overlapForceType = overlapForceTypes[3]

forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor]

## Simulation time parameters: 
timeStepSize = 10^(-5)
T = 0.05
# T = 20*timeStepSize
timeInterval = (0.0, T)

NumberOfSimulations = Int(1.2 * 10^4)
NumberOfSampleTimes = 6        # must be 2 at least
sampleTimes = [T * k / (NumberOfSampleTimes - 1) for k = 0:NumberOfSampleTimes-1]
# sampleTimesRange = 0:T/(NumberOfSampleTimes-1):T

## Simulation name 
# simulationName = "drift-$(floor(Int64, log10(areaForceFactor)))-$(floor(Int64,log10(edgeForceFactor)))-$(floor(Int64, log10(interiorAngleForceFactor)))-$(floor(Int64,log10(overlapForceFactor)))"
simulationName = "midSim1-masterOverlap"


## Space Discretisation for heatmap 
NumberOfHeatGridPoints = 30
HeatStepSize = 2 * domainL / NumberOfHeatGridPoints
# [-L, -L + HeatStepSize, -L + 2*HeatStepSize, ..., -L + NumberOfHeatGridPoints*HeatStepSize]
HeatGrid = [-domainL + k * HeatStepSize for k in 0:NumberOfHeatGridPoints]
