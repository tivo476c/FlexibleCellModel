## General setup:

domainL = 2.0                               # 0.5 / 40.0
domain = (-domainL, domainL)                # domain where cells can move: [-5.0, 5.0]^2 
NumberOfCellWallPoints = 4                  # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 1                        # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells
D = 1                                       # diffusitivity constant 
radius = 1.0                              # cell radius 

## Force scalings: 
hardness = 0.5                                       # tells how hard the cells are hardness =1 for hard cells, hardness =0 for soft cells, or something in between 

areaForceFactor = (1 - hardness) * 1e6 + hardness * 1e10
edgeForceFactor = (1 - hardness) * 1e1 + hardness * 1e5
interiorAngleForceFactor = (1 - hardness) * 1e0 + hardness * 1e2
overlapForceFactor = 2e4
# forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor]
forceScalings = [0, edgeForceFactor, 0, 0]

overlapForceTypes = ["bachelorThesis", "radiusBilliard", "combination"]
overlapForceType = overlapForceTypes[2]

## Simulation time parameters: 
timeStepSize = 10^(-5)
T = timeStepSize*5
timeInterval = (0.0, T)

NumberOfSimulations = 1
NumberOfSampleTimes = 6          # must be 2 at least
sampleTimes = [T * k / (NumberOfSampleTimes - 1) for k = 0:NumberOfSampleTimes-1]
# sampleTimesRange = 0:T/(NumberOfSampleTimes-1):T

## Simulation name 
# date = today()
# currentTime = Dates.format(now(), "HH-MM")
# currentTime = "13-21"
simulationName = "testEdge2-decreaseEdges-N_v$NumberOfCellWallPoints"

## Space Discretisation for heatmap 
NumberOfHeatGridPoints = 50
HeatStepSize = 2 * domainL / NumberOfHeatGridPoints
# [-L, -L + HeatStepSize, -L + 2*HeatStepSize, ..., -L + NumberOfHeatGridPoints*HeatStepSize]
HeatGrid = [-domainL + k * HeatStepSize for k in 0:NumberOfHeatGridPoints]
