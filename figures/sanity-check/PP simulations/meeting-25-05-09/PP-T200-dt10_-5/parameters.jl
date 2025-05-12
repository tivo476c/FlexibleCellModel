using Dates

## General setup:

domainL = 0.5
domain = (-domainL, domainL)                # domain where cells can move: [-5.0, 5.0]^2 
NumberOfCellWallPoints = 0                  # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 400                          # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells
D = 1                                       # diffusitivity constant 
radius = 0.005                               # cell radius 

## Force scalings: 
areaForceFactor = 0.0
edgeForceFactor = 0.0
interiorAngleForceFactor = 0.0
overlapForceFactor = 0.0
overlapForceTypes = ["bachelorThesis", "billiard", "combination", "radiusBilliard"]
overlapForceType = overlapForceTypes[4]
boundaryPushForceFactor = 0.01

forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor, boundaryPushForceFactor]


## Simulation time parameters: 
T = 200.0
timeInterval = (0.0, T)
timeStepSize = 10^(-5)

## sampleTimes = [k/10.0 for k = 0:100]
NumberOfSimulations = 100
NumberOfSampleTimes = 11
sampleTimes = [T*k/(NumberOfSampleTimes-1) for k = 0:NumberOfSampleTimes-1]

## Simulation name 
# date = today()
# date = "2025-05-08"
# # currentTime = Dates.format(now(), "HH-MM")
# currentTime = "12-00"
# simulationName = string("PP-SIM_", date, "_", currentTime)
simulationName = "PP-T200-dt10_-5"

## Space Discretisation for heatmap 
NumberOfHeatGridPoints = 50 
HeatStepSize = 2*domainL / NumberOfHeatGridPoints
# [-L, -L + HeatStepSize, -L + 2*HeatStepSize, ..., -L + NumberOfHeatGridPoints*HeatStepSize]
HeatGrid = [-domainL + k * HeatStepSize for k in 0:NumberOfHeatGridPoints]
