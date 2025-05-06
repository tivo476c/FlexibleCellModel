using Dates

## General setup:

domainL = 5.0
domain = (-domainL, domainL)                # domain where cells can move: [-5.0, 5.0]^2 
NumberOfCellWallPoints = 0                  # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 20                          # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells
D = 100                          # diffusitivity constant 
radius = 0.05           # cell radius 

## Force scalings: 
areaForceFactor = 0.0
edgeForceFactor = 0.0
interiorAngleForceFactor = 0.0
overlapForceFactor = 0.0
overlapForceTypes = ["bachelorThesis", "billiard", "combination", "radiusBilliard"]
overlapForceType = overlapForceTypes[4]
boundaryPushForceFactor = 0.0

forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor, boundaryPushForceFactor]


## Simulation time parameters: 
T = 10.0
timeInterval = (0.0, T)
timeStepSize = 10^(-4)

## sampleTimes = [k/10.0 for k = 0:100]
NumberOfSimulations = 10
NumberOfSampleTimes = 10
sampleTimes = [T*k/(NumberOfSampleTimes-1) for k = 0:NumberOfSampleTimes-1]
NumberOfSampleTimes = length(sampleTimes)

## Simulation name 
date = today()
# date = "2025-04-21"
currentTime = Dates.format(now(), "HH-MM")
# currentTime = "09-36"
simulationName = string("HSCM-SIM_", date, "_", currentTime)

## Space Discretisation for heatmap 
NumberOfHeatGridPoints = 50 
HeatStepSize = 2*domainL / NumberOfHeatGridPoints
# [-L, -L + HeatStepSize, -L + 2*HeatStepSize, ..., -L + NumberOfHeatGridPoints*HeatStepSize]
HeatGrid = [-domainL + k * HeatStepSize for k in 0:NumberOfHeatGridPoints]
