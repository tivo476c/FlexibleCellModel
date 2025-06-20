## General setup:

domainL = 0.0125 # 0.5 / 40.0
domain = (-domainL, domainL)                # domain where cells can move: [-5.0, 5.0]^2 
NumberOfCellWallPoints = 20                  # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 2                           # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells
D = 1                                       # diffusitivity constant 
radius = 0.2                               # cell radius 
radius = 0.005                               # cell radius 

## Force scalings: 
areaForceFactor = 1e5
edgeForceFactor = 1e1
interiorAngleForceFactor = 1e0
overlapForceFactor = 1e5
overlapForceTypes = ["bachelorThesis", "billiard", "combination", "radiusBilliard"]
overlapForceType = overlapForceTypes[1]

forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor]

## Simulation time parameters: 
T = 10^(-5) * 100
timeInterval = (0.0, T)
timeStepSize = 10^(-6)

# NumberOfSimulations = 10^4
NumberOfSampleTimes = 101          # must be 2 at least
sampleTimes = [T * k / (NumberOfSampleTimes - 1) for k = 0:NumberOfSampleTimes-1]
# sampleTimesRange = 0:T/(NumberOfSampleTimes-1):T

## Simulation name 
# date = today()
# currentTime = Dates.format(now(), "HH-MM")
currentTime = "13-21"
# simulationName = string("PP-SIM_", date, "_", currentTime)
simulationName = "DF-5-1-0-$overlapForceFactor-dt$timeStepSize-finer"

## Space Discretisation for heatmap 
NumberOfHeatGridPoints = 30
HeatStepSize = 2 * domainL / NumberOfHeatGridPoints
# [-L, -L + HeatStepSize, -L + 2*HeatStepSize, ..., -L + NumberOfHeatGridPoints*HeatStepSize]
HeatGrid = [-domainL + k * HeatStepSize for k in 0:NumberOfHeatGridPoints]
