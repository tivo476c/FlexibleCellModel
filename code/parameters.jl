## General setup:

domainL = 0.0125 # 0.5 / 40.0
domain = (-domainL, domainL)                # domain where cells can move: [-5.0, 5.0]^2 
NumberOfCellWallPoints = 6                  # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 2                           # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells
D = 1                                       # diffusitivity constant 
radius = 0.005                              # cell radius 

## Force scalings: 
hardness = 1                                       # tells how hard the cells are hardness =1 for hard cells, hardness =0 for soft cells, or something in between 

areaForceFactor = (1-hardness)*1e6 + hardness*1e10
edgeForceFactor = (1-hardness)*1e3 + hardness*1e5
interiorAngleForceFactor = (1-hardness)*1e0 + hardness*1e2
overlapForceFactor = 1e5 
overlapForceTypes = ["bachelorThesis", "radiusBilliard", "combination"]
overlapForceType = overlapForceTypes[3]

forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor]

## Simulation time parameters: 
T = 10^(-6) * 20
timeInterval = (0.0, T)
timeStepSize = 10^(-6)

# NumberOfSimulations = 10^4
NumberOfSampleTimes = 21          # must be 2 at least
sampleTimes = [T * k / (NumberOfSampleTimes - 1) for k = 0:NumberOfSampleTimes-1]
# sampleTimesRange = 0:T/(NumberOfSampleTimes-1):T

## Simulation name 
# date = today()
# currentTime = Dates.format(now(), "HH-MM")
currentTime = "13-21"
# simulationName = string("PP-SIM_", date, "_", currentTime)
simulationName = "drift-$(floor(Int64, log10(areaForceFactor)))-$(floor(Int64,log10(edgeForceFactor)))-$(floor(Int64, log10(interiorAngleForceFactor)))-$(floor(Int64,log10(overlapForceFactor)))"

## Space Discretisation for heatmap 
NumberOfHeatGridPoints = 30
HeatStepSize = 2 * domainL / NumberOfHeatGridPoints
# [-L, -L + HeatStepSize, -L + 2*HeatStepSize, ..., -L + NumberOfHeatGridPoints*HeatStepSize]
HeatGrid = [-domainL + k * HeatStepSize for k in 0:NumberOfHeatGridPoints]
