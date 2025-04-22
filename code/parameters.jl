using Dates

# General setup:

domain = (-5.0, 5.0)                # domain where cells can move: [-5.0, 5.0]^2 
NumberOfCellWallPoints = 3               # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 400                         # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells
D = 100                          # diffusitivity constant 
radius = 0.05           # cell radius 

# Force scalings: 
# [70,2,1.5,110]
areaForceFactor = 0.0
edgeForceFactor = 0.0
interiorAngleForceFactor = 0.0
overlapForceFactor = 1.0
overlapForceTypes = ["bachelorThesis", "billiard", "combination", "radiusBilliard"]
overlapForceType = overlapForceTypes[4]
boundaryPushForceFactor = 1.0

forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor, boundaryPushForceFactor]


# Simulation parameters: 

timeInterval = (0.0, 10.0)
timeStepSize = 10^(-4)
lengthOfSolVec = floor((timeInterval[2] - timeInterval[1]) / timeStepSize)   # maybe length is one more
sampleTimeStepSize = floor(lengthOfSolVec / 10.0)

#sampleTimes = collect(0:sampleTimeStepSize:lengthOfSolVec) 
#sampleTimes = collect(0:100:lengthOfSolVec) 
#sampleTimes = Int.(sampleTimes)              # Convert each element to Int
#sampleTimes[1] = 1 

# sampleTimes = [k/10.0 for k = 0:100]
sampleTimes = [k/10.0 for k = 0:100]

NumberOfSimulations = 39  # TODO: change to 10^4 or something like this (test how many!!!)
NumberOfSampleTimes = length(sampleTimes)


date = today()
# date = "2025-04-21"
# currentTime = Dates.format(now(), "HH-MM")
currentTime = "09-36"

simulationName = string("HSCM-SIM_", date, "_", currentTime)

# Space Discretisation

grid = [-5.0 + k * 0.25 for k in 0:40]
discreteSpaceArray = zeros(Int, 40, 40)
function giveCoordinates(x::Int, y::Int)
    x_min = -5.0 + (x - 1) * 0.25
    x_max = x_min + 0.25
    y_max = 5.0 - (y - 1) * 0.25
    y_min = y_max - 0.25

    return (x_min, x_max), (y_min, y_max)
end





