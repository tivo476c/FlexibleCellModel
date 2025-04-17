# General setup:

domain = (-5.0, 5.0)                # domain where cells can move: [-5.0, 5.0]^2 
NumberOfCellWallPoints = 0               # number of wall points per cell [OLD NAME: "N"]
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
overlapForceFactor = 0.0
overlapForceTypes = ["bachelorThesis", "billiard", "combination", "radiusBilliard"]
overlapForceType = overlapForceTypes[4]
boundaryPushForceFactor = 10.0

forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor, boundaryPushForceFactor]


# Simulation parameters: 

timeInterval = (0.0, 0.05)
timeStepSize = 10^(-5)
lengthOfSolVec = floor((timeInterval[2] - timeInterval[1]) / timeStepSize)   # maybe length is one more
sampleTimeStepSize = floor(lengthOfSolVec / 10.0)

#sampleTimes = collect(0:sampleTimeStepSize:lengthOfSolVec) 
#sampleTimes = collect(0:100:lengthOfSolVec) 
#sampleTimes = Int.(sampleTimes)              # Convert each element to Int
#sampleTimes[1] = 1 

# sampleTimes = [1,1000,3000,5000]
sampleTimes = [0.0, 0.01, 0.02, 0.05]

NumberOfSimulations = 10^3  # TODO: change to 10^4 or something like this (test how many!!!)
NumberOfSampleTimes = length(sampleTimes)
#simulationName = string("run1-NoSim", NumberOfSimulations, "-T",timeInterval[2]) 
simulationName = string("ppmodel-16-04-25-lappi1")

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





