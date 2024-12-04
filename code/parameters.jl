# General setup:

domain = (-5.0, 5.0)                # domain where cells can move: [-5.0, 5.0]^2 
N = 40                              # number of wall points per cell 
M = 16                               # number of cells 
D = 5.0                             # diffusitivity constant 


# Force scalings: 
# [70,2,1.5,110]
areaForceFactor = 70.0 
edgeForceFactor = 2.0
interiorAngleForceFactor = 0.0

overlapForceFactor = 110.0
overlapForceTypes = ["bachelorThesis", "billiard", "combination"]
overlapForceType = overlapForceTypes[3]
scalingBachelor = 0.5

boundaryPushForceFactor = 110.0
forceScalings = [areaForceFactor, edgeForceFactor, interiorAngleForceFactor, overlapForceFactor, boundaryPushForceFactor]


# Simulation parameters: 

timeInterval = (0.0, 10.0)
timeStepSize = 2.0e-8 
lengthOfSolVec = Int((timeInterval[2] - timeInterval[1]) / timeStepSize)    # maybe length is one more

sampleTimeStepSize = floor(lengthOfSolVec/10.0)                             # sample times for creating heatmaps
#sampleTimes = collect(0:sampleTimeStepSize:lengthOfSolVec) 
sampleTimes = collect(0:5000:lengthOfSolVec)                    
sampleTimes[1] = 1 

saveAtTimes = 0.5

NumberOfSimulations = 10 
NumberOfSampleTimes = length(sampleTimes) 
simulationName = string("testrun", NumberOfSimulations, "-T",timeInterval[2]) 

# Space Discretisation

grid = [-5.0 + k * 0.25 for k in 0:40]
discreteSpaceArray = zeros(Int, 40, 40)
function giveCoordinates(x::Int, y::Int) 
    x_min = -5.0 + (x-1) * 0.25
    x_max = x_min + 0.25 
    y_max =  5.0 - (y-1) * 0.25
    y_min = y_max - 0.25 

    return (x_min, x_max), (y_min, y_max)
end 





