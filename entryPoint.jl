using Pkg

"""
This file is the entry point of my FlexibleCellModel code. 
It gets executed, whenever parameters.jl gets executed.
It uses its parameter configuration to start the according simulation. 
"""

println("Started entryPoint.jl")
println("Loading all packages")


#######  OPTION 1: JUST RUN A METHOD WITHOUT PARALLELIZING
include("code/parameters.jl")
include("code/cell_functionalities.jl")
include("code/computeOverlap.jl")
include("code/energies.jl")
include("code/figure-creating-code/heatmap.jl")
include("code/simulationFunctionalities.jl")

tspan = timeInterval
simPath = joinpath(homedir(), "simulations", simulationName)
locationsPath = joinpath(simPath, "locations")
heatMapsPath = joinpath(simPath, "heatmaps")
gifPath = joinpath(simPath, string(simulationName, ".gif"))
energyDiaPath = joinpath(simPath, "energies-$simulationName.png")
p = [timeStepSize, D]


# runSimulation_locations()                                   # for producing location data and heatmaps

#######  OPTION 2: PARALLELIZED 
# if wanna use parallelized run, inclusions happen in startParallelizedRun.jl: 

# using Pkg 
# Pkg.activate(".")
# Pkg.instantiate()
# include("code/startParallelizedRun.jl")


#######  OPTION 3: PLAYGROUND 
createCrossSectionAllPlots(["AAAsoftSim-bachelorOverlap", "AAAmidSim3-bachelorOverlap", "AAAtest3-hard-DF-FIRSTWORKING"])
# createCrossSectionAllPlots("AAAmidSim3-bachelorOverlap")
# createCrossSectionAllPlots("AAAtest3-hard-DF-FIRSTWORKING")
# doAsphericityCheck(simList=["AAAsoftSim-bachelorOverlap", "AAAmidSim3-bachelorOverlap", "AAAtest3-hard-DF-FIRSTWORKING"])



######################################## now doing the plotting 

# saveAsphericityData( simList=["AAAsoftSim-bachelorOverlap", "AAAmidSim3-bachelorOverlap", "AAAtest3-hard-DF-FIRSTWORKING"], Nsims=100)

# doAsphericityCheck(simList=["AAAsoftSim-bachelorOverlap", "AAAmidSim3-bachelorOverlap", "AAAtest3-hard-DF-FIRSTWORKING"])
# createAsphericityPlot()