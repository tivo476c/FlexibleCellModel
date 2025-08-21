using Pkg
"""
This file is the entry point of my FlexibleCellModel code. 
It gets executed, whenever parameters.jl gets executed.
It uses its parameter configuration to start the according simulation. 
"""
# Pkg.develop(path="C:/Users/voglt/Desktop/FlexibleCellModel")

println("Started entryPoint.jl")
println("Loading all packages")


###  OPTION 1: JUST RUN A METHOD WITHOUT PARALLELIZING
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

runShow_overlap()
# @time runSimulation_locations()

###  OPTION 2: PARALLELIZED 
# if wanna use parallelized run, inclusions happen in startParallelizedRun.jl: 

# using Pkg 
# Pkg.activate(".")
# Pkg.instantiate()
# include("code/startParallelizedRun.jl")
