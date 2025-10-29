"""
This file does:
* simulate the point particles thingy like in chapman12 
* simulate the the model for our flexible cell model
"""

println("Started startParallelizedRun.jl")


using Distributed

NuProcs = 60
if nprocs() == 1
    addprocs(NuProcs)
end
println("Using $(nprocs()) processors")
println("Loading parameters")
# load all parameters 
@everywhere begin

    # including cell functionalities with decreasing grade of fundamentality 
    include("parameters.jl")
    println("parameters.jl done")
    include("cell_functionalities.jl")
    println("cell_functionalities.jl done")
    include("computeOverlap.jl")
    println("computeOverlap.jl done")
    include("energies.jl")
    println("energies.jl done")
    include("figure-creating-code/heatmap.jl")
    println("figure-creating-code/heatmap.jl done")
    include("simulationFunctionalities.jl")  
    println("simulationFunctionalities.jl done")

    tspan = timeInterval
    simPath = joinpath(homedir(), "simulations", simulationName)
    locationsPath = joinpath(simPath, "locations")
    heatMapsPath = joinpath(simPath, "heatmaps")
    gifPath = joinpath(simPath, string(simulationName, ".gif"))
    p = [timeStepSize, D, 0, 0, 0]

end

# @time runSimulation_locations()
println("Starting simulation: $simulationName")

# run wanted function 
@time runSimulation_locations()

# delete workers 
rmprocs(workers())







