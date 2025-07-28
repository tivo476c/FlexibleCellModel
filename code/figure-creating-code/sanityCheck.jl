"""
This file does:
* simulate the point particles thingy like in chapman12 
* simulate the the model for our flexible cell model
"""

println("Started sanityCheck.jl \nHere a paralized run can be implemented.")


using Distributed

NuProcs = 3
if nprocs() == 1
    addprocs(NuProcs)
end
println("Using $(nprocs()) processors")
println("Loading parameters")
# load all parameters 
@everywhere begin   
    
    # including cell functionalities with decreasing grade of fundamentality 
    include("../parameters.jl")
    include("../cell_functionalities.jl") 
    include("../computeOverlap.jl") 
    
    include("../energies.jl") 
    
    include("heatmap.jl")
    
    include("../simulationFunctionalities.jl")



    tspan = timeInterval
    simPath = joinpath(homedir(), "simulations", simulationName)
    locationsPath = joinpath(simPath, "locations")
    heatMapsPath = joinpath(simPath, "heatmaps")
    gifPath = joinpath(simPath, string(simulationName, ".gif"))
    gifPath = joinpath(simPath, string(simulationName, ".gif"))
    p = [timeStepSize, D]

end

# @time runSimulation_locations()
println("Starting simulation: $simulationName")

@time runShow_overlap()








