"""
This file does:
* simulate the point particles thingy like in chapman12 
* simulate the the model for our flexible cell model
"""

using Distributed
NuProcs = 3
# addprocs(NuProcs)

# load all parameters 
@everywhere begin   
    
    include("heatmap.jl")
    include("../simulationFunctionalities.jl")
    include("../parameters.jl")
    include("../energies.jl") 
    include("../cell_functionalities.jl") 

    tspan = timeInterval
    simPath = joinpath(homedir(), "simulations", simulationName)
    locationsPath = joinpath(simPath, "locations")
    heatMapsPath = joinpath(simPath, "heatmaps")
    gifPath = joinpath(simPath, string(simulationName, ".gif"))
    p = [timeStepSize, D]

end

@time runSimulation_locations()








