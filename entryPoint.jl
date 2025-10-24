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
# createCrossSectionAllPlots("AAAsoftSim-bachelorOverlap")
# createCrossSectionAllPlots("AAAmidSim3-bachelorOverlap")
# createCrossSectionAllPlots("AAAtest3-hard-DF-FIRSTWORKING")
# createCrossSectionAllPlots(["AAAsoftSim-bachelorOverlap", "AAAmidSim3-bachelorOverlap", "AAAtest3-hard-DF-FIRSTWORKING"])

# doAsphericityCheck()
#create path

aspPath = joinpath(simPath, "asphericity")
if !isdir(aspPath)
    mkdir(aspPath)
end 

A_d, E_d, I_d = computeDesiredStates_circleCells()
p = timeStepSize, D, A_d, E_d, I_d 
u0 = initializeCells(radius)

cellProblem = SDEProblem(energies!, brownian_DF!, u0, timeInterval, p, noise_rate_prototype=zeros(2 * M * N, 2 * M))
@time sol = solve(cellProblem,
    EM(),
    dt=timeStepSize,
)
extractedSol = extractSolution(sol)
# extractedSol.t[timeStep in {1,...,6}] ... time steps 
# extractedSol.u[timeStep in {1,...,6}] ... cell configs to that time 

# collect asphericity data 
AllAsphericities = []
for timesteps = 1:length(extractedSol.t)
    allCells = solutionToCells(extractedSol.u[timesteps])

    Asphericities_1time = []
    for cell in allCells    
        # cell.x for all x coords and cell.y for all y coords 
        cell_area = areaPolygon(cell.x, cell.y)
        cell_perimeter = sum(computeEdgeLengths(cell))
        asphericity = 4*pi*cell_area/(cell_perimeter^2)
        push!(Asphericities_1time, asphericity)
    end 
    push!(AllAsphericities, Asphericities_1time)
end
# AllAsphericities[timeStep in {1,...,6}][cell in {1,...,400}]

begin
    xLowerBound = 0.50
    xUpperBound = 0.95
    Nintervals = 10
    xStepSize = (xUpperBound - xLowerBound) / Nintervals
    factor = xStepSize^(-1)
    aspValues = xLowerBound:xStepSize:xUpperBound
    for timestep = 1:length(extractedSol.t)
        xScale  = ["$(aspValues[i]) - $(aspValues[i+1])" for i=1:length(aspValues)-1]  
        yCounts = zeros(length(aspValues)-1)
        for cellAsp in AllAsphericities[timestep]
            yIntervalIndex = ceil(Int, factor*(cellAsp - xLowerBound)) 
            yCounts[yIntervalIndex] += 1
        end 
        aspPlot = bar(xScale, yCounts, 
        label=false,
        title="Aspherictiy",
        xlabel="asphericity",
        ylabel="Number of cells",
        xrotation = 45,
        dpi=500,
        )
        barname = "asphericity-chart-time$timestep.png"
        savefig(joinpath(aspPath, barname))
    end 
end