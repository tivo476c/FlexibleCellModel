include("../parameters.jl")
include("../energies.jl")
include("../simulationFunctionalities.jl")
include("heatmap.jl")

using Distributions

function giveCentreNormalDistrInDomain(radius, mean=0.0, std_dev=1.0)
    dist = Normal(mean, std_dev)
    while true
        numX, numY = rand(dist), rand(dist)
        if -5 + radius <= numX <= 5 - radius && -5 + radius <= numY <= 5 - radius
            return [numX, numY]
        end
    end
end

function isFeasible(newCell, oldCells, radius) 
    for c in oldCells
        if(norm(newCell - c) < 2*radius)
            return false 
        end 
    end 
    return true 
end 

begin
    include("../parameters.jl")
    include("../energies.jl")
    include("heatmap.jl")
    
    function initializeChapman(radius)
        """
        Initializes a starting vector u0 for the solver to work with it. 
        Cell by cell gets inserted at a feasible location in the domain. 
        A position of newCell is feasible if distance(newCell, oldCell) > 2*radius for all oldCells.
        This function first just computes a list that holds all cell centres. From this list, all cell wallpoints get computed and outputted in u0. 
        """
        savedCentres = []
        for i = 1:M 
            newCentre = giveCentreNormalDistrInDomain(radius)
            while(! isFeasible(newCentre, savedCentres, radius))
                newCentre = giveCentreNormalDistrInDomain(radius)
            end 
            push!(savedCentres, newCentre)
        end 
        xCoords = Float64[]
        yCoords = Float64[]
        for centre in savedCentres
            discreteCell = cellToDiscreteCell(circleCell(centre, radius), N)

            xCoords = vcat(xCoords, discreteCell.x)
            yCoords = vcat(yCoords, discreteCell.y)
        end 
            
        return vcat(xCoords,yCoords)
    end 
    
    simPath = joinpath(homedir(), "simulations", simulationName)
    locationsPath = joinpath(simPath, "locations")
    heatMapsPath = joinpath(simPath, "heatmaps")
    
    C = circleCell([0.0,0.0], radius)
    cDF = cellToDiscreteCell(C, N) 
    A1 = ones(M) * areaPolygon(cDF.x, cDF.y) # ∈ R^N
    E1 = ones(N*M)              # ∈ (R^N)^M
    I1 = ones(N*M)              # ∈ (R^N)^M
    e = computeEdgeLengths(cDF)
    ia = computeInteriorAngles(cDF)
    
    # initial condition 
    
    for i = 1:M
        E1[(i-1)*N+1 : i*N] = e
        I1[(i-1)*N+1 : i*N] = ia
    end 
    
    u0 = initializeChapman(radius)
    
    tspan = timeInterval
    p = [timeStepSize, D]
    prob_cell1 = SDEProblem(energies, brownian, u0, tspan, p, noise=WienerProcess(0., 0.))     
    @time sol = solve(prob_cell1, EM(), dt=timeStepSize, saveat = collect(0:15*timeStepSize:timeInterval[2]))    
    
    # HEATMAP stuff begins 
    
    mkpath(simPath) 
    cp(joinpath(homedir(), "OneDrive", "Desktop", "Master-Thesis", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    mkpath(locationsPath)
    mkpath(heatMapsPath)
    
    # save one simulation as give 
    gifPath = joinpath(simPath, string(simulationName, ".gif"))
    createSimGif(gifPath, sol)
    @distributed for i = 1:NumberOfSimulations
        @time sol = solve(prob_cell1, EM(), dt=timeStepSize, saveat=sampleTimes)
        createLocationFile(sol, i, locationsPath)
    end 

    matrices = makeMatrices()

    createHeatmaps(matrices)

end

