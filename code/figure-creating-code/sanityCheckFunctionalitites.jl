include("../energies.jl")
include("../simulationFunctionalities.jl")
include("../parameters_pointParticles.jl")
include("heatmap.jl")

using Distributions, Distributed


function giveCentreNormalDistrInDomain(radius, mean=0.0, deviation=1.0)
    """
        finds a point (x,y) in [-5+radius, 5-radius]^2 that gets initialized according to the Normal distribution N(mean, deviation)
    """
    dist = Normal(mean, deviation)
    while true
        numX, numY = rand(dist), rand(dist)
        if -5 + radius <= numX <= 5 - radius && -5 + radius <= numY <= 5 - radius
            return [numX, numY]
        end
    end
end

function isFeasible(newCell, oldCells, radius) 
    """
    Returns true, iff |newCell - oldCell| > 2*radius for all oldCells     
    
    args:
            * newCell: midpoint of the new cell 
            * oldCells: list of midpoints of all already initialized cells 
            * radius: cell radius of all cells (that get initialized as circles in this model) 
    """

    for c in oldCells
        if(norm(newCell - c) < 2*radius)
            return false 
        end 
    end 
    return true 
end 

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

function InitializePointParticles()
    """
    It does not matter which normally distributed coordinate gets used for x and which for y coordinate. 
        -> just compute NumberOfCells x coordinates, than all y coordinates 
    """
    res = zeros(Float64, 2*NumberOfCells)
    dist = Normal(0.0, 10.0)
    for i = 1:2*NumberOfCells
        newCoord = rand(dist)

        while  newCoord < -5 || newCoord > 5 
            newCoord = rand(dist)
        end        
        res[i] = newCoord
    end
    return res 
end 

function doAPointParticleSimulationRun(simRun)
    println("simrun ", simRun)
    u0 = InitializePointParticles()
    prob_pointParticles = SDEProblem(energies, brownian, u0, tspan, p, noise=WienerProcess(0., 0.))     
    sol = solve(prob_pointParticles, EM(), dt=timeStepSize, saveat = sampleTimes)  
    createLocationFile(sol, simRun, locationsPath)
end 