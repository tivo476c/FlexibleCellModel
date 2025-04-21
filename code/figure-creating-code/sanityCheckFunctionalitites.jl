include("../energies.jl")
include("../simulationFunctionalities.jl")
include("../parameters.jl")
include("heatmap.jl")

using Distributions, Distributed


function giveCentreNormalDistrInDomain(radius; mean=0.0, deviation=0.9^2)
    """
        finds a point (x,y) in [-5+radius, 5-radius]^2 that gets initialized according to the Normal distribution N(mean, deviation)
    """
    dist = MvNormal([mean, mean], deviation * I(2))
    while true
        numX, numY = rand(dist)
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
        if (norm(newCell - c) < 2 * radius)
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
        while (!isFeasible(newCentre, savedCentres, radius))
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

    return vcat(xCoords, yCoords)
end

function InitializePointParticles(radius)
    """
    It does not matter which normally distributed coordinate gets used for x and which for y coordinate. 
        -> just compute NumberOfCells x coordinates, than all y coordinates 
    """
    savedCentres = []
    for i = 1:M
        newCentre = giveCentreNormalDistrInDomain(radius)
        while (!isFeasible(newCentre, savedCentres, radius))
            newCentre = giveCentreNormalDistrInDomain(radius)
        end
        push!(savedCentres, newCentre)
    end
    xCoords = Float64[]
    yCoords = Float64[]
    for centre in savedCentres
        push!(xCoords, centre[1])
        push!(yCoords, centre[2])
    end

    return vcat(xCoords, yCoords)
end

function doAPointParticleSimulationRun(simRun)
    println("simrun ", simRun)
    u0 = InitializePointParticles(radius)
    prob_pointParticles = SDEProblem(energies, brownian, u0, tspan, p, noise=WienerProcess(0.0, 0.0))
    sol = solve(prob_pointParticles, EM(), dt=timeStepSize, saveat=sampleTimes)
    createLocationFile(sol, simRun, locationsPath)
end

function doAHSCMSimulationRun(simRun)
    println("simrun ", simRun)
    u0 = initializeChapman(radius)
    prob_pointParticles = SDEProblem(energies, brownian, u0, tspan, p, noise=WienerProcess(0.0, 0.0))
    sol = solve(prob_pointParticles, EM(), dt=timeStepSize, saveat=sampleTimes)
    createLocationFile(sol, simRun, locationsPath)
end