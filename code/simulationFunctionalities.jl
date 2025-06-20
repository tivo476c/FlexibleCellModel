include("energies.jl")
include("cell_functionalities.jl")
# include("parameters.jl")
include("figure-creating-code/heatmap.jl")

using Distributions, Distributed, Printf, Sockets

function giveCentreNormalDistrInDomain(radius; mean=0.0, deviation=0.09^2)
    """
        finds a point (x,y) in [-5+radius, 5-radius]^2 that gets initialized according to the Normal distribution N(mean, deviation)
    """
    dist = MvNormal([mean, mean], deviation * I(2))
    for _ = 1:1000
        numX, numY = rand(dist)
        if -domainL + radius <= numX <= domainL - radius && -domainL + radius <= numY <= domainL - radius
            return [numX, numY]
        end
    end
    println("ERROR: could not find a point in the domain in giveCentreNormalDistrInDomain in 1000 iterations")
    return
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

function initializeCells(radius)
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
        # while (!isFeasible(newCentre, savedCentres, radius))
        #     newCentre = giveCentreNormalDistrInDomain(radius)
        # end
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

function createSimGif(gifPath::String,
    sol;
    title="",
    xlab="x",
    ylab="y",
    fps=3,
    dpi=300)

    if N == 0
        if radius == 0
            # PP cells 
            animSDE = @animate for i = 1:length(sol.t)

                u = sol.u[i]
                time = sol.t[i]

                X, Y = solutionToCells(u)

                timelabel = "t = $(@sprintf("%.5f", time))"

                scatter(X, Y,
                    markersize=2,
                    aspect_ratio=:equal,
                    opacity=0.25,
                    dpi=dpi,
                    title=title,
                    label=false,
                    xlims=domain,
                    ylims=domain,
                    xguidefontsize=13,
                    xlabel=timelabel)
            end
        else
            # radius > 0 -> HSCM cells 

            animSDE = @animate for i = 1:length(sol.t)

                u = sol.u[i]
                time = sol.t[i]

                X, Y = solutionToCells(u)

                circle = circleCell([X[1], Y[1]], radius)
                discreteCircleCell = cellToDiscreteCell(circle, 20)

                plot(discreteCircleCell.x, discreteCircleCell.y,
                    seriestype=:shape,
                    aspect_ratio=:equal,
                    opacity=0.25,
                    dpi=dpi,
                    title=title,
                    label=false,
                    xlims=domain,
                    ylims=domain,
                    xguidefontsize=13,
                    xlabel="t = $(@sprintf("%.5f", time))")

                for i = 2:NumberOfCells

                    circle = circleCell([X[i], Y[i]], radius)
                    discreteCircleCell = cellToDiscreteCell(circle, 20)

                    plot!(discreteCircleCell.x, discreteCircleCell.y,
                        seriestype=:shape,
                        aspect_ratio=:equal,
                        opacity=0.25,
                        dpi=dpi,
                        title=title,
                        label=false,
                        xlims=domain,
                        ylims=domain,
                        xguidefontsize=13,
                        xlabel="t = $(@sprintf("%.5f", time))")
                end

            end

        end
    else
        # DF cells 
        animSDE = @animate for i = 1:length(sol.t)

            u = sol.u[i]
            time = sol.t[i]

            X, Y = solutionToCells(u)               # now each cell is: [X[...], Y[...]]


            plot(X[1], Y[1],
                seriestype=:shape,
                aspect_ratio=:equal,
                opacity=0.25,
                dpi=dpi,
                title=title,
                label=false,
                xlims=domain,
                ylims=domain,
                xguidefontsize=13,
                xlabel="t = $(@sprintf("%.6f", time))")

            for i = 2:NumberOfCells

                plot!(X[i], Y[i],
                    seriestype=:shape,
                    aspect_ratio=:equal,
                    opacity=0.25,
                    dpi=dpi,
                    title=title,
                    label=false,
                    xlims=domain,
                    ylims=domain,
                    xguidefontsize=13,
                    xlabel="t = $(@sprintf("%.6f", time))")
            end

        end

    end

    gif(animSDE, gifPath, fps=fps)
end

struct smallSolution
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
end

function extractSolution(sol)

    if N == 0
        solSize = 2 * M
    else
        solSize = 2 * M * N
    end
    res = smallSolution(sampleTimes, [zeros(length(solSize)) for _ = 1:NumberOfSampleTimes])

    i = 1
    for t in sampleTimes
        idxInRealSol = giveClosestTimeIdx(sol, t)
        res.u[i] = sol.u[idxInRealSol]
        i += 1
    end

    return res
end

function giveClosestTimeIdx(sol, wantTime)

    if wantTime == 0.0
        return 1
    else
        for i = 2:length(sol.t)
            if sol.t[i] >= wantTime
                return i - 1
            end
        end
    end

    println("warning in giveClosestTime: wantTime", wantTime, " is larger than latest solution time")
    return length(sol.t)
end

function simulateExplicitEuler(u0)

    currentState = u0
    stateLength = length(u0)
    time = 0
    res = [zeros(stateLength) for _ = 1:length(sampleTimes)]

    NumberOfTimeSteps = floor(T / timeStepSize)
    NumberOfTimeStepsPerSave = floor(NumberOfTimeSteps / (length(sampleTimes) - 1))

    res[1] = u0
    resCount = 2
    NumberOfTimeStep = 0

    while resCount <= length(res)
        time += timeStepSize
        currentState += timeStepSize * energies(zeros(stateLength), currentState, p, time)
        NumberOfTimeStep += 1
        if mod(NumberOfTimeStep, NumberOfTimeStepsPerSave) == 0
            res[resCount] = currentState
            resCount += 1
        end
    end

    return res

end

function do1SimulationRun(simRun)
    println("simrun ", simRun)
    if N == 0
        u0 = InitializePointParticles(radius)
    else
        u0 = initializeCells(radius)
    end

    cellProb = SDEProblem(energies!, brownian!, u0, tspan, p)
    @time sol = solve(cellProb,
        EM(),
        callback=CallBack_reflectiveBC_cellOverlap,
        dt=timeStepSize,
    )
    extractedSol = extractSolution(sol)
    createLocationFile(extractedSol, simRun, locationsPath)
end

function doSimulationRuns_countLocations(currentProcss, NuSims)

    println("Starting simulations $(Int64((currentProcss-1)*NuSims+1))-$(Int64(currentProcss*NuSims)) on this core.")
    res = [zeros(Int64, NumberOfHeatGridPoints, NumberOfHeatGridPoints) for _ in 1:NumberOfSampleTimes]

    for counter = 1:nuSims

        if N == 0
            u0 = InitializePointParticles(radius)
        else
            u0 = initializeCells(radius)
        end

        cellProb = SDEProblem(energies!, brownian!, u0, tspan, p)
        sol = solve(cellProb,
            EM(),
            callback=CallBack_reflectiveBC_cellOverlap,
            dt=timeStepSize,
        )
        extractedSol = extractSolution(sol)

        for sampleTime = 1:NumberOfSampleTimes

            # extract centres from solution 
            X, Y = solutionToCells(extractedSol.u[sampleTime])
            centresX = zeros(NumberOfCells)
            centresY = zeros(NumberOfCells)
            if NumberOfCellWallPoints == 0
                centresX = X
                centresY = Y
            else
                for i = 1:NumberOfCells
                    centresX[i] = sum(X[i]) / Float64(NumberOfCellWallPoints)
                    centresY[i] = sum(Y[i]) / Float64(NumberOfCellWallPoints)
                end
            end

            # add centres to matrices
            for i = 1:NumberOfCells
                coords = [centresX[i], centresY[i]]
                row, column = getMatrixIndex(coords)
                res[sampleTime][row, column] += 1
            end

        end

    end

    return res # each matrix holds the heatmap counts for one sampleTime of all simulations 
end

function computeDesiredStates_circleCells()
    C = circleCell([0.0, 0.0], radius)
    cDF = cellToDiscreteCell(C, N)
    A_d = ones(M) * areaPolygon(cDF.x, cDF.y) # ∈ R^N
    E_d = ones(N * M)              # ∈ (R^N)^M
    I_d = ones(N * M)              # ∈ (R^N)^M
    e = computeEdgeLengths(cDF)
    ia = computeInteriorAngles(cDF)

    for i = 1:M
        E_d[(i-1)*N+1:i*N] = e
        I_d[(i-1)*N+1:i*N] = ia
    end

    return A_d, E_d, I_d
end

function runSimulation_locations()
    """
    Runs a full simulation that results in heatmaps over NumberOfSimulations simulation runs.
    In this function, the locations of the cell centres from all simulations at all sample times are saved in a .txt file. 
    """
    ### create paths 
    println("creating paths")
    mkpath(simPath)
    cp(joinpath(homedir(), "OneDrive", "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    mkpath(locationsPath)
    mkpath(heatMapsPath)

    if N != 0
        A_d, E_d, I_d = computeDesiredStates_circleCells()
    end

    ### 1st save one simulation as gif 
    println("save one sim as gif")

    if N == 0
        u0 = InitializePointParticles(radius)
    else
        u0 = initializeCells(radius)
    end

    cellProblem = SDEProblem(energies!, brownian!, u0, timeInterval, p)
    @time sol = solve(cellProblem,
        EM(),
        callback=CallBack_reflectiveBC_cellOverlap,
        dt=timeStepSize,
    )
    extractedSol = extractSolution(sol)
    createSimGif(gifPath, extractedSol)

    ### 2nd: CREATE ALL POINT LOCATIONS FOR ALL SIMULATIONS 
    results = pmap(do1SimulationRun, 1:NumberOfSimulations)

    ### 3rd: CREATE THE HEATMAP FROM ALL SIMULATION DATA 
    heatmatrices = makeMatrices()
    createHeatmaps(heatmatrices)
end

function runSimulation(NuProcs)
    """
    Runs a full simulation that results in heatmaps over NumberOfSimulations simulation runs.
    In this function, no locations are saved, but directly counted in the heat matrices and then thrown away.  
    """

    ## create paths 
    println("creating paths")
    mkpath(simPath)
    cp(joinpath(homedir(), "OneDrive", "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    mkpath(heatMapsPath)

    if N != 0
        A_d, E_d, I_d = computeDesiredStates_circleCells()
    end

    # 1st save one simulation as gif 
    println("save one sim as gif")

    if N == 0
        u0 = InitializePointParticles(radius)
    else
        u0 = initializeCells(radius)
    end

    cellProblem = SDEProblem(energies!, brownian!, u0, timeInterval, p)
    @time sol = solve(cellProblem,
        EM(),
        callback=CallBack_reflectiveBC_cellOverlap,
        dt=timeStepSize,
    )
    extractedSol = extractSolution(sol)
    createSimGif(gifPath, extractedSol)

    ## 2nd: CREATE ALL POINT LOCATIONS FOR ALL SIMULATIONS 
    @everywhere matrices = [zeros(Int64, NumberOfHeatGridPoints, NumberOfHeatGridPoints) for _ in 1:NumberOfSampleTimes]

    # distribute NumberOfSimulations 
    missingSims = mod(NumberOfSimulations, NuProcs)
    NuSims = floor(NumberOfSimulations / NuProcs)
    # then: NumberOfSimulations = NuSims*NuProcs + missingSims

    resultingMatrices = pmap(currentProcss -> doSimulationRuns_countLocations(currentProcss, NuSims), 1:NuProcs)
    matrices = doSimulationRuns_countLocations(missingSims)
    for res in resultingMatrices
        for t = 1:NumberOfSampleTimes
            matrices[t] += res[t]
        end
    end

    ### 3rd: CREATE THE HEATMAP FROM ALL SIMULATION DATA 
    createHeatmaps(matrices)

end

function runShow_overlap()
    """
    Runs a light simulation that produces a gif showing the process of resolving an overlap between two cells. 
    """

    ## create paths 
    println("creating paths")
    mkpath(simPath)
    if gethostname() == "treuesStueck"      # home pc xd 
        cp(joinpath(homedir(), "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    else # laptop 
        cp(joinpath(homedir(), "OneDrive", "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    end

    ## 1st save one simulation as gif 
    println("save one sim as gif")

    # u0 = [-0.00375, 0.15, 0.0, 0.0]  # point particle initialisation 

    # c1 = rectangleCell(Rectangle(-0.003, -0.001, -0.005, 0.005), NumberOfCellWallPoints)
    # c2 = rectangleCell(Rectangle(0.001, 0.003, -0.005, 0.005), NumberOfCellWallPoints)
    c1 = cellToDiscreteCell(circleCell([-0.004, 0.0], radius), NumberOfCellWallPoints)
    c2 = cellToDiscreteCell(circleCell([0.004, 0.0], radius), NumberOfCellWallPoints)
    u0 = [c1.x; c2.x; c1.y; c2.y]
    # A_d = 0.4 * sin(0.1 * pi) # TODO: change back
    A_d = 7.853981633974483e-5 * 2 
    E_d = 0.0015643446504023087 * ones(N)
    I_d = 0.9 * pi * ones(N)
    p = timeStepSize, D, A_d, E_d, I_d
    cellProblem = SDEProblem(energies!, nomotion!, u0, timeInterval, p)
    @time sol = solve(cellProblem,
        EM(),
        dt=timeStepSize,
    )
    extractedSol = extractSolution(sol)
   
    println("len(extractedSol) = $(length(extractedSol.t))")
    createSimGif(gifPath, extractedSol; title=simulationName)

    # for overlapScaling in [1, 10, 10^2, 10^3, round(timeStepSize^(-1))]
    # forceScalings[4] = overlapScaling 
    # cellProblem = SDEProblem(energies!, nomotion!, u0, timeInterval, p) 
    # @time sol = solve(cellProblem, 
    #                 EM(), 
    #                 dt=timeStepSize, 
    #                 )
    # extractedSol = extractSolution(sol)

    # gifPath = joinpath(simPath, string("$simulationName-scaling_$overlapScaling.gif"))
    # createSimGif(gifPath, extractedSol; title="$simulationName-scaling_$overlapScaling") 
    # end 

end