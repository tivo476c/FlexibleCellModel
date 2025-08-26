include("energies.jl")
include("cell_functionalities.jl")
# include("parameters.jl")
include("figure-creating-code/heatmap.jl")

using Distributions, Distributed, LaTeXStrings, Printf, Sockets

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

mutable struct forceArrows
    label::String
    vectors::Vector{Float64}
end

function getForces(u; A_d=[], E_d=[], I_d=[], overlap=false)
    res = forceArrows[]
    if A_d != []
        areaForceArr = forceArrows("area force", timeStepSize * forceScalings[1] * areaForce(u, A_d))
        # areaForceArr = forceArrows("area force", timeStepSize * areaForce(u, A_d))
        push!(res, areaForceArr)
    end
    if E_d != []
        edgeForceArr = forceArrows("edge force", timeStepSize * forceScalings[2] * edgeForce(u, E_d))
        # edgeForceArr = forceArrows("edge force", timeStepSize * edgeForce(u, E_d))
        push!(res, edgeForceArr)
    end
    if I_d != []
        angleForceArr = forceArrows("angle force", timeStepSize * forceScalings[3] * interiorAngleForce(u, I_d))
        # angleForceArr = forceArrows("angle force", timeStepSize * interiorAngleForce(u, I_d))
        push!(res, angleForceArr)
    end
    if overlap
        overlapForceArr = forceArrows("overlap force", timeStepSize * forceScalings[4] * overlapForce(u))
        # overlapForceArr = forceArrows("overlap force", timeStepSize * overlapForce(u))
        push!(res, overlapForceArr)
    end

    overallVecs = zeros(2 * M * N)
    for u in res
        overallVecs .+= u.vectors
    end

    overallForceArr = forceArrows("overall force", 1e-4 * overallVecs)
    push!(res, overallForceArr)

    return res
end

function addEnergyForces(plt, u;
    A_d=[],
    E_d=[],
    I_d=[],
    overlap=false)

    forceArrows = getForces(u; A_d=A_d, E_d=E_d, I_d=I_d, overlap=true)

    x_starts = u[1:N*M]
    y_starts = u[N*M+1:2*M*N]

    colors = Dict("area force" => 1, "edge force" => 2, "angle force" => 3, "overlap force" => 4, "overall force" => 5)

    for force in [forceArrows[5]]

        # compute ends of area force vectors 
        force_x_ends = x_starts + force.vectors[1:N*M]
        force_y_ends = y_starts + force.vectors[N*M+1:2*N*M]

        # x_vec = force_x_ends .- x_starts
        # y_vec = force_y_ends .- y_starts

        # quiver!(plt, x_starts, y_starts, quiver=(x_vec, y_vec), arrow=false, color=colors[force.label])
        scatter!(plt, [NaN], [NaN], label=force.label, color=colors[force.label])    # needed for label 
        for i in eachindex(x_starts)
            plot!(plt, [x_starts[i], force_x_ends[i]], [y_starts[i], force_y_ends[i]], lw=1, label=false)
        end
    end
end



function createSimGif(gifPath::String,
    sol;
    A_d=[],
    E_d=[],
    I_d=[],
    title="",
    xlab="x",
    ylab="y",
    fps=1,
    dpi=300)

    println("entered createSimGif")

    if N == 0
        if radius == 0
            # PP cells 
            animSDE = @animate for i = 1:length(sol.t)

                u = sol.u[i]
                time = sol.t[i] - timeStepSize

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
                time = sol.t[i] - timeStepSize

                X, Y = solutionToXY(u)

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
        animSDE = @animate for i = 1:NumberOfSampleTimes

            u = sol.u[i]
            time = sampleTimes[i] - timeStepSize
            X, Y = solutionToXY(u)               # now each cell is: [X[...], Y[...]]
            title = simulationName
            plt = plot(X[1], Y[1],
                seriestype=:shape,
                aspect_ratio=:equal,
                opacity=0.25,
                dpi=dpi,
                title=title,
                label=false,
                xlims=domain,
                ylims=domain,
                xguidefontsize=13,
                # xlabel="t = $(@sprintf("%.1e", time))")
                xlabel="t = $(round(Int, time/timeStepSize))e-05")

            # scatter!(plt, X[1], Y[1])

            for i = 2:NumberOfCells

                plot!(plt, X[i], Y[i],
                    seriestype=:shape,
                    aspect_ratio=:equal,
                    opacity=0.25,
                    label=false,
                )
            end

            savefig(plt, joinpath(simPath, "t$(round(Int, time/timeStepSize)).png"))

            # addEnergyForces(plt, u; A_d=A_d, E_d=E_d, I_d=I_d, overlap=true)

        end

    end

    gif(animSDE, gifPath, fps=fps)
end

function createEnergyDiagram(diaPath::String,
    sol;
    A_d=0,
    E_d=0,
    I_d=0,
    title="Energy diagram",
    xlab="time",
    ylab="energy",
    dpi=500)

    areaVector = zeros(NumberOfSampleTimes)
    edgeVector = zeros(NumberOfSampleTimes)
    angleVector = zeros(NumberOfSampleTimes)
    overlapVector = zeros(NumberOfSampleTimes)

    for i = 1:NumberOfSampleTimes
        C = solutionToCells(sol.u[i])
        for j = 1:length(C)
            areaVector[i] += areaEnergyCell(C[j], A_d)
            edgeVector[i] += edgeEnergyCell(C[j], E_d)
            angleVector[i] += angleEnergyCell(C[j], I_d)
        end
        overlapVector[i] = overlapEnergy(sol.u[i])
    end

    sampleTimes2 = sampleTimes[2:end]
    areaVector2 = areaVector[2:end]
    edgeVector2 = edgeVector[2:end]
    angleVector2 = angleVector[2:end]
    overlapVector2 = overlapVector[2:end]
    plot(sampleTimes2, areaVector2, label="Area energy", title=title, xlab=xlab, ylab=ylab, dpi=dpi)
    plot!(sampleTimes2, edgeVector2, label="Edge energy")
    plot!(sampleTimes2, angleVector2, label="Interior angle energy")
    plot!(sampleTimes2, overlapVector2, label="Overlap energy")

    savefig(diaPath)
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

    # println("warning in giveClosestTime: wantTime $wantTime is larger than latest solution time $(sol.t[end])")
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
        A_d, E_d, I_d = computeDesiredStates_circleCells()
        p = timeStepSize, D, A_d, E_d, I_d
    end

    cellProb = SDEProblem(energies!, brownian_DF!, u0, timeInterval, p, noise_rate_prototype=zeros(2 * M * N, 2 * M))
    @time sol = solve(cellProb,
        EM(),
        # callback=CallBack_reflectiveBC_cellOverlap,
        dt=timeStepSize,
    )
    # fakeSol = smallSolution([0.0], [u0])
    extractedSol = extractSolution(sol)
    createLocationFile(extractedSol, simRun, locationsPath)
    # createLocationFile(fakeSol, 1, locationsPath)
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

        cellProb = SDEProblem(energies!, brownian_DF, u0, tspan, p)
        sol = solve(cellProb,
            EM(),
            callback=CallBack_reflectiveBC_cellOverlap,
            dt=timeStepSize,
        )
        extractedSol = extractSolution(sol)

        for sampleTime = 1:NumberOfSampleTimes

            # extract centres from solution 
            X, Y = solutionToXY(extractedSol.u[sampleTime])
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
    A_d = areaPolygon(cDF.x, cDF.y)         # ∈ R
    E_d = computeEdgeLengths(cDF)           # ∈ R^N
    I_d = computeInteriorAngles(cDF)        # ∈ R^N    

    return A_d, E_d, I_d
end

function runSimulation_locations()
    """
    Runs a full simulation that results in heatmaps over NumberOfSimulations simulation runs.
    In this function, the locations of the cell centres from all simulations at all sample times are saved in a .txt file. 
    """
    ## create paths 
    # println("creating paths")
    # mkpath(simPath)
    # if gethostname() == "treuesStueck"      # home pc xd 
    #     cp(joinpath(homedir(), "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    # elseif gethostname() == "iwr-25177394"
    #     cp(joinpath(homedir(), "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)    # use correct uniPC path
    # else # laptop 
    #     cp(joinpath(homedir(), "OneDrive", "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    # end
    # mkpath(heatMapsPath)
    # mkpath(locationsPath)

    # if N != 0
    #     A_d, E_d, I_d = computeDesiredStates_circleCells()
    #     p = timeStepSize, D, A_d, E_d, I_d 
    # end

    # # 1st save one simulation as gif 
    # println("save one sim as gif")

    # if N == 0
    #     u0 = InitializePointParticles(radius)
    # else
    #     u0 = initializeCells(radius)
    # end

    # cellProblem = SDEProblem(energies!, brownian_DF!, u0, timeInterval, p, noise_rate_prototype=zeros(2 * M * N, 2 * M))
    # @time sol = solve(cellProblem,
    #     EM(),
    #     # callback=CallBack_reflectiveBC_cellOverlap,
    #     dt=timeStepSize,
    # )
    # extractedSol = extractSolution(sol)
    # createSimGif(gifPath, extractedSol)

    ### 2nd: CREATE ALL POINT LOCATIONS FOR ALL SIMULATIONS 
    # results = pmap(do1SimulationRun, 1:NumberOfSimulations)

    ### 3rd: CREATE THE HEATMAP FROM ALL SIMULATION DATA 
    heatmatrices = makeMatrices()
    # smoothenMatrix!(heatmatrices, 30)
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

    # if N != 0
    #     A_d, E_d, I_d = computeDesiredStates_circleCells()
    #     p = timeStepSize, D, A_d, E_d, I_d 
    # end

    # 1st save one simulation as gif 
    println("save one sim as gif")

    if N == 0
        u0 = InitializePointParticles(radius)
    else
        u0 = initializeCells(radius)
    end

    cellProblem = SDEProblem(energies!, brownian_DF, u0, timeInterval, p)
    @time sol = solve(cellProblem,
        EM(),
        # callback=CallBack_reflectiveBC_cellOverlap,
        dt=timeStepSize,
    )
    createSimGif(gifPath, sol)

    # extractedSol = extractSolution(sol)
    # createSimGif(gifPath, extractedSol)

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
    # println("save one sim as gif")

    ###########################################################################################
    #### deforming overlap force config 
    # c1 = cellToDiscreteCell(circleCell([-0.0025, 0.0], radius), 6) 
    # c2 = cellToDiscreteCell(circleCell([0.0025, 0.0], radius), 6; rotation=pi/6.0) 
    # u0 = [c1.x; c2.x; c1.y; c2.y]

    ### area force config 
    # c1 = DiscreteCell(0.7.*[0.01, 0.0, -0.01, -0.01, 0.0, 0.01], 0.1.*[0.0025, 0.005, 0.0025, -0.0025, -0.005, -0.0025])
    # u0 = [c1.x; c1.y]

    ### edge force config 
    # c1 = DiscreteCell([0.01, 0.0005, -0.0005, -0.01, -0.0005, 0.0005], [0.0, 0.004, 0.004, 0.0, -0.004, -0.004])
    # u0 = [c1.x; c1.y]

    ### interior angle force config 
    # c1 = DiscreteCell([0.003, 0.009, -0.003, -0.009, -0.003, 0.009], [0.0, 0.003, 0.003, 0.0, -0.003, -0.003])
    # u0 = [c1.x; c1.y]

    #### deforming overlap force config 
    c1 = cellToDiscreteCell(circleCell([-0.006, 0.0], radius), 6)
    c2 = cellToDiscreteCell(circleCell([0.006, 0.0], radius), 6; rotation=pi / 6.0)
    u0 = [c1.x; c2.x; c1.y; c2.y]


    ###########################################################################################
    A_d = circleArea(radius, N)
    E_d = circleEdgeLengths(radius, N)
    I_d = circleInteriorAngles(N)
    p = timeStepSize, D, A_d, E_d, I_d
    cellProblem = SDEProblem(energies!, nomotion!, u0, timeInterval, p)
    @time sol = solve(cellProblem,
        EM(),
        dt=timeStepSize,
    )
    extractedSol = extractSolution(sol)
    createSimGif(gifPath, extractedSol; A_d=A_d, E_d=E_d, I_d=I_d, title=simulationName)
    # createSimGif(gifPath, extractedSol; A_d=A_d, E_d=E_d, I_d=I_d, title="", fps=1)
    createEnergyDiagram(energyDiaPath, extractedSol; A_d=A_d, E_d=E_d, I_d=I_d)

    # for intAngleScale in [5e0 ,1e1, 5e1,1e2]

    # for overlapScaling in [1e3, 1e4, 1e5]
    #     for overlapScaling in [2e4]
    #         forceScalings[4] = overlapScaling

    #         intanglstring = @sprintf("%.1e", forceScalings[3])
    #         overlapstring = @sprintf("%.1e", forceScalings[4])
    #         nameSim = string("drift_9_4_$(intanglstring)_$(overlapstring)")
    #         simPath = joinpath(homedir(), "simulations", nameSim)
    #         gifPath = joinpath(simPath, string(nameSim, ".gif"))
    #         energyDiaPath = joinpath(simPath, "energies-$nameSim.png")
    #         mkpath(simPath)
    #         if gethostname() == "treuesStueck"      # home pc xd 
    #             cp(joinpath(homedir(), "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    #         else # laptop 
    #             cp(joinpath(homedir(), "OneDrive", "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    #         end

    #         cellProblem = SDEProblem(energies!, nomotion!, u0, timeInterval, p)
    #         @time sol = solve(cellProblem,
    #             EM(),
    #             dt=timeStepSize,
    #         )
    #         extractedSol = extractSolution(sol)
    #         createSimGif(gifPath, extractedSol; title=nameSim)
    #         createEnergyDiagram(energyDiaPath, extractedSol; A_d=A_d, E_d=E_d, I_d=I_d)
    #     end
    # end

end
