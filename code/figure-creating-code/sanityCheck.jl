"""
This file does:
* simulate the point particles thingy like in chapman12 
* simulate the the model for our flexible cell model
"""

include("../parameters.jl")
include("sanityCheckFunctionalitites.jl")
include("heatmap.jl")
# addprocs(5)
addprocs(3)

println("hahahihi")
begin
    """
    This loop computes one Heatmap from NoSims amount of simulation 
    """

    ## get all needed paths 
    @everywhere begin
        include("../parameters.jl")
        include("sanityCheckFunctionalitites.jl")

        ### 1st: DO PRIOR WORK 
        tspan = timeInterval
        simPath = joinpath(homedir(), "simulations", simulationName)
        locationsPath = joinpath(simPath, "locations")
        heatMapsPath = joinpath(simPath, "heatmaps")
        gifPath = joinpath(simPath, string(simulationName, ".gif"))
        p = [timeStepSize, D]
    end

    ## create paths 
    println("Creating paths")
    mkpath(simPath)
    cp(joinpath(homedir(), "OneDrive", "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    mkpath(locationsPath)
    mkpath(heatMapsPath)


    ## INITIALIZATION 
    C = circleCell([0.0, 0.0], radius)
    cDF = cellToDiscreteCell(C, N)
    # A1 = ones(M) * areaPolygon(cDF.x, cDF.y) # ∈ R^N
    # E1 = ones(N * M)              # ∈ (R^N)^M
    # I1 = ones(N * M)              # ∈ (R^N)^M
    # e = computeEdgeLengths(cDF)
    # ia = computeInteriorAngles(cDF)

    # for i = 1:M
    #     E1[(i-1)*N+1:i*N] = e
    #     I1[(i-1)*N+1:i*N] = ia
    # end

    ## save one simulation as gif 
    println("Save a simulation as .gif")
    u0 = initializeChapman(radius)
    prob_pointParticles = SDEProblem(energies, brownian, u0, tspan, p, noise=WienerProcess(0.0, 0.0))
    sol = solve(prob_pointParticles, EM(), dt=timeStepSize, saveat=sampleTimes)
    createSimGif(gifPath, sol)

    tspan = timeInterval
    p = [timeStepSize, D]

    ## NOW START ALL SIMULATION RUNS 
    results = pmap(doAHSCMSimulationRun, 1:NumberOfSimulations)

    matrices = makeMatrices()
    createHeatmaps(matrices)

end



include("sanityCheckFunctionalitites.jl")
include("../parameters.jl")
include("../simulationFunctionalities.jl")
# addprocs(6)
addprocs(3)
println("hallo")
begin
    ##### PARALLELIZED CREATION OF POINT PARTICLE HEAT MAP 
    @everywhere begin
        include("heatmap.jl")
        include("sanityCheckFunctionalitites.jl")
        include("../parameters.jl")
        include("../energies.jl")
        include("../simulationFunctionalities.jl")


        ### 1st: DO PRIOR WORK 
        tspan = timeInterval
        simPath = joinpath(homedir(), "simulations", simulationName)
        locationsPath = joinpath(simPath, "locations")
        heatMapsPath = joinpath(simPath, "heatmaps")
        gifPath = joinpath(simPath, string(simulationName, ".gif"))
        p = [timeStepSize, D]
    end

    ## create paths 
    println("creating paths")
    mkpath(simPath)
    cp(joinpath(homedir(), "OneDrive", "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    mkpath(locationsPath)
    mkpath(heatMapsPath)

    ## save one simulation as gif 
    println("save one sim as gif")

    u0 = InitializePointParticles(radius)

    prob_pointParticles = SDEProblem(energies!, brownian!, u0, tspan, p) 
    @time sol = solve(prob_pointParticles, 
                      EM(), 
                      callback=CallBack_reflectiveBC, 
                      dt=timeStepSize, 
                    #   saveat=sampleTimes, 
                    #   tstops=sampleTimes,
                    #   dense=false
                    )
    extractedSol = extractSolution(sol)
    createSimGif(gifPath, extractedSol) 

    ### 2nd: CREATE ALL POINT LOCATIONS FOR ALL SIMULATIONS 
    results = pmap(doAPointParticleSimulationRun, 1:NumberOfSimulations)

    ### 3rd: CREATE THE HEATMAP FROM ALL SIMULATION DATA 
    matrices = makeMatrices()
    createHeatmaps(matrices)
end

# do own explicit euler 
begin 
    include("heatmap.jl")
    include("sanityCheckFunctionalitites.jl")
    include("../parameters.jl")
    include("../simulationFunctionalities.jl")
    include("../energies.jl")

    tspan = timeInterval
    simPath = joinpath(homedir(), "simulations", simulationName)
    locationsPath = joinpath(simPath, "locations")
    heatMapsPath = joinpath(simPath, "heatmaps")
    gifPath = joinpath(simPath, string(simulationName, ".gif"))
    p = [timeStepSize, D]

    ## create paths 
    println("creating paths")
    mkpath(simPath)
    cp(joinpath(homedir(), "OneDrive", "Desktop", "FlexibleCellModel", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    mkpath(locationsPath)
    mkpath(heatMapsPath)

    matrices = [zeros(Int64, NumberOfHeatGridPoints, NumberOfHeatGridPoints) for _ in 1:NumberOfSampleTimes]
    for i = 1:NumberOfSimulations
        u0 = InitializePointParticles(radius)
        res = simulateExplicitEuler(u0)
        addSolToMatrices(res, matrices)
    end 
    createHeatmaps(matrices)
end






