include("energies.jl")
include("heatmap.jl")
include("parameters.jl")



function createSimGif(  gifPath::String, 
    sol, 
    title="", 
    xlab="x",
    ylab="y",
    fps=5,
    dpi=300)

    # is length(sol) == T / dt or length(sol) == T / saveAtTimes ???

    animSDE = @animate for t = 1:length(sol)
        time = t[1]
        # time = t[1]*timeStepSize+1
        # time = t[1]*saveAtTimes+1
        X, Y = solutionToCells(sol[time])
        x = X[1]
        y = Y[1]

        timelabel = "t = $(@sprintf("%.2f", time*timeStepSize))"  

        plot(x, y, 
            seriestype = :shape, 
            aspect_ratio = :equal, 
            opacity = 0.25, 
            dpi = dpi, 
            title = title,
            label = false,
            xlims = domain,
            ylims = domain,
            xguidefontsize=13,
            xlabel = timelabel)

        for i = 2:M 
            plot!(X[i],Y[i], seriestype=:shape, opacity=.25, label = false)
        end  

    end

    gif(animSDE, gifPath, fps = fps)

end 

function nullnummer(du, u, p, t) 
    return zeros(length(du))
end 

begin
    include("energies.jl")
    include("parameters.jl")
    include("heatmap.jl")

    simPath = joinpath(homedir(), "simulations", simulationName)
    locationsPath = joinpath(simPath, "locations")
    heatMapsPath = joinpath(simPath, "heatmaps")

    C = circleCell([-1.875,1.875], 0.3)
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

    u0 = zeros(2*M*N) 

    for i = 0:3
        for j = 0:3
            c = moveC(cDF, j*1.25, -i*1.25)
            u0[ 1+N*(j+4i): N*(1+j+4*i)] = c.x
            u0[ 1+N*(j+4i+M): N*(1+j+4*i+M) ] = c.y
        end     
    end 

    tspan = timeInterval
    Δt = timeStepSize
    p = [Δt, D]
    prob_cell1 = SDEProblem(energies, brownian, u0, tspan, p, noise=WienerProcess(0., 0.))        
    @time sol = solve(prob_cell1, EM(), dt=timeStepSize)
    
    #sol = solve(prob_cell1, EM(), dt=Δt, saveat=Δt)
    #println("saveat = deltaT: ", length(sol))
    #sol = solve(prob_cell1, EM(), dt=Δt, saveat=saveAtTimes)
    #println("saveat = saveAtTimes: ", length(sol))
    
    
    # HEATMAP stuff begins 
    
    mkpath(simPath) 
    #if !isfile(joinpath(simPath, "parameters.jl"))
    cp(joinpath(homedir(), "OneDrive", "Desktop", "Master-Thesis", "code", "parameters.jl"), joinpath(simPath, "parameters.jl"), force=true)
    #end 
    mkpath(locationsPath)
    mkpath(heatMapsPath)
    
    # save one simulation as give 
    gifPath = joinpath(simPath, string(simulationName, ".gif"))
    createSimGif(gifPath, sol)
    steps = timeStepSize
    for i = 1:NumberOfSimulations
        sol = solve(prob_cell1, RKMil(), dt=Δt)
        createLocationFile(sol, i, locationsPath)
    end 

    matrices = makeMatrices()

    createHeatmaps(matrices)


end


