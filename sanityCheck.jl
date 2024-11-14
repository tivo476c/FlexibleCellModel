begin
    include("energies.jl")
    include("SaveCoordsModule.jl")

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
    prob_cell1 = SDEProblem( energies, brownian, u0, tspan, p, noise=WienerProcess(0., 0.))        
    steps = Δt



    
    sol = solve(prob_cell1, EM(), dt=Δt, saveat=steps)

    simPath = joinpath(homedir(), "simulations", simulationName)
    locationsPath = joinpath(simPath, "locations")
    heatMapsPath = joinpath(simPath, "heatmaps")

    mkpath(simPath) 
    if !isfile(joinpath(simPath, "Parameters.jl"))
        cp(joinpath(pwd(), "OneDrive", "Desktop", "BA-Code", "Parameters.jl"), joinpath(simPath, "Parameters.jl"))
    end 
    mkpath(locationsPath)
    mkpath(heatMapsPath)

    for i = 1:NumberOfSimulations
        sol = solve(prob_cell1, EM(), dt=Δt, saveat=steps)
        createLocationFile(sol, i)
    end 

    matrices = makeMatrices()

    createHeatmaps(matrices)
end 




#=
#------------------------------------- ANIMATION

xdomain =(-6,6)
ydomain = (-6,6)
titte = "Sanity Check"

animSDE = @animate for t = 1:floor(length(sol)/4.0)

    time = Int(4*t)
    x = sol[2][1:N]  
    y = sol[time][ M*N+1 : M*N + N]
    X,Y = solutionToCells(sol[time]) 

    xlab = string("t = ", time)
    plot(X[1], Y[1], seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=xdomain, ylims = ydomain, xguidefontsize=15, xlabel = xlab, title = titte)

    for i = 2:M

        plot!(X[i],Y[i], seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false )

    end 

end

#=
animSDE = @animate for t = 1:floor(length(sol)/4.0)

    time = Int(4*t)
    x = sol[2][1:N]  
    y = sol[time][ M*N+1 : M*N + N]

    xlab = string("t = ", time)
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=xdomain, ylims = ydomain, xguidefontsize=15, xlabel = xlab, title = titte)

    for i = 2:M

        x = sol[time][ (i-1)*N+1  : i*N]
        y = sol[time][ (i-1)*N+1+N*M  : i*N+N*M]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false )

    end 

end
=#

gif(animSDE, "findParameters.mp4", fps = 3)
=#
