include("energies.jl")
include("simulation.jl")
include("parameters.jl")

using Plots 
# SETUP: 2 cells overlapping 

C = circleCell([-1.5,0.0], 2.0)
cDF1 = cellToDiscreteCell(C, N) 
cDF2 = moveC(cDF1, 3.0, 0.0)

A1 = ones(M) * areaPolygon(cDF1.x, cDF1.y) # ∈ R^N
E1 = ones(N*M)              # ∈ (R^N)^M
I1 = ones(N*M)              # ∈ (R^N)^M
e = computeEdgeLengths(cDF1)
ia = computeInteriorAngles(cDF1)

# initial condition 

for i = 1:M
    E1[(i-1)*N+1 : i*N] = e
    I1[(i-1)*N+1 : i*N] = ia
end 

u0 = vcat(cDF1.x, cDF2.x, cDF1.y, cDF2.y)

p = [timeStepSize, D]
# sus:
prob_cell1 = SDEProblem( energies, brownian, u0, timeInterval, p, noise=WienerProcess(0., 0.))        

overlapProblem = ODEProblem(energies, u0, timeInterval)

gifPath = joinpath(homedir(), "showOverlap", "overlap-combination.gif")
#createSimGif(gifPath, prob_cell1)

# sol = solve(prob_cell1, EM(), dt=timeStepSize, saveat=saveAtTimes)

sol = solve(overlapProblem)


animSDE = @animate for t = 1:length(sol)
    # animSDE = @animate for t = 0:length(sol)

    time = t[1]
    # time = t[1]*timeStepSize+1
    # time = t[1]*saveAtTimes+1
    X, Y = solutionToCells(sol[time])
    x = X[1]
    y = Y[1]

    plot(x, y, 
    seriestype = :shape, 
    aspect_ratio = :equal, 
    opacity = 0.25, 
    dpi = 300,
    xlims = domain,
    ylims = domain)

    for i = 2:M 
        plot!(X[i],Y[i], seriestype=:shape, opacity=.25)
    end  
end


gif(animSDE, gifPath, fps = 5)





