include("energies.jl")
include("parameters.jl")

# ------------------ OVERLAP PLOT 

C = Cell([3.75, 5.0], x -> 3.0*(1 - 0.3*cos(2*x)))
cDF = cellToDiscreteCell(C, N) 
shapePlot(cDF)
cDF2 = moveC(cDF, 2.5, 0.0)
u0 = [cDF.x; cDF2.x; cDF.y; cDF2.y]


p = [timeStepSize, D]
prob_cell1 = SDEProblem( energies, nomotion, u0, timeInterval, p, noise=WienerProcess(0., 0.))        # how to correctly pass Δt and D in p ?
sol = solve(prob_cell1, RKMil(), dt=timeStepSize)
domain = (1.0, 9.0)

animSDE = @animate for t ∈ 0: round(Int64, length(sol)/4 - 1 )  

    time = t[1]*4+1
    x = sol[time][1:N]
    y = sol[time][ M*N+1 : N*(M+1)]
    ct = DiscreteCell(x,y)
    #lab = string("time: ", time)
    area = areaPolygon(x,y)
    a1 = round(  intAngle2( vertex(ct, 2), vertex(ct, 3),vertex(ct, 4) ) / π * 180 , digits=1)
    a2 = round(  intAngle2( vertex(ct, 5), vertex(ct, 6),vertex(ct, 1) ) / π * 180 , digits=1)

    xlab = string( "ι3(C(", time ,")) = ", a1, "°; ι6(C(", time ,")) = ", a2, "°"  )
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=domain, ylims = domain, xguidefontsize=13, xlabel = xlab )
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[time][ (i-1)*N+1  : i*N]
        y = sol[time][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 
    
end

gif(animSDE, fps = 10)