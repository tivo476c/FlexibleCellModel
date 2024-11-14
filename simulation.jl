using DifferentialEquations, StochasticDiffEq

# include("energies.jl")
include("parameters.jl")

"""
    createSimGif

solves problem and creates and saves gif at gifPath
"""
function createSimGif(  gifPath::String, 
                        problem::SDEProblem, 
                        title="", 
                        label="", 
                        xlab="x",
                        ylab="y",
                        fps=5,
                        dpi=300)

    sol = solve(problem, EM(), dt=timeStepSize, saveat=saveAtTimes)
    # is length(sol) == T / dt or length(sol) == T / saveAtTimes ???

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
        dpi = dpi, 
        title = title,
        label = label,
        xlims = domain,
        ylims = domain,
        xguidefontsize=13,
        xlabel = xlab,
        ylabel = ylab)
    
        for i = 2:M 
            plot!(X[i],Y[i], seriestype=:shape, opacity=.25)
        end  

    end
    
    gif(animSDE, fps = fps)
    savefig(gifPath)

end 