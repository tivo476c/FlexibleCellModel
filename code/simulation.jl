using DifferentialEquations, StochasticDiffEq, Plots

# include("energies.jl")
include("parameters.jl")


function createSimGif(  gifPath::String, 
                        sol, 
                        title="", 
                        label="", 
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
    
    gif(animSDE, gifPath, fps = fps)

end 