using Printf

function createSimGif(  gifPath::String, 
    sol, 
    title="", 
    xlab="x",
    ylab="y",
    fps=5,
    dpi=300)

    # is length(sol) == T / dt or length(sol) == T / saveAtTimes ???

    animSDE = @animate for i = 1:length(sol) 
        # time = t[1]*timeStepSize+1
        # time = t[1]*saveAtTimes+1

        u = sol.u[i]
        time = sol.t[i]

        X, Y = solutionToCells(u)
        x = X[1]
        y = Y[1]

        timelabel = "t = $(@sprintf("%.4f", time))"  

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