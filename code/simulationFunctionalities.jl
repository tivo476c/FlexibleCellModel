include("energies.jl")
include("parameters.jl")

using Printf

function createSimGif(  gifPath::String, 
    sol, 
    title="", 
    xlab="x",
    ylab="y",
    fps=5,
    dpi=300)

    # is length(sol) == T / dt or length(sol) == T / saveAtTimes ???
    if N != 0
        # HSCM cells 
        animSDE = @animate for i = 1:length(sol.t) 
            # time = t[1]*timeStepSize+1
            # time = t[1]*saveAtTimes+1

            u = sol.u[i]
            time = sol.t[i]

            X, Y = solutionToCells(u)
            x = X[1]
            y = Y[1]

            timelabel = "t = $(@sprintf("%.2f", time))"  

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
    else 
        # Point particles 
        animSDE = @animate for i = 1:length(sol.t) 
            # time = t[1]*timeStepSize+1
            # time = t[1]*saveAtTimes+1
           
            u = sol.u[i]
            time = sol.t[i]

            X, Y = solutionToCells(u)
           
            timelabel = "t = $(@sprintf("%.4f", time))"  

            scatter(X, Y,
                aspect_ratio = :equal, 
                dpi = dpi, 
                title = title,
                label = false,
                xlims = domain,
                ylims = domain,
                xguidefontsize=13,
                xlabel = timelabel)

        end

        gif(animSDE, gifPath, fps = fps)
    end 

end 

struct smallSolution
    t::Vector{Float64}  
    u::Vector{Vector{Float64}}
end 

function extractSolution(sol)

    if N==0
        solSize = 2*M
    else 
        solSize = 2*M*N
    end
    res = smallSolution(sampleTimes, [zeros(length(solSize)) for _=1:NumberOfSampleTimes])

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
            if sol.t[i] > wantTime
                return i-1
            end 
        end
    end

    println("warning in giveClosestTime: wantTime", wantTime," is larger than latest solution time")
    return length(sol.t)
end 

function simulateExplicitEuler(u0)

    currentState = u0
    stateLength = length(u0)
    time = 0 
    res = [zeros(stateLength) for _ = 1:length(sampleTimes)]
    
    NumberOfTimeSteps = floor(T / timeStepSize)
    NumberOfTimeStepsPerSave = floor(NumberOfTimeSteps/(length(sampleTimes)-1))

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