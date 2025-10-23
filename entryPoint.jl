using Pkg
"""
This file is the entry point of my FlexibleCellModel code. 
It gets executed, whenever parameters.jl gets executed.
It uses its parameter configuration to start the according simulation. 
"""

println("Started entryPoint.jl")
println("Loading all packages")


#######  OPTION 1: JUST RUN A METHOD WITHOUT PARALLELIZING
include("code/parameters.jl")
include("code/cell_functionalities.jl")
include("code/computeOverlap.jl")
include("code/energies.jl")
include("code/figure-creating-code/heatmap.jl")
include("code/simulationFunctionalities.jl")

tspan = timeInterval
simPath = joinpath(homedir(), "simulations", simulationName)
locationsPath = joinpath(simPath, "locations")
heatMapsPath = joinpath(simPath, "heatmaps")
gifPath = joinpath(simPath, string(simulationName, ".gif"))
energyDiaPath = joinpath(simPath, "energies-$simulationName.png")
p = [timeStepSize, D]


# runSimulation_locations()                                   # for producing location data and heatmaps

#######  OPTION 2: PARALLELIZED 
# if wanna use parallelized run, inclusions happen in startParallelizedRun.jl: 

# using Pkg 
# Pkg.activate(".")
# Pkg.instantiate()
# include("code/startParallelizedRun.jl")


#######  OPTION 3: PLAYGROUND 

# runShow_overlap()                                         # for producing overlap figures and gifs 
sims = ["AAAsoftSim-bachelorOverlap", "AAAmidSim3-bachelorOverlap", "AAAtest3-hard-DF-FIRSTWORKING"]

allPlots = []
for simname in sims
    simulationName = simname
    heatmatrices = makeMatrices()
    plots = createCrossSection(heatmatrices; simname=simname)
    push!(allPlots,plots)
end 

combinedPlots = []
YdataSoft = []
YdataMid  = []
YdataHard = []
for t=1:6
    time = 0.01*t 
    push!(YdataSoft, allPlots[1][t].series_list[1][:y])
    push!(YdataMid , allPlots[2][t].series_list[1][:y])
    # push!(YdataHard, Plots.series_list(allPlots[3][t][1])[:y])
end 

println(" allPlots[1][t].series_list[1][:y] = $(allPlots[1][5].series_list[1][:y])")
println(" allPlots[2][t].series_list[1][:y] = $(allPlots[2][5].series_list[1][:y])")
YdataSoft == YdataMid 

begin
    t = 0.01
    caption = string("x\n t = ", @sprintf("%.2f", t))
    plot(heatcells, 
        YdataSoft[5],
        title="Cross section density",
        label="h = 0.0",
        color=1,
        xlimits=domain,
        ylimits=(0.0,4.0),
        xlab=caption,
        ylab="Density",
        dpi=500)

    plot!(heatcells, 
          YdataMid[5],
          label="h = 0.5"
        )

    plot!(heatcells, 
          YdataHard[5],
          label="h = 1.0"
        )


end

begin
    heatcells = HeatGrid[1:end-1] .+ 0.5 * HeatStepSize

    for t=1:6 
        time = 0.01*(t-1)
        caption = string("x\n t = ", @sprintf("%.2f", time))
        onetimeplot = plot(heatcells, 
                    YdataMid[t],
                    title="Cross section density",
                    label="h = 0.5",
                    color=2,
                    xlimits=domain,
                    ylimits=(0.0,4.0),
                    xlab=caption,
                    ylab="Density",
                    dpi=500)

        # plot!(onetimeplot, heatcells, YdataMid[t],
        #     label="h = 0.5",
        #     color=2,
        # )

        # plot!(onetimeplot, heatcells, YdataHard[t],
        #     label="h = 1.0",
        #     color=3,
        # )
        savefig(onetimeplot, joinpath(homedir(), "simulations", "crosssections-collection", "cross_time$(time).png"))
        push!(combinedPlots, onetimeplot)

    end 

    anim = @animate for p in combinedPlots
        plot(p)
    end
    # save as GIF
    gifname="combined-cross-section-evolution.gif"
    gif(anim, joinpath(homedir(), "simulations", "crosssections-collection", gifname), fps=1)
end


