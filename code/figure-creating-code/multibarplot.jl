using StatsPlots


xLowerBound = 0.50
xUpperBound = 0.95
Nintervals = 10
xStepSize = (xUpperBound - xLowerBound) / Nintervals
factor = xStepSize^(-1)
aspValues = xLowerBound:xStepSize:xUpperBound
xScale  = repeat(["$(aspValues[i]) - $(aspValues[i+1])" for i=1:Nintervals], outer=3)  

asp1 = 50*ones(10)
asp2 = 30*ones(10)
asp3 = 70*ones(10)
yData = hcat(asp1, asp2, asp3) 
allLabels = repeat(["h = 0", "h = 0.5", "h = 1"], inner=10)


groupedbar(xScale, yData, 
            group = allLabels,
            xlabel = "Asphericity",
            ylabel = "Number of cells",
            ylims = (0, maximum(yData) * 1.2),
            xrotation = 45, 
            bar_position =:dodge, 
            bar_width = 0.85, 
            legend =:topright,
            legendfontsize = 12,
            xguidefontsize = 12,     # X label font size
            yguidefontsize = 12,     # Y label font size
            xtickfontsize = 10,      # X tick label font size
            ytickfontsize = 10,      # Y tick label font size
            framestyle = :box,
            dpi=500,
            size=(900,800),
            )
aspPath = joinpath(homedir(), "simulations", "asphericities")
barname = "test.png"
savefig(joinpath(aspPath, barname))