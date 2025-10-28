using Plots
using Printf

include("../cell_functionalities.jl")
include("../energies.jl")

Nvertices = 1000
radius = 0.005    
cell = cellToDiscreteCell(circleCell([0.0,0.0], 0.005), Nvertices)
area = pi * radius^2
perimeter = 2 * pi * radius
asphericity = 4 * pi * area / (perimeter^2)

lims = (-0.007, 0.007)
xlab = "\nx\n\narea = $(@sprintf("%.2e", area)),\nperimeter = $(@sprintf("%.2e", perimeter)),\nasphericity = $(@sprintf("%.2e", asphericity))"
plot(
    cell.x,
    cell.y,
    label=false,
    xlab=xlab,
    ylab="y", 
    seriestype=:shape,
    aspectratio=:equal,
    xlims=lims,
    ylims=lims, 
    xticks=[-0.005, 0.0, 0.005],
    yticks=[-0.005, 0.0, 0.005],
    dpi=500,
    framestyle =:box,
    legendfontsize = 10,
    # xguidefontsize = 16,     # X label font size
    # yguidefontsize = 16,     # Y label font size
    xtickfontsize = 10,      # X tick label font size
    ytickfontsize = 10,      # Y tick label font size
    size=(400,400),
)

aspPath = joinpath(homedir(), "simulations", "asphericities")
plotName = "asp_circle.png"
savefig(joinpath(aspPath, plotName))


