using Plots
using Printf

include("../cell_functionalities.jl")
include("../energies.jl")

cell = cellToDiscreteCell(circleCell([0.0,0.0], 0.005), 3)
area = areaPolygon(cell.x, cell.y)
perimeter = sum(computeEdgeLengths(cell))
asphericity = 4 * pi * area / (perimeter^2)

lims = (-0.007, 0.007)
label = "area = $(@sprintf("%.2f", area)) \nperimeter = $(@sprintf("%.2f", perimeter)) \nasphericity = $(@sprintf("%.2f", asphericity))"
plot(
    cell.x,
    cell.y,
    label=label,
    xlab="x",
    ylab="y", 
    seriestype=:shape,
    xlims=lims,
    ylims=lims, 
    dpi=500,
)


