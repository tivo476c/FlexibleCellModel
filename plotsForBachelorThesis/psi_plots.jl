include("../cell_functionalities.jl")
#include("../BakeTheLists4.jl")
include("../BakeTheLists5.jl")

using LinearAlgebra, Plots, ColorSchemes 

psi(x, y) = 0.5 * (1 + tanh(500*(1 - norm([x, y], 2) )))

function indicator(x,y)

    if(norm([x,y], 2) <= 1)
        return 1
    else 
        return 0
    end 

end

function approximation(h:: Float64)
    N = floor(2.0 / h) 
    res = 0
    m0 = [-1.0, -1.0] - 0.5*h*[1.0, 1.0]
    for i = 1:N 
        for j = 1:N
            
            m   = m0 + h*[i, j]
            res = res + psi(m[1], m[2])

        end
    end 

    return h^2 * res 

end

function realArea(x)
    return pi 
end

r1(x) = 1 - 0.2 * cos(2*x)
r2(x) = 0.5*sin(2*x)+1

function psi1(x, y)

    ϕ = computeAngle([0.0, 0.0], [x,y]) 
    s = 500
    return 0.5*(1 + tanh(s*(r1(ϕ)- norm([x,y],2))))

end 
function psi2(x, y)

    ϕ = computeAngle([0.0, 0.0], [x,y]) 
    s = 500
    return 0.5*(1 + tanh(s*(r2(ϕ)- norm([x,y],2))))

end 

psi3(x,y) = psi1(x,y) * psi2(x,y)

X = (-1.5 : 0.001 : 1.5)  
Y = X

# 3D Plots of the figure which should illustrate how the integral of psi1 * psi2 computes the overlapping area _______________
surface(X, Y, psi3, seriescolor=cgrad([:orange, :blue], [0.1, 0.3, 0.8]), 
        xtickfontsize=10,
        ytickfontsize=10,
        ztickfontsize=11,
        dpi=500)
#contourf for both
contour(X, Y, psi1, seriescolor=cgrad([:orange, :blue], [0.1, 0.3, 0.8]), dpi=500, aspect_ratio=:equal,
xtickfontsize=11,
ytickfontsize=11,
ztickfontsize=11)
savefig("DisplayContOverlap2_contour.png")


#psi plots that demonstrate, how s changes the psifunction 
psi5(x, y) = 0.5 * (1 + tanh(5*(1 - norm([x, y], 2) )))
psi20(x, y) = 0.5 * (1 + tanh(5*(20 - norm([x, y], 2) )))
psi500(x, y) = 0.5 * (1 + tanh(5*(500 - norm([x, y], 2) )))
function indicator(x,y)

    if(norm([x,y], 2) <= 1)
        return 1
    else 
        return 0
    end 

end
X = (-2.0 : 0.005 : 2.0) 
psi.(X,X)
using Plots, ColorSchemes
Y = X
surface(X, Y, indicator, seriescolor=cgrad([:orange, :blue], [0.1, 0.3, 0.8]))
#contour( X, Y, indicator, aspect_ratio=:equal)
#contour( X, Y, psi500, aspect_ratio=:equal)



# 2D Plots ___________________________________________________________________________________________________________________
f_unitcircle(x) = 0.5 * (1 + tanh(50*(1 - abs(x) )))

# plot for psi function for different s *****************************************************
plot!(f_unitcircle, label = "s = 50",
linewidth=3,
tickfontsize=13,
        dpi=500, legendfont=font(14)
        )
savefig("preSmoothIndicator")

# plot for error estimation of summed up midpoint rule ***************************************
plot( 0.001:0.001:1, realArea, label="Unit Circle Area", xlabel = "h", legendfont=font(14), 
    linewidth=2,
    guidefontsize=16,
    tickfontsize=13
)
plot!(approximation, label="Approximated Area", 
    linewidth=2,
    tickfontsize=13
    )
savefig("quadrature_approximation.png")


# cell plots 
C = Cell([0.0, 0.0], x -> 1 - 0.2 * cos(2*x))
D = Cell([0.0, 0.0], x -> 0.5*sin(2*x)+1)
multiShapePlot(MutableLinkedList{Cell}(C,D), MutableLinkedList{String}("C", "ζ"))
savefig("DisplayContOverlap1.png")


c1 = pointsToArray( outerWall( C, 0.01 ) )  
plot(c1[1,:], c1[2,:], aspect_ratio=:equal, label = false,
    xlabel = "x",
    ylabel = "y",
    lw = 2, 
    guidefontsize=12,
    xlims = (0.0, 1.0),
    ylims = (0.0, 1.5),
    ticks = false,
    dpi=500)

savefig("DisplayContOverlap4.png") 
plot!()


# figure rOverlap 
C = cellToDiscreteCell(peanutCell([1.0, 0.0],0.3), 10)
D = cellToDiscreteCell(  Cell([0.0, 0.0], x -> 1 - 0.2 * cos(2*x)), 10)
multiShapePlot(MutableLinkedList{DiscreteCell}(C,D), MutableLinkedList{String}("C", "ζ"))
R1 = rectangle(C)
R2 = rectangle(D)
R3 = rectangleOverlap(C, D)
plot!(xRect(R1), yRect(R1), label="R_C", seriestype=:shape, opacity=0.1)
plot!(xRect(R2), yRect(R2), label="R_ξ", seriestype=:shape, opacity=0.1)
plot!(xRect(R3), yRect(R3), label="R_overlap", seriestype=:shape, opacity=0.2)
savefig("Roverlap")


# figure critical edges 

C = DiscreteCell([ 5.0, 4.0, 1.0, 1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 7.0, 8.0, 7.0, 1.0, 1.0, 4.0 ], [0.0, 1.0, 1.0, 3.0, 3.0, 4.0, 5.0, 5.0, -5.0, -5.0, -4.0, -3.0, -3.0,-1.0,  -1.0])
D = DiscreteCell([0.0, 2.0, 4.0], [4.1, -2.0, 4.0])
E = cellToDiscreteCell(peanutCell([2,5, 2.8], 2.0), 10)
F = DiscreteCell([0.0, 2.0, 4.0, 2.0], [3.0, 0.0, 3.0, 6.0])
G = DiscreteCell([-1.5, 1.5, 1.5, -1.5], [-5.5, -5.5, 5.5, 5.5])
R = Rectangle(0.5,6.0,-2.0,4.0)
multiShapePlot(MutableLinkedList{DiscreteCell}(C), MutableLinkedList{String}("Cell"))
plot!(xRect(R), yRect(R), label="Rectangle", seriestype=:shape, opacity=0.2)
#=
liste = findAllCriticalEdgeLists(C, R) activate BakeTheLists4.jl
for i = 1:length(liste)

    for e ∈ liste[i]
        scatter!([C.x[e.i]], [C.y[e.i]], color=i+2, label=false)
    end 

end 
plot!()

savefig("finding_criticalEdges")
=#

### figure at the end of chapter 2 to illustrate that the discrete overlaps are getting computed correctly 
C1 = C
#C2 = moveC(G, 0.5, 0.0)
#C2 = moveC(G, 3.0, 0.0)
C2 = moveC(G, 5.5, 0.0)


plot(C1.x, C1.y, seriestype=:shape, opacity = 0.2, aspect_ratio=:equal, 
     label = false, xlims=(-2.0,10.0), ylims=(-6.0,9.0), dpi=500, legend=:topright, size = (800, 800))
plot!(C2.x, C2.y, seriestype=:shape, opacity = 0.2, aspect_ratio=:equal, label = false)
overlaps = getOverlap(C1, C2) 
lastO = last(overlaps)
pop!(overlaps)
areaoverlap = areaPolygon(lastO.x, lastO.y)
for o ∈ overlaps

    areaoverlap += areaPolygon(o.x, o.y)
    plot!(o.x, o.y, seriestype=:shape, opacity = 0.2, aspect_ratio=:equal, label = false, color = 3)

end

lab = "overlapping area = " * string(areaoverlap)
lab = lab[1:23]
plot!(lastO.x, lastO.y, seriestype=:shape, opacity = 0.2, aspect_ratio=:equal, label = lab, color = 3, legendfont=font(10))
savefig("discreteOverlapEx3")
# ----

susO = overlaps[1]


xRes 
