include("../../cell_functionalities.jl")
include("../BakeTheLists4.jl")

C = DiscreteCell([ 5.0, 4.0, 1.0, 1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 7.0, 8.0, 7.0, 1.0, 1.0, 4.0 ], [0.0, 1.0, 1.0, 3.0, 3.0, 4.0, 5.0, 5.0, -5.0, -5.0, -4.0, -3.0, -3.0,-1.0,  -1.0])
D = DiscreteCell([8.0, 7.0, 5.0, 1.5, 2.5], [5.0, -4.5, -2.0, -4.5, 2.0])
R = rectangleOverlap(C,D)
liste = findAllCriticalEdgeLists(C, R)

multiShapePlot(MutableLinkedList{DiscreteCell}(C, D), MutableLinkedList{String}("C", "ζ"))
plot!(xRect(R), yRect(R), seriestype=:shape, opacity=0.3, label="R_overlap")
for i = 1:length(liste)

    for e ∈ liste[i]
        scatter!([C.x[e.i]], [C.y[e.i]], color=i+2, label=false)
    end 

end 

plot!()