include("cell_functionalities.jl")
using DataStructures
include("BakeTheLists4.jl")

C = DiscreteCell([ 5.0, 4.0, 1.0, 1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 7.0, 8.0, 7.0, 1.0, 1.0, 4.0 ], [0.0, 1.0, 1.0, 3.0, 3.0, 4.0, 5.0, 5.0, -5.0, -5.0, -4.0, -3.0, -3.0,-1.0,  -1.0])
R = Rectangle(1.5, 8.0, -5.0, 5.0)
liste = findAllCriticalEdgeLists(C, R)

shapePlot(C, "Cell")
plot!(xRect(R), yRect(R), seriestype=:shape, opacity=0.3, label="R_overlap")
for i = 1:length(liste)

    for e ∈ liste[i]
        scatter!([C.x[e.i]], [C.y[e.i]], color=i+2, label=false)
    end 

end 

plot!()

I = MutableLinkedList{Float64}(2.0)
J = MutableLinkedList{Float64}(5.0)
push!(I, J)

i = [1,2,5]
j = i - 1
lastindex(i)

C1 = DiscreteCell([ 5.0, 4.0, 1.0, 1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 7.0, 8.0, 7.0, 1.0, 1.0, 4.0 ], [0.0, 1.0, 1.0, 3.0, 3.0, 4.0, 5.0, 5.0, -5.0, -5.0, -4.0, -3.0, -3.0,-1.0,  -1.0])
plot(C1.x, C1.y, seriestype=:shape, aspect_ratio=:equal)
areaPolygon(C3.x, C3.y)

C2 = addX(C1, -5.0)
C3 = addY(C2, -3.0)
C4 = addY(C1, -6.0)