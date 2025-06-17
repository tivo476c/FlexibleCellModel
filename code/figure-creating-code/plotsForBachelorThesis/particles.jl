using DataStructures
using Plots
using LinearAlgebra
include("../../cell_functionalities.jl")


domain = (0.0,10.0)
P = MutableLinkedList{Vector{Float64}}()

for i = 1:100

    free = false 
    x = rand(Float64, 2) * 10

    while(!free)

        x = rand(Float64, 2) * 10
        free = true  
        for p ∈ P
            if( norm(p-x, 2) <= 0.3 )
                free = false 
                break 
            end 
            if( ! (0.2 < x[1] < 9.8) || !(0.2 < x[2] < 9.8))
                free = false 
                break 
            end 
        end 

    end 

    push!(P,x)

end 

scatter( [P[1][1]], [P[1][2]], color = 1, markersize=6, xlims= domain, ylims= domain, label=false, aspect_ratio=:equal, dpi=500)
for i = 2:100
    scatter!([P[i][1]], [P[i][2]], color = 1, markersize=6, xlims= domain, ylims= domain, label=false, aspect_ratio=:equal, dpi=500)
end 

#savefig("SphereModel")
plot!() 


r = 0.5
f(x) = 0.5*sin(2*x) + 0.5 + r
p(x) = -3*x^2 / π^2 + 6*x/π + 1
g(x) = f(x) * (sin(x)+2) * 1.5

c1 = Cell([4.5, 3.5], x-> f(x) * (sin(x)+2) * 1.5)

N = 200
z = zeros(2,N)

for i = 1:N
    angle = 2*π*i / N
    z[1, i] = c1.r(angle) * cos(angle) + c1.centre[1];
    z[2, i] = c1.r(angle) * sin(angle)+ c1.centre[2];
end

plot(z[1, :], z[2, :], seriestype=:shape, color = 1, xlims= domain, ylims= domain, label=false, aspect_ratio=:equal, dpi=500)

psi_c1(x,y) = psiFunction(c1,x,y) * 2 - 1

contour((0.0:0.01:10.0), (0.0:0.01:10.0), psi_c1,  xlims= domain, ylims= domain, dpi=500, aspect_ratio=:equal )
savefig("PhaseFieldModel")

d1 = DiscreteCell([3.0, 5.0, 7.0, 7.0,5.0, 3.0], [3.0, 1.0, 3.0, 7.0, 9.0, 7.0])
plot(d1.x, d1.y, seriestype=:shape, color = 1, xlims= domain, ylims= domain, label=false, aspect_ratio=:equal, dpi=500)
#savefig("VertexModel")




