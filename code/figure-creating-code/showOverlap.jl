include("../energies.jl")
# include("../parameters.jl")

using Plots

# SETUP: 2 cells overlapping
C = circleCell([-1.5,0.0], 2.0)
cDF1 = cellToDiscreteCell(C, N) 

A_d = ones(M) * areaPolygon(cDF1.x, cDF1.y) # ∈ R^N
E_d = ones(N*M)                             # ∈ (R^N)^M
I_d = ones(N*M)                             # ∈ (R^N)^M
e = computeEdgeLengths(cDF1)
ia = computeInteriorAngles(cDF1)

xNeedle1 = [0.01, 0.002, -0.002, -0.01, -0.002, 0.002]
yNeedle1 = [0.0, 0.002706329386826, 0.002706329386826, 0.0, -0.002706329386826, -0.002706329386826]
cNeedle1 = DiscreteCell(xNeedle1, yNeedle1)


a= 0.0008 
b = 0.00245
xNeedle2 = [0.01, 0, -0.01, -0.01, 0, 0.01]
yNeedle2 = [a,b,a,-a,-b,-a]
cNeedle2 = DiscreteCell(xNeedle2, yNeedle2)
# while true
#     global b 
#     global cNeedle2
#     global yNeedle2
#     b -= 0.00001
#     yNeedle2 = [a,b,a,-a,-b,-a]
#     cNeedle2 = DiscreteCell(xNeedle2, yNeedle2)
#     area = areaPolygon(xNeedle2, yNeedle2)
#     if area < 6.495190528383289e-5
#         println("b = $b")
#         break
#     end 
# end

# initial condition plot

println("area = $(areaPolygon(cNeedle2.x, cNeedle2.y))")
# plt = plot(cNeedle2.x, cNeedle2.y,
#            seriestype=:shape,
#            aspect_ratio=:equal,
#            dpi=200
#           )
