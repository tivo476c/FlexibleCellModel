using LinearAlgebra
using LaTeXStrings
include("../parameters.jl")
include("../energies.jl")

################### PARAMETERS 
NumberOfCellWallPoints = 6                  # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 2                           # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells

################### CELL CONFIG 
# initial config
# c1 = cellToDiscreteCell(circleCell([-0.006, 0.0], radius), 6) 
# c2 = cellToDiscreteCell(circleCell([0.006, 0.0], radius), 6; rotation=pi/6.0) 

# t1 
c1 = DiscreteCell([0.0012360679774997899, -0.0012639320225002098, -0.006263932022500209, -0.00876393202250021, -0.0062639320225002125, -0.0012639320225002137] .+ 0.001 , [0.0, 0.004330127018922193, 0.004330127018922193, 6.123233995736766e-19, -0.004330127018922192, -0.004330127018922195])
c2 = DiscreteCell([0.008094059041422404, 0.0037639320225002102, -0.0005661949964219821, -0.0005661949964219839, 0.0037639320225002094, 0.008094059041422403] .- 0.001, [0.0024999999999999996, 0.005, 0.002500000000000002, -0.0024999999999999988, -0.005, -0.0025000000000000022])

# first wrong config 
# c1 = DiscreteCell([0.0018788942990846757, 1.1131757133992795e14, -0.004027864045000419, -0.0065278640450004205, -0.004027864045000423, 1.1131757133992764e14], [1.626303258728257e-19, 6.426922977864271e13, 0.004330127018922193, 6.123233995736766e-19, -0.004330127018922192, -6.42692297786425e13])
# c2 = DiscreteCell([0.0058579910639226145, 2.792564706188567e12, -0.00125564214596432, -0.001255642145964321, -8.377694118565703e12, 0.005857991063922613], [0.0024999999999999996, -4.836863954542246e12, 0.0024192503572867674, -0.002419250357286765, -1.4510591863626748e13, -0.0025000000000000022])



plt = plot(xlims=(-0.009, 0.009), ylims=(-0.006, 0.006), dpi=500)

################### PLOTTING  
u = [c1.x; c2.x; c1.y; c2.y]
# u = [c1new.x; c2new.x; c1new.y; c2new.y]

X, Y = solutionToXY(u)
plot!(plt, X[1], Y[1],
    label=L"C_1",
    seriestype=:shape,
    aspect_ratio=:equal,
    opacity=0.25,
    dpi=500,
    # xlims=domain,
    # ylims=domain,
)


for i = 2:NumberOfCells
    plot!(plt, X[i], Y[i],
        label=L"C_2",
        seriestype=:shape,
        aspect_ratio=:equal,
        opacity=0.25,
    )
end

################### function 
k = 1

overlaps, vertexLists = getOverlap(c1, c2)

# now add dynamic for vertices that are not part of the overlap cell 
for vertexListIndex in eachindex(vertexLists)
    vertexList = vertexLists[vertexListIndex]
    overlap = overlaps[vertexListIndex]
    K = length(overlap.x)
    area = areaPolygon(overlap.x, overlap.y)
    # indices_c1, indices_c2 = collectOverlapIndices(o, c1, c2)           # collect all vertices, that are part of the overlap. these are the vertices which get a force applied
    areaGradientOverlap = areaGradientCell(overlap; NVertices=K)

    for indV in eachindex(vertexList)
        intersec = vertexList[indV]
        if isa(intersec, Intersection)
            # TODO: test whether first edge index "i" is always for c1 and "j" for c2  
            # u1, u2 ∈ c1; v1, v2 ∈ c2
            # v1, u1 are the overlap OUTside vertices of c1, c2 
            # v2, u2 are the overlap  INside vertices of c1, c2 
            scatter!(plt, [intersec.x], [intersec.y], label=false, color=3)
            u1, u2, outsideInd_u, insideInd_u = getInside_OutsideVertices(intersec.i, c1, vertexList)
            v1, v2, outsideInd_v, insideInd_v = getInside_OutsideVertices(intersec.j, c2, vertexList)
            scatter!(plt, [v1[1], u1[1]], [v1[2], u1[2]], label=false, color="white")    # outside vertices
            scatter!(plt, [v2[1], u2[1]], [v2[2], u2[2]], label=false, color="black")     #  inside vertices
           
        end

    end
end

intersectionstring = L"\vec{w}_{%$(1)}"

intersections = [
    (L"\vec{w}_{1}", 0.001, 0.0033),
    (L"\vec{w}_{2}", 0.001, -0.0033),
    (L"\vec{v}_{1}^{1}", 0.003, 0.0),
    (L"\vec{v}_{2}^{1}", 0.00045, 0.0047),
    (L"\vec{v}_{6}^{1}", 0.00045, -0.0047),
    (L"\vec{v}_{2}^{2}", 0.0035, 0.0055),
    (L"\vec{v}_{3}^{2}", -0.0023, 0.0028),
    (L"\vec{v}_{4}^{2}", -0.0023, -0.0023),
    (L"\vec{v}_{5}^{2}", 0.0035, -0.0052),
]

for (s, px, py) in intersections
    # color = occursin("-", s) ? :red : :green
    annotate!(px, py, text(s, :black, 14))
end

plot!(plt)