using LinearAlgebra

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
c1 = DiscreteCell([0.0012360679774997899, -0.0012639320225002098, -0.006263932022500209, -0.00876393202250021, -0.0062639320225002125, -0.0012639320225002137], [0.0, 0.004330127018922193, 0.004330127018922193, 6.123233995736766e-19, -0.004330127018922192, -0.004330127018922195])
c2 = DiscreteCell([0.008094059041422404, 0.0037639320225002102, -0.0005661949964219821, -0.0005661949964219839, 0.0037639320225002094, 0.008094059041422403], [0.0024999999999999996, 0.005, 0.002500000000000002, -0.0024999999999999988, -0.005, -0.0025000000000000022])

# first wrong config 
# c1 = DiscreteCell([0.0018788942990846757, 1.1131757133992795e14, -0.004027864045000419, -0.0065278640450004205, -0.004027864045000423, 1.1131757133992764e14], [1.626303258728257e-19, 6.426922977864271e13, 0.004330127018922193, 6.123233995736766e-19, -0.004330127018922192, -6.42692297786425e13])
# c2 = DiscreteCell([0.0058579910639226145, 2.792564706188567e12, -0.00125564214596432, -0.001255642145964321, -8.377694118565703e12, 0.005857991063922613], [0.0024999999999999996, -4.836863954542246e12, 0.0024192503572867674, -0.002419250357286765, -1.4510591863626748e13, -0.0025000000000000022])



plt = plot(xlims=domain, ylims=domain)
################### function 
k = 1
r1x = zeros(N)
r1y = zeros(N)
r2x = zeros(N)
r2y = zeros(N)
overlaps, vertexLists = getOverlap(c1, c2)
# each overlap is now a list of vertices 
for o ∈ overlaps

    K = length(o.x)
    area = areaPolygon(o.x, o.y)
    v1, v2 = collectOverlapIndices(o, c1, c2)           # collect all vertices, that are part of the overlap. these are the vertices which get a force applied

    gradO = areaGradientCell(o; NVertices=K)
    for ind ∈ v1

        i, j = ind
        r1x[i] = -0.5 * area^(k - 1) * gradO[j]
        r1y[i] = -0.5 * area^(k - 1) * gradO[j+K]

    end
    for ind ∈ v2

        i, j = ind
        r2x[i] = -0.5 * area^(k - 1) * gradO[j]
        r2y[i] = -0.5 * area^(k - 1) * gradO[j+K]
        println("")
        println("r2y[$i] = $(r2y[i])")
        println("")
    end
end

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
            u1, u2, outsideInd_u, insideInd_u = getInside_OutsideVertices(intersec.i, c1, vertexList)
            v1, v2, outsideInd_v, insideInd_v = getInside_OutsideVertices(intersec.j, c2, vertexList)
            scatter!(plt, [v1[1], u1[1]], [v1[2], u1[2]], label=false, color="blue")    # outside vertices
            scatter!(plt, [v2[1], u2[1]], [v2[2], u2[2]], label=false, color="red")     #  inside vertices
            I = [1 0; 0 1]
            #TODO:  
            f = cross2d(v1-u1, v2-v1)
            g = cross2d(u2-u1, v2-v1)
            t = f / g
            dtu1 = ([-(v2[2] - v1[2]), v2[1] - v1[1]]*g - f*[-(v2[2] - v1[2]),   v2[1] - v1[1] ]) / g^2
            dtu2 = (                                    - f*[  v2[2] - v1[2] , -(v2[1] - v1[1])]) / g^2
            dtv1 = ([ -u1[2] + v2[2] , u1[1] - v2[1]]*g - f*[  u2[2] - u1[2] , -(u2[1] - u1[1])]) / g^2
            dtv2 = ([-(v1[2] - u1[2]), v1[1] - u1[1]]*g - f*[-(u2[2] - u1[2]),   u2[1] - u1[1] ]) / g^2
            
            println("dtv1 = ([ -u1[2] + v2[2] , u1[1] - v2[1]]*g - f*[  u2[2] - u1[2] , -(u2[1] - u1[1])]) / g^2")
            println("f = $f")
            println("g = $g")
            println("u1 = $u1, u2 = $u2")
            println("v1 = $v1, v2 = $v2")
            dwu1 = (1 - t)*I + (u2-u1) * dtu1'         
            dwu2 = t*I + (u2-u1) * dtu2'                      
            dwv1 = (u2-u1) * dtv1'                      
            dwv2 = (u2-u1) * dtv2'                        

            du2t = [r1x[insideInd_u], r1y[insideInd_u]]           # the change thats already applied to u2 
            dv2t = [r2x[insideInd_v], r2y[insideInd_v]]           # the change thats already applied to v2 
            areaGradient_i = [areaGradientOverlap[indV], areaGradientOverlap[indV+K]]
            dOi = 0.5 * area^(k - 1) * areaGradient_i             # grad_intersection Overlap 

            R = -dOi - dwu2*du2t - dwv2*dv2t 
            
            println("dwv1 = (u2-u1) * dtv1'")
            println("(u2-u1)=$(u2-u1), dtv1' = $(dtv1')")
            println("dwv1 = $dwv1")
            println("c1 = $c1")
            println("c2 = $c2")
            # dwv1_inverted = inv(dwv1)
            # dwu1_inverted = inv(dwu1)

            # dv1t = 0.5 * dwv1_inverted * R 
            # du1t = 0.5 * dwu1_inverted * R 
            dv1t = 0.5 * pinv(dwv1) * R 
            du1t = 0.5 * pinv(dwu1) * R 

            r1x[outsideInd_u] += dv1t[1]
            r1y[outsideInd_u] += dv1t[2]
            r2x[outsideInd_v] += du1t[1]
            r2y[outsideInd_v] += du1t[2]
            println("")
            println("r1x[$outsideInd_u] = $(r1x[outsideInd_u])")
            println("r1y[$outsideInd_u] = $(r1y[outsideInd_u])")
            println("r2x[$outsideInd_v] = $(r2x[outsideInd_v])")
            println("r2y[$outsideInd_v] = $(r2y[outsideInd_v])")
            println("")
        end 

    end 
end

r1x = forceScalings[4] * r1x
r1y = forceScalings[4] * r1y
r2x = forceScalings[4] * r2x
r2y = forceScalings[4] * r2y

c1new = DiscreteCell(c1.x + r1x*timeStepSize, c1.y + r1y*timeStepSize)
c2new = DiscreteCell(c2.x + r2x*timeStepSize, c2.y + r2y*timeStepSize)

################### PLOTTING  
u = [c1.x; c2.x; c1.y; c2.y]
# u = [c1new.x; c2new.x; c1new.y; c2new.y]

X, Y = solutionToXY(u) 
plot!(plt, X[1], Y[1],
                seriestype=:shape,
                aspect_ratio=:equal,
                opacity=0.25,
                dpi=500,
                label=false,
                # xlims=domain,
                # ylims=domain,
                )


for i = 2:NumberOfCells
    plot!(plt, X[i], Y[i],
        seriestype=:shape,
        aspect_ratio=:equal,
        opacity=0.25,
        label=false,
    )
end

plot!(plt)
