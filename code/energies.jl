include("cell_functionalities.jl")
include("computeOverlap.jl")
include("parameters.jl")

using DifferentialEquations, StochasticDiffEq, Distributions, DataStructures

# Printf, StochasticDiffEq, Distributions, DataStructures, LinearAlgebra, OrdinaryDiffEq, Plots, ColorSchemes


#-------------------------------------- FINAL ENERGY

function energies!(du, u, p, t)

    if N != 0
        res = zeros(2 * N * M)
    else
        res = zeros(2 * M)
    end

    # F1 is vector of all current edge lenghts, needed for edge Force and interior angle Force 
    # J1 is vector of all current interior angles, needed for interior angle Force 
    if forceScalings[1] != 0    # area force
        A1 = zeros(M)
        for i = 1:M
            c = DiscreteCell(u[N*(i-1)+1:N*i], u[N*(i-1+M)+1:N*(i+M)])
            polygon = Vector{Vector{Float64}, N}(undef) 
            for j = 1:N
                x, y = c.x[j], c.y[j]
            end 
            push!(polygon, [x,y])
            A1[i] = areaPolygon(polygon)
        end
    end

    if forceScalings[2] != 0    # edge force
        F1 = zeros(M * N)
        for i = 1:M
            c = DiscreteCell(u[N*(i-1)+1:N*i], u[N*(i-1+M)+1:N*(i+M)])
            F1[(i-1)*N+1:i*N] = computeEdgeLengths(c)
        end

    end
    if forceScalings[3] != 0    # interior angle force 
        J1 = zeros(M * N)
        for i = 1:M
            c = DiscreteCell(u[N*(i-1)+1:N*i], u[N*(i-1+M)+1:N*(i+M)])
            J1[(i-1)*N+1:i*N] = computeInteriorAngles(c)
        end

    end

    if (forceScalings[1] != 0)
        res += forceScalings[1] * areaForce(u, A1)
    end
    if (forceScalings[2] != 0)
        res += forceScalings[2] * edgeForce(u, E1, F1)
    end
    if (forceScalings[3] != 0)
        res += forceScalings[3] * interiorAngleForce(u, I1, J1)
    end
    if (forceScalings[4] != 0)
        res += forceScalings[4] * overlapForce(u)
    end

    for i = 1:length(du)
        du[i] = res[i]
    end

    return du

end
#-------------------------------------- BROWNIAN MOTION 

function brownian!(du, u, p, t)

    Δt, D = p
    # x = rand(Normal(0.0, 1.0), 2 * NumberOfCells)

    if NumberOfCellWallPoints == 0
        du .= sqrt(2 * D) 
    elseif NumberOfCellWallPoints != 0
        # TODO: TO BE TESTED
        du = zeros(2*M*N, 2*M)
        for i = 1:2*M
            lineIdx = N*(i-1)+1:i*N
            du[lineIdx,i] = sqrt(2 * D) 
        end
    end

end

function nomotion(du, u, p, t)
    du = zeros(2 * M * N)
end


#-------------------------------------- BOUNDARY CONDITION

function apply_BC(u,t, integrator)
    return minimum(u) < -domainL || maximum(u) > domainL
end 

function reflectiveBC!(integrator)
    u = integrator.u
    for i in eachindex(u)
        if u[i] < -domainL
            u[i] = -domainL + (-domainL - u[i])
        elseif u[i] > domainL
            u[i] = domainL - (u[i] - domainL)
        end
    end
end 

CallBack_reflectiveBC = DiscreteCallback(apply_BC, reflectiveBC!)

function apply_BC_overlap(u,t, integrator)
    return true
end 

function reflectiveBC_overlap!(integrator)
    u = integrator.u
    
    # reflective BC
    for i in eachindex(u)
        if u[i] < -domainL
            u[i] = -domainL + (-domainL - u[i])
        elseif u[i] > domainL
            u[i] = domainL - (u[i] - domainL)
        end
    end

    # Overlap check 
    for i = 1:NumberOfCells
        
        u_i = [u[i], u[i + NumberOfCells]]
        for j = i+1:NumberOfCells
            u_j = [u[j], u[j + NumberOfCells]]

            distance = norm(u_i - u_j)
            if distance < 2*radius 

                pushVec = (u_i - u_j) / distance * (2*radius - distance)
                u[i] += pushVec[1]
                u[i+NumberOfCells] += pushVec[2]
                u[j] -= pushVec[1]
                u[j+NumberOfCells] -= pushVec[2]

            end 
        end 
    end 
end 

CallBack_reflectiveBC_cellOverlap = DiscreteCallback(apply_BC_overlap, reflectiveBC_overlap!)




#-------------------------------------- SET Force FUNCTIONS

#------------------- AREA Force -> add A1 ∈ R^M to initialisation 

# we need a vector of all wanted cell volumes A1 ∈ R^M in order to implement the area Force 

function areaGradient(c, a1=0.0)

    N = length(c.x)
    #A = areaPolygon( c.x, c.y ) 
    #factor = a1 - A  
    res = zeros(2 * N)
    res[1] = c.y[2] - c.y[N]
    res[N+1] = c.x[N] - c.x[2]
    res[N] = c.y[1] - c.y[N-1]
    res[2*N] = c.x[N-1] - c.x[1]
    for i = 2:N-1

        res[i] = c.y[i+1] - c.y[i-1]
        res[i+N] = c.x[i-1] - c.x[i+1]

    end

    return res

end


# computes all vertice forces for a single cell  
function areaForceCell(c, a1)

    A = areaPolygon(c.x, c.y)
    factor = 0.5 * (a1 - A)
    res = zeros(2 * N)

    res[1] = c.y[2] - c.y[N]
    res[N+1] = c.x[N] - c.x[2]
    res[N] = c.y[1] - c.y[N-1]
    res[2*N] = c.x[N-1] - c.x[1]
    for i = 2:N-1

        res[i] = c.y[i+1] - c.y[i-1]
        res[i+N] = c.x[i-1] - c.x[i+1]

    end

    return factor * res

end

# final area Force that can be used in foo 
function areaForce(u, A1)

    res = zeros(2 * M * N)

    for i = 1:M

        c = DiscreteCell(u[N*(i-1)+1:N*i], u[N*(i-1+M)+1:N*(i+M)])
        a = areaForceCell(c, A1[i])
        for j = 1:N

            res[N*(i-1)+j] = a[j]
            res[N*(i-1+M)+j] = a[j+N]

        end


    end

    return 0.2 * res

end

#------------------- EDGE Force -> add E1 ∈ R^(M*N) to initialisation 

# we need a vector of all wanted edge lengths E1 ∈ R^(M*N) in order to implement the edge Force 
# E1[ (i-1)*N + 1 : i*N] are the edge lengths of cell i 

function computeEdgeLengths(c::DiscreteCell)

    N = length(c.x)
    res = zeros(N)

    for i = 1:N-1
        res[i] = norm([c.x[i], c.y[i]] - [c.x[i+1], c.y[i+1]], 2)
    end
    res[N] = norm([c.x[N], c.y[N]] - [c.x[1], c.y[1]], 2)

    return res
    # res[1] = length(vertex_1, vertex_2) | res[i] = length(vertex_i, vertex_i+1) | res[N] = length(vertex_N, vertex_1)

end



# computes all Edge Energies for a single cell 
function edgeForceCell(c, e1, f1)

    res = zeros(2 * N)
    #f1 = computeEdgeLengths(c) 

    # 1st and last x value 
    res[1] = (e1[1] / f1[1] - 1) * (c.x[1] - c.x[2]) + (e1[N] / f1[N] - 1) * (c.x[1] - c.x[N])
    res[N] = (e1[N] / f1[N] - 1) * (c.x[N] - c.x[1]) + (e1[N-1] / f1[N-1] - 1) * (c.x[N] - c.x[N-1])
    # 1st and last y value 
    res[N+1] = (e1[1] / f1[1] - 1) * (c.y[1] - c.y[2]) + (e1[N] / f1[N] - 1) * (c.y[1] - c.y[N])
    res[2*N] = (e1[N] / f1[N] - 1) * (c.y[N] - c.y[1]) + (e1[N-1] / f1[N-1] - 1) * (c.y[N] - c.y[N-1])

    # all remaining values 
    for i = 2:(N-1)

        res[i] = (e1[i] / f1[i] - 1) * (c.x[i] - c.x[i+1]) + (e1[i-1] / f1[i-1] - 1) * (c.x[i] - c.x[i-1])
        res[i+N] = (e1[i] / f1[i] - 1) * (c.y[i] - c.y[i+1]) + (e1[i-1] / f1[i-1] - 1) * (c.y[i] - c.y[i-1])

    end

    return res

end



# computes all vertice forces for all cells 
function edgeForce(u, E1, F1)

    res = zeros(2 * M * N)

    for i = 1:M

        c = DiscreteCell(u[N*(i-1)+1:N*i], u[N*(i-1+M)+1:N*(i+M)])
        edges = E1[(i-1)*N+1:i*N]
        f1 = F1[(i-1)*N+1:i*N]
        e = edgeForceCell(c, edges, f1)
        for j = 1:N

            res[N*(i-1)+j] = e[j]
            res[N*(i-1+M)+j] = e[j+N]

        end


    end

    return res


end


#------------------- INTERIOR ANGLE Force -> add I1 ∈ R^(M*N) to initialisation 

# we need a vector of all wanted interior angles I1 ∈ R^(M*N) in order to implement the edge Force 
# I1[ (i-1)*N + 1 : i*N] are the edge lengths of cell i 

#= old stuff 
function computeCenter(c::DiscreteCell)

    x = 0.0
    y = 0.0
    for i = 1:N 
        x += c.x[i]
        y += c.y[i]
    end 

    return [x,y] / N 

end 

function intAngle(x,y,z, centre) 

    d1 = norm(x-y, 2)
    d2 = norm(y-z, 2)
    d3 = norm(z-x, 2) 

    r1 = norm(centre-x, 2)
    r2 = norm(centre-y, 2)
    r3 = norm(centre-z, 2)

    d = (d1^2 + d2^2 - d3^2) / (2*d1*d2)

    if( d <= -1)
        return π 
    elseif(d >= 1)
        return 0.0 
    else 
        res = acos(d) 
    end 

    if( r2*1.01 < min(r1,r3) ) 
        return 2*π - res 
    else 
        return res 
    end 

end 
=#

# turn in vectors that span the according angle 

function arctan2(x, y)
    if x > 0
        return atan(y / x)
    elseif x < 0
        if y > 0
            return atan(y / x) + pi
        elseif y < 0
            return atan(y / x) - pi
        elseif y == 0
            return pi
        end
    elseif x == 0
        if y > 0
            return pi / 2.0
        elseif y < 0
            return -pi / 2.0
        elseif y == 0
            println("(0,0) cant be put in arctan2(x,y)")
        end
    end
end

function intAngleMT(v_prev, v_curr, v_next)

    v1 = v_prev - v_curr
    v2 = v_next - v_curr

    a1 = arctan2(v1[1], v1[2])
    a2 = arctan2(v2[1], v2[2])
    res = mod(a1 - a2, (2 * π))
    return res
end


function intAngle2(x, y, z)

    v1 = x - y
    v2 = z - y
    return mod(-atan(v1[1] * v2[2] - v1[2] * v2[1], dot(v1, v2)), (2 * π))
end

function intAngle2(v1, v2)
    return mod(-atan(v1[1] * v2[2] - v1[2] * v2[1], dot(v1, v2)), (2 * π))
end

function intAngle3(x, y, z)
    v1 = x - y
    v2 = z - y
    return mod(atan(v1[2], v1[1]) - atan(v2[2], v2[1]), (2 * π))
end

function computeInteriorAngles(c::DiscreteCell)

    V = Vector{Vector{Float64}}(undef, N)
    for i = 1:N
        V[i] = vertex(c, i)
    end

    res = zeros(N)
    res[1] = intAngle2(V[N], V[1], V[2])
    res[N] = intAngle2(V[N-1], V[N], V[1])
    for i = 2:N-1
        res[i] = intAngle2(V[i-1], V[i], V[i+1])
    end

    return res

end

function d_xi_interiorAngle(x, y, z)

    h = 0.01
    return (intAngle2(x, [y[1] + h, y[2]], z) - intAngle2(x, y, z)) / h

end

function d_yi_interiorAngle(x, y, z)

    h = 0.01
    return (intAngle2(x, [y[1], y[2] + h], z) - intAngle2(x, y, z)) / h

end

function d_xi_interiorAngle2(x, y, z)

    v1 = x - y
    v2 = z - y


    return v1[2] / (v1[1]^2 + v1[2]^2) - v2[2] / (v2[1]^2 + v2[2]^2)

end

function d_yi_interiorAngle2(x, y, z)

    v1 = x - y
    v2 = z - y

    return v2[1] / (v2[1]^2 + v2[2]^2) - v1[1] / (v1[1]^2 + v1[2]^2)
end



function interiorAngleForceCell(c, I, J)

    res = zeros(2 * N)
    res[1] = (I[1] - J[1]) * d_xi_interiorAngle(vertex(c, N), vertex(c, 1), vertex(c, 2))
    res[N] = (I[N] - J[N]) * d_xi_interiorAngle(vertex(c, N - 1), vertex(c, N), vertex(c, 1))
    res[N+1] = (I[1] - J[1]) * d_yi_interiorAngle(vertex(c, N), vertex(c, 1), vertex(c, 2))
    res[2*N] = (I[N] - J[N]) * d_yi_interiorAngle(vertex(c, N - 1), vertex(c, N), vertex(c, 1))

    for i = 2:N-1

        res[i] = (I[i] - J[i]) * d_xi_interiorAngle(vertex(c, i - 1), vertex(c, i), vertex(c, i + 1))
        res[i+N] = (I[i] - J[i]) * d_yi_interiorAngle(vertex(c, i - 1), vertex(c, i), vertex(c, i + 1))

    end

    return -res

end

function interiorAngleForceCell2(c, I, J)

    res = zeros(2 * N)

    res[1] = (I[1] - J[1]) * d_xi_interiorAngle2(vertex(c, N), vertex(c, 1), vertex(c, 2))
    res[N] = (I[N] - J[N]) * d_xi_interiorAngle2(vertex(c, N - 1), vertex(c, N), vertex(c, 1))
    res[N+1] = (I[1] - J[1]) * d_yi_interiorAngle2(vertex(c, N), vertex(c, 1), vertex(c, 2))
    res[2*N] = (I[N] - J[N]) * d_yi_interiorAngle2(vertex(c, N - 1), vertex(c, N), vertex(c, 1))

    for i = 2:N-1

        res[i] = (I[i] - J[i]) * d_xi_interiorAngle2(vertex(c, i - 1), vertex(c, i), vertex(c, i + 1))
        res[i+N] = (I[i] - J[i]) * d_yi_interiorAngle2(vertex(c, i - 1), vertex(c, i), vertex(c, i + 1))

    end

    return res

end

function interiorAngleForceCell_MT1(c, I, J)
    """
    given: 
        * 1 DF cell c 
        * list of desired cell interior angles for that cell I 
        * list of current cell interior angles for that cell J 
    """

    res = zeros(2 * NumberOfCellWallPoints) # res[1:NumberOfCellWallPoints] for x coordinates, res[NumberOfCellWallPoints+1:2*NumberOfCellWallPoints] for y coordinates
    for k = 1:NumberOfCellWallPoints

        currentIdx = k
        if currentIdx == 1
            prevIdx = NumberOfCellWallPoints
        else
            prevIdx = currentIdx - 1
        end

        if currentIdx == NumberOfCellWallPoints
            nextIdx = 1
        else
            nextIdx = currentIdx + 1
        end

        v_prev = vertex(c, prevIdx)
        v_curr = vertex(c, currentIdx)
        v_next = vertex(c, nextIdx)

        # assign x dynamic for vertex k 
        # res[k] -= (J[prevIdx]    - I[prevIdx]   ) * (-1.0/norm(v_curr - v_prev, 2)^2*(v_curr[2] - v_prev[2]) )                
        res[k] -= (J[currentIdx] - I[currentIdx]) * (1.0 / norm(v_curr - v_prev, 2)^2 * (v_prev[2] - v_curr[2]))
        res[k] -= (J[currentIdx] - I[currentIdx]) * (-1.0 / norm(v_curr - v_next, 2)^2 * (v_next[2] - v_curr[2]))
        # res[k] -= (J[nextIdx]    - I[nextIdx]   ) * ( 1.0/norm(v_curr - v_next, 2)^2*(v_curr[2] - v_next[2]) )

        # assign y dynamic for vertex k 
        # res[NumberOfCellWallPoints + k] -= (J[prevIdx] - I[prevIdx])*      (- 1.0/norm(v_curr - v_prev, 2)^2*(v_prev[1] - v_curr[1]) )                
        res[NumberOfCellWallPoints+k] -= (J[currentIdx] - I[currentIdx]) * (1.0 / norm(v_curr - v_prev, 2)^2 * (v_curr[1] - v_prev[1]))
        res[NumberOfCellWallPoints+k] -= (J[currentIdx] - I[currentIdx]) * (-1.0 / norm(v_curr - v_next, 2)^2 * (v_curr[1] - v_next[1]))
        # res[NumberOfCellWallPoints + k] -= (J[nextIdx] - I[nextIdx])*      (  1.0/norm(v_curr - v_next, 2)^2*(v_next[1] - v_curr[1]) )

    end

    return res ./ 20

end

function interiorAngleForce(u, I1, J1)

    res = zeros(2 * M * N)
    for i = 1:M

        c = DiscreteCell(u[N*(i-1)+1:N*i], u[N*(i-1+M)+1:N*(i+M)])

        i1 = I1[(i-1)*N+1:i*N]
        j1 = J1[(i-1)*N+1:i*N]
        # a = interiorAngleForceCell2(c, i1, j1)
        a = interiorAngleForceCell_MT1(c, i1, j1)
        for j = 1:N

            res[N*(i-1)+j] = a[j]
            res[N*(i-1+M)+j] = a[j+N]

        end


    end

    return res

end


#------------------- OVERLAP Force (the notorious)

# pairs the vertex indice j from the overlap with the according vertex indice i in c1(v1) / c2(v2)
function collectOverlapIndices(o, c1, c2)

    v1 = Set()
    v2 = Set()
    for j = 1:length(o.x)

        vert = [o.x[j], o.y[j]]

        for i = 1:N

            p = vertex(c1, i)
            if (norm(p - vert, 2) <= 0.00002)
                push!(v1, [i, j])
                break
            end

            q = vertex(c2, i)
            if (norm(q - vert, 2) <= 0.00002)
                push!(v2, [i, j])
                break
            end

        end

    end

    return v1, v2

end

function overlapForceCells(c1, c2, overlapForceType=overlapForceType)
    """
    Returns the force caused by the cell overlap between the 2 given DF cells. 
    Also selects which overlap force is selected: 
        overlapforce ∈ {"bachelorThesis", "billiard", "combination"}
    """

    if (overlapForceType == "bachelorThesis")
        return bachelorOverlapForceCells(c1, c2)
    elseif (overlapForceType == "billiard")
        return billiardOverlapForceCells(c1, c2)
    elseif (overlapForceType == "combination")
        return combinationOverlapForceCells(c1, c2)
    elseif (overlapForceType == "radiusBilliard")
        return radiusOverlapForceCells(c1, c2)
    else
        print("error: 3rd argument overlapForceType must be in {bachelorThesis, billiard, combination}")
        return []
    end
end

function bachelorOverlapForceCells(c1, c2)

    r1x = zeros(N)
    r1y = zeros(N)
    r2x = zeros(N)
    r2y = zeros(N)
    overlaps = getOverlap(c1, c2)
    for o ∈ overlaps
        K = length(o.x)
        area = areaPolygon(o.x, o.y)
        # collect all vertices, that are part of the overlap. these are the vertices which get a force applied
        v1, v2 = collectOverlapIndices(o, c1, c2)

        gradO = areaGradient(o)
        for ind ∈ v1

            i, j = ind
            r1x[i] = -0.5 * area * gradO[j]
            r1y[i] = -0.5 * area * gradO[j+K]

        end
        for ind ∈ v2

            i, j = ind
            r2x[i] = -0.5 * area * gradO[j]
            r2y[i] = -0.5 * area * gradO[j+K]

        end
    end
    return r1x, r1y, r2x, r2y

end

function radiusBilliardOverlapForce(u)
    # currently the size of overlap is not considered 
    res = zeros(2*M)

    for i = 1:M 
        centreI = [u[i], u[i+M]]
        for j = i+1:M 
            centreJ = [u[j], u[j+M]]
            
            distance = norm(centreI - centreJ)
            if distance < 2*radius 

                pushVec = (centreI - centreJ) / distance * (2*radius - distance)
                res[i] += pushVec[1]
                res[i+M] += pushVec[2]
                res[j] -= pushVec[1]
                res[j+M] -= pushVec[2]

            end 

        end 
    end 

    return res 
end

function computeCenter(c::DiscreteCell)

    x = sum(c.x)
    y = sum(c.y)
    return [x, y] / N

end

function combinationOverlapForceCells(c1, c2, scalingBachelor=0.5)
    # scalingBachelor must be in [0,1]
    r1xBach, r1yBach, r2xBach, r2yBach = bachelorOverlapForceCells(c1, c2)
    r1xBill, r1yBill, r2xBill, r2yBill = billiardOverlapForceCells(c1, c2)
    scalingBilliard = 1 - scalingBachelor

    r1x = scalingBachelor * r1xBach + scalingBilliard * r1xBill
    r1y = scalingBachelor * r1yBach + scalingBilliard * r1yBill
    r2x = scalingBachelor * r2xBach + scalingBilliard * r2xBill
    r2y = scalingBachelor * r2yBach + scalingBilliard * r2yBill

    #resBachelor = bachelorOverlapForceCells(c1,c2)
    #resBilliard = billiardOverlapForceCells(c1,c2)

    return r1x, r1y, r2x, r2y
end

function overlapForce(u)

    if N == 0 
        # we can just do radiusBilliard
        return radiusBilliardOverlapForce(u)
    end 

    cells = Vector{DiscreteCell}(undef, M)
    for i = 1:M
        cells[i] = DiscreteCell(u[N*(i-1)+1:N*i], u[N*(i-1+M)+1:N*(i+M)])
    end

    res = zeros(2 * M * N)

    for i = 1:M
        for j = (i+1):M

            r1x, r1y, r2x, r2y = overlapForceCells(cells[i], cells[j])
            res[N*(i-1)+1:N*i] += r1x
            res[N*(i-1+M)+1:N*(i+M)] += r1y
            res[N*(j-1)+1:N*j] += r2x
            res[N*(j-1+M)+1:N*(j+M)] += r2y

        end
    end

    return res

end

