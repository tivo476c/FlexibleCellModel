include("cell_functionalities.jl")
include("computeOverlap.jl")
# include("parameters.jl")

using DifferentialEquations, StochasticDiffEq, Distributions, DataStructures

# Printf, StochasticDiffEq, Distributions, DataStructures, LinearAlgebra, OrdinaryDiffEq, Plots, ColorSchemes


#-------------------------------------- FINAL ENERGY

function energies!(du, u, p, t)
    """
     Deterministic part of our SDE, includes area, edge, interior angle and overlap force. 
     
     Args:
        du ... change that gets applied to u: u += dt*du + brownian 
        u  ... current state of the cell system [c1.x, ..., cM.x, c1.y, ..., cM.y] (in R^(2MN))
        p = dt, D, A_d, E_d, I_d, where
            dt ... time step size (in R) 
            D  ... diffusion coefficient (in R) 
            A_d ... desired cell areas for all cells (in R)
            E_d ... desired edge lengths of a cell (in R^N)
            I_d ... desired interior angles of a cell (in R^N)
        t  ... current time 
    """
    _, _, A_d, E_d, I_d = p
    if N != 0
        res = zeros(2 * N * M)
    else
        res = zeros(2 * M)
    end

    if (forceScalings[1] != 0)
        res += forceScalings[1] * areaForce(u, A_d; k=2)
    end
    if (forceScalings[2] != 0)
        res += forceScalings[2] * edgeForce(u, E_d; k=2)
    end
    if (forceScalings[3] != 0)
        res += forceScalings[3] * interiorAngleForce(u, I_d)
    end
    if (forceScalings[4] != 0)
        res += forceScalings[4] * overlapForce(u)
    end

    res[1] = 0
    res[3:5] = [0, 0, 0]
    res[7:8] = [0, 0]

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
        du = zeros(2 * M * N, 2 * M)
        for i = 1:2*M
            lineIdx = N*(i-1)+1:i*N
            du[lineIdx, i] = sqrt(2 * D)
        end
    end

end

function nomotion!(du, u, p, t)
    du = zeros(2 * M * N)
end
#-------------------------------------- BOUNDARY CONDITION

function apply_BC(u, t, integrator)
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

function apply_BC_overlap(u, t, integrator)
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

        u_i = [u[i], u[i+NumberOfCells]]
        for j = i+1:NumberOfCells
            u_j = [u[j], u[j+NumberOfCells]]

            distance = norm(u_i - u_j)
            if distance < 2 * radius

                pushVec = (u_i - u_j) / distance * (2 * radius - distance)
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

#------------------- AREA Force -> add A_d ∈ R to initialisation 

function areaEnergyCell(c, A_d; k=2)
    """
    Returns area energy of a cell c. 
    """
    return areaForceFactor / k * abs(A_d - areaPolygon(c.x, c.y))^k
end

function areaGradientCell(c)
    """
    Returns for a cell c a vector of length (R^2)^N, that holds for each vertex its area gradient.  
    First N entries are for the x coords of the N vertices and the second N entries for the y coords. 
    """

    #A = areaPolygon( c.x, c.y ) 
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
function areaForceCell(c, A_d; k=2)
    """
    Returns for a cell c a vector of length (R^2)^N, that holds for each vertex the area force that gets applied in each time step of a simulation.  
    First N entries are for the x coords of the N vertices and the second N entries for the y coords. 
    """

    A = areaPolygon(c.x, c.y)
    factor = 0.5 * sign(A_d - A) * abs(A_d - A)^(k - 1)

    return factor * areaGradientCell(c)

end

# final area Force that can be used in foo 
function areaForce(u, A_d; k=2)
    """
    Returns for a cell system u a vector of length (R^2N)^M, that holds for each vertex the area force that gets applied in each time step of a simulation.  
    """
    res = zeros(2 * M * N)
    cells = getCellsFromU(u)
    for i = 1:M

        a = areaForceCell(cells[i], A_d; k=k)
        for j = 1:N
            res[N*(i-1)+j] = a[j]
            res[N*(i-1+M)+j] = a[j+N]
        end

    end
    return res

end



#------------------- EDGE Force -> add E_d ∈ R^(M*N) to initialisation 

# we need a vector of all wanted edge lengths E_d ∈ R^(M*N) in order to implement the edge Force 
# E_d[ (i-1)*N + 1 : i*N] are the edge lengths of cell i 

function computeEdgeLengths(c::DiscreteCell)
    """
    Returns vector R^N of all edge lengths of the given cell. 
    """
    res = zeros(N)
    for i = 1:N-1
        res[i] = norm([c.x[i], c.y[i]] - [c.x[i+1], c.y[i+1]], 2)
    end
    res[N] = norm([c.x[N], c.y[N]] - [c.x[1], c.y[1]], 2)
    return res
end

function edgeEnergyCell(c, E_d; k=2)
    """
    Returns edge energy of cell c. 
    """
    edgeLengths = computeEdgeLengths(c)
    res = 0
    for i = 1:N
        res += edgeForceFactor / k * abs(E_d[i] - edgeLengths[i])^k
    end
    return res
end

# computes all Edge Energies for a single cell 
function edgeForceCell(c, E_d; k=2)
    """
    Returns the edge force vectors for all cell vertices of c. 
    """
    res = zeros(2 * N)
    E = computeEdgeLengths(c)
    # print("edges: E = $E, desired edges: E_d = $E_d")
    # all remaining values 
    for i = 1:N

        if i == 1
            next = 2
            prev = N
        elseif i == N
            next = 1
            prev = N - 1
        else
            next = i + 1
            prev = i - 1
        end

        res[i] = sign(E_d[prev] - E[prev]) / E[prev] * abs(E_d[prev] - E[prev])^(k - 1) * (c.x[i] - c.x[prev])
        +sign(E_d[i] - E[i]) / E[i] * abs(E_d[i] - E[i])^(k - 1) * (c.x[i] - c.x[next])
        res[i+N] = sign(E_d[prev] - E[prev]) / E[prev] * abs(E_d[prev] - E[prev])^(k - 1) * (c.y[i] - c.y[prev])
        +sign(E_d[i] - E[i]) / E[i] * abs(E_d[i] - E[i])^(k - 1) * (c.y[i] - c.y[next])
    end

    return res

end

# computes all vertice forces for all cells 
function edgeForce(u, E_d; k=2)
    """
    Returns all edge force vectors for all vertices in the cell system u.  
    """
    res = zeros(2 * M * N)
    cells = getCellsFromU(u)
    for i = 1:M
        e = edgeForceCell(cells[i], E_d; k=k)
        for j = 1:N
            res[N*(i-1)+j] = e[j]
            res[N*(i-1+M)+j] = e[j+N]
        end
    end
    return res
end


#------------------- INTERIOR ANGLE Force -> add I_d ∈ R^(M*N) to initialisation 

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

function computeIntAngles(c)

    V = Vector{Vector{Float64}}(undef, N)
    for i = 1:N
        V[i] = vertex(c, i)
    end
    res = zeros(N)
    for i = 1:N
        if i == 1
            prev = N
            next = 2
        elseif i == N
            prev = N - 1
            next = 1
        else
            prev = i - 1
            next = i + 1
        end
        res[i] = intAngleMT(V[prev], V[i], V[next])
    end
    return res

end

function distAngles(a1, a2) 
    """
    Returns the distance between the two angles a1, a2 in [0, 2*pi) considering the periodicity of the interval (0 = 2*pi). 
    """
    d = mod(a1-a2, 2*pi)
    return minimum( [d, 2*pi - d] )

end

function angleEnergyCell(c, I_d; k=2)
    res = 0
    intAngles = computeIntAngles(c)
    for i = 1:N
        res += interiorAngleForceFactor / k * distAngles(I_d[i], intAngles[i])^k
    end
    return res
end

function interiorAngleForceCell_MT1(c, I_d)
    """
    given: 
        * 1 DF cell c 
        * list of desired cell interior angles for that cell I 
    """

    intAngles = computeIntAngles(c)
    res = zeros(2 * NumberOfCellWallPoints) # res[1:NumberOfCellWallPoints] for x coordinates, res[NumberOfCellWallPoints+1:2*NumberOfCellWallPoints] for y coordinates
    for i = 1:NumberOfCellWallPoints

        if i == 1
            prev = NumberOfCellWallPoints
        else
            prev = i - 1
        end

        if i == NumberOfCellWallPoints
            next = 1
        else
            next = i + 1
        end

        v_prev = vertex(c, prev)
        v_curr = vertex(c, i)
        v_next = vertex(c, next)

        # assign x dynamic for vertex k 
        # res[i] += sign(I_d[prev] - intAngles[prev]) * distAngles(I_d[prev], intAngles[prev])^0 * (-1.0 / norm(v_curr - v_prev, 2)^2 * (v_curr[2] - v_prev[2]))
        res[i] += sign(I_d[i] - intAngles[i]) * distAngles(I_d[i], intAngles[i])^0 * (1.0 / norm(v_curr - v_prev, 2)^2 * (v_prev[2] - v_curr[2]))
        res[i] += sign(I_d[i] - intAngles[i]) * distAngles(I_d[i], intAngles[i])^0 * (-1.0 / norm(v_curr - v_next, 2)^2 * (v_next[2] - v_curr[2]))
        # res[i] += sign(I_d[next] - intAngles[next]) * distAngles(I_d[next], intAngles[next])^0 * (1.0 / norm(v_curr - v_next, 2)^2 * (v_curr[2] - v_next[2]))

        # assign y dynamic for vertex k 
        # res[NumberOfCellWallPoints+i] += sign(I_d[prev] - intAngles[prev]) * distAngles(I_d[prev], intAngles[prev])^0 * (-1.0 / norm(v_curr - v_prev, 2)^2 * (v_prev[1] - v_curr[1]))
        res[NumberOfCellWallPoints+i] += sign(I_d[i] - intAngles[i]) * distAngles(I_d[i], intAngles[i])^0 * (1.0 / norm(v_curr - v_prev, 2)^2 * (v_curr[1] - v_prev[1]))
        res[NumberOfCellWallPoints+i] += sign(I_d[i] - intAngles[i]) * distAngles(I_d[i], intAngles[i])^0 * (-1.0 / norm(v_curr - v_next, 2)^2 * (v_curr[1] - v_next[1]))
        # res[NumberOfCellWallPoints+i] += sign(I_d[next] - intAngles[next]) * distAngles(I_d[next], intAngles[next])^0 * (1.0 / norm(v_curr - v_next, 2)^2 * (v_next[1] - v_curr[1]))

    end

    return res

end

function interiorAngleForce(u, I_d)

    res = zeros(2 * M * N)
    C = getCellsFromU(u)
    for i = 1:M
        a = interiorAngleForceCell_MT1(C[i], I_d)
        for j = 1:N
            res[N*(i-1)+j] = a[j]
            res[N*(i-1+M)+j] = a[j+N]
        end
    end

    return res ./ 30

end


#------------------- OVERLAP Force (the notorious)
# pairs the vertex indice j from the overlap with the according vertex indice i in c1(v1) / c2(v2)

function overlapEnergy(u; k=1) 
    """
    Computes bachelor overlap energy of the cell system.  
    """
    if overlapForceFactor == 0 
        return 0 
    end 

    C = solutionToCells(u) 
    overlapEnergy = 0 
    for i = 1:length(C) 
        for j = i+1:length(C) 

            overlaps = getOverlap(C[i], C[j])
            for o ∈ overlaps
                area = areaPolygon(o.x, o.y)
                overlapEnergy += overlapForceFactor/k * area ^ k 
            end 
            
        end 
    end 
end 
    
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

        gradO = areaGradientCell(o)
        for ind ∈ v1

            i, j = ind
            r1x[i] = -0.5 * 2 * gradO[j]
            r1y[i] = -0.5 * 2 * gradO[j+K]
            # r1x[i] = -0.5 * area * gradO[j]
            # r1y[i] = -0.5 * area * gradO[j+K]

        end
        for ind ∈ v2

            i, j = ind
            r2x[i] = -0.5 * 2 * gradO[j]
            r2y[i] = -0.5 * 2 * gradO[j+K]
            # r2x[i] = -0.5 * area * gradO[j]
            # r2y[i] = -0.5 * area * gradO[j+K]

        end
    end
    return r1x, r1y, r2x, r2y

end

function radiusBilliardOverlapForce(u)
    # currently the size of overlap is not considered 
    res = zeros(2 * M)

    for i = 1:M
        centreI = [u[i], u[i+M]]
        for j = i+1:M
            centreJ = [u[j], u[j+M]]

            distance = norm(centreI - centreJ)
            if distance < 2 * radius

                # pushVec = (centreI - centreJ) / distance * (2 * radius - distance)  # pushVec b
                # pushVec = (centreI - centreJ) / distance                            # pushVec c
                pushVec = (centreI - centreJ)                                       # pushVec d
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

