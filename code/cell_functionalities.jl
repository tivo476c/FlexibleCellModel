using DataStructures
using Plots
using LinearAlgebra
using OrdinaryDiffEq

import Base: +
function +(v::Vector{Float64}, x::Float64)
    res = deepcopy(v) 
    for i = 1:length(v) 
        res[i] += x 
    end 
    return res 
end 

import Base: -
function -(v::Vector{Float64}, x::Float64)
    res = deepcopy(v) 
    for i = 1:length(v) 
        res[i] -= x 
    end 
    return res 
end 


[1.0, 1.2, 1.4] - 1.0


# structure for CRF cells 
mutable struct Cell
    centre::Vector{Float64}
    r::Function
end

# structure for DF cells 
mutable struct DiscreteCell
    x::Vector{Float64}
    y::Vector{Float64}
end

# structur to save rectangles
mutable struct Rectangle
    x0:: Float64
    x1:: Float64
    y0:: Float64
    y1:: Float64
end

# adds vertex at first index to DF cell 
function addFirst(c::DiscreteCell, p::Vector{Float64})::DiscreteCell
    return DiscreteCell( [p[1];C.x], [p[2]; C.y])
end 

# adds vertex at last index to DF cell 
function addLast(c::DiscreteCell, p::Vector{Float64})::DiscreteCell
    return DiscreteCell( [C.x; p[1]], [C.y; p[2]])
end 

function vertex(c::DiscreteCell, i::Int64)
    return [c.x[i], c.y[i]]
end 



# moves a DF cell in direction (x,y)
function moveC(cell::DiscreteCell, x::Float64, y::Float64)
    
    resX = copy(cell.x)
    resY = copy(cell.y)
    for i=1:length(cell.x)
        resX[i] += x 
        resY[i] += y 
    end 

    return DiscreteCell(resX, resY)

end 

# builds DF cell out of 2 vectors x and y 
function XYToDiscreteCell(x::Vector{Float64}, y::Vector{Float64}):: DiscreteCell
    if(length(x) == length(y))
        return DiscreteCell(x,y)
    else
       println("The lengths of x and y must match in XYToDiscreteCell.") 
       return nothing   
    end
end 

# turns CRF cell c to DF cell with N vertices 
function cellToDiscreteCell(c::Cell, N::Int64)::DiscreteCell

    x = Vector{Float64}(undef, N)
    y = Vector{Float64}(undef, N)
    for i = 1:N
        ϕ = 2 * pi / N * (i-1)
        x[i] = c.centre[1] + c.r(ϕ) * cos(ϕ) 
        y[i] = c.centre[2] + c.r(ϕ) * sin(ϕ) 
    end 
    return DiscreteCell(x,y)
    
end 



# returns the x and y coordinates of a DF cell 
function DiscreteCellToXY(c::DiscreteCell)
    if(length(c.x) == length(c.y))
        return c.x, c.y 
    else
        println("The lengths of x and y must match in DiscreteCellToXY.") 
        return nothing
    end
end

#erstellt Zelle mit konstantem Radius r > 0 und Mittelpunkt M
function circleCell(M:: Vector{Float64}, r:: Float64) :: Cell 
    return Cell(M, x -> r )
end 

# erstellt eine DF Zelle die das rechteck r darstellt mit N vertices
function rectangleCell( r::Rectangle, N::Int64 )

    if(N < 4 || isodd(N))
        return nothing 
    end 

    hx = r.x1 - r.x0 
    hy = r.y1 - r.y0 

    #amount of points on one horizontal side: 

    # N = 2*(nx + ny)
    #nx = round(Int64, hx * N * 0.5 / hy) 
    #ny = round(Int64, 0.5*(N - 2*nx)) 

    N4 = round(Int64, N/4)

    Δx = hx / (N4 )
    Δy = hy / (N4 )
    

    x1 = collect(r.x0:Δx:r.x1-Δx) 
    x2 = ones(N4) * r.x1 
    x3 = collect(r.x1:-Δx:r.x0+Δx) 
    x4 = ones(N4) * r.x0 

    y1 = ones(N4) * r.y0 
    y2 = collect(r.y0:Δy:r.y1-Δy)
    y3 = ones(N4) * r.y1
    y4 = collect( r.y1:-Δy:r.y0+Δy  )

    return DiscreteCell( [x1;x2;x3;x4], [y1;y2;y3;y4] )

end 



#erstellt unsere Peanutcell mit mindestradius r 
function peanutCell(M:: Vector{Float64}, r:: Float64) :: Cell 
    return Cell(M, x -> 0.5*sin(2*x) + 0.5 + r)
end

#Bestimmt den Punkt der Außenwand der Zelle, zugehörig zum Winkel phi und dessen Radius
function outerWallPoint(c::Cell, phi:: Float64):: Vector{Float64}
    x1 = c.centre[1] + cos(phi) * c.r(phi)
    x2 = c.centre[2] + sin(phi) * c.r(phi)
    return [x1, x2]
end 

#returns angle at given point (x,y)
function computeAngle(centre:: Vector{Float64}, x:: Vector{Float64}):: Float64
    
    if centre == x 
        return 0
    end
    
    dx = x[1] - centre[1] 
    dy = x[2] - centre[2] 
    
    if dx == 0
        if dy > 0 
            return 0.5*π
        else 
            return 1.5*π
        end
    end 
    
    if dy == 0
        if dx > 0
            return 0
        else
            return π
        end
    end
    
    if dy > 0
        
        if dx > 0
            return atan(dy/dx)
        else
            return atan(-dx/dy) + 0.5*π
        end
        
        
    else 
        
        if dx > 0
            return atan(-dx/dy) + 1.5*π
        else
            return atan(dy/dx) + π
        end
        
    end
end

# computes the psifunction value of cell c at point (x,y)
function psiFunction(c::Cell, x::Float64, y::Float64)::Float64
    s = 3
    ϕ = computeAngle(c.centre, [x,y])
    return 0.5*(1 + tanh(s*(c.r(ϕ) - norm( [x,y] - c.centre, 2 ))))
end 


#gibt Liste aller äußeren Punkte einer Zellwand mit der Winkelschrittweite dw (0 < dw < 2*pi) wieder
function outerWall(c:: Cell, dw:: Float64):: MutableLinkedList{Vector{Float64}} 

    res = MutableLinkedList{Vector{Float64}}(outerWallPoint(c,dw))

    #println(2*pi/dw)
    
    for w = 2*dw : dw : 2*pi+dw*0.99
        p = outerWallPoint(c,w)
        append!(res, p)
    end     
    
    return res;

end

# returns a x and y vector for wall points on a CRF cell 
function cellToXY(c:: Cell, dw:: Float64)

    N = round(Int64, 2*π / dw) 
    x = Vector{Float64}(undef, N)
    y = Vector{Float64}(undef, N)
    list = outerWall(c, dw)

    for i = 1:N 
        x[i] = list[i][1]
        y[i] = list[i][2]
    end

    return x, y

end

#bildet eine MutableLinkedList von Vectoren [x,y] auf ein Array ab 
function pointsToArray(P:: MutableLinkedList{Vector{Float64}}):: Array{Float64,2} 

    l = length(P)
    res = zeros(2,l)
    for i = 1:1:l
        res[:,i] = P[i]
    end 

    return res;
end 

#plots a cell
function plotCell(c:: Cell, dw:: Float64)

    c1 = pointsToArray( outerWall( c, dw ) )  
    plot(c1[1,:], c1[2,:], aspect_ratio=:equal)

end 

#adds a cell plot to an already existing plot
function plotCell!(c:: Cell, dw:: Float64)

    c1 = pointsToArray( outerWall( c, dw ) )  
    plot!(c1[1,:], c1[2,:], aspect_ratio=:equal)

end 

#creates a shape plot of a cell (CRF)
function shapePlot(c1:: Cell, lab::String ="")
    N = 200
    z = zeros(2,N)
    
    for i = 1:N
        angle = 2*π*i / N
        z[1, i] = c1.r(angle) * cos(angle) + c1.centre[1]
        z[2, i] = c1.r(angle) * sin(angle)+ c1.centre[2]
    end

    return plot(z[1, :], z[2, :], seriestype=:shape, aspect_ratio=:equal, label=lab)
     
end 

#creates a shape plot of a cell (discrete)
function shapePlot(c1:: DiscreteCell, lab::String ="")
    
    return plot(c1.x, c1.y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, label=lab, legendfont=font(14), marker =:circle)
     
end 


#creates a shape plot of multiple cells from a MutableLinkedList{Cell}
function multiShapePlot(c:: MutableLinkedList{DiscreteCell}, labels::MutableLinkedList{String}=nothing)

    if(length(c) == 0) 
        return nothing
    end
    if(length(c) > length(labels))
        println("Use more labels in multiShapePlot")
        return nothing
    end

    p = plot(c[1].x, c[1].y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, label=labels[1], legendfont=font(14))

    for i = 2:lastindex(c)
        plot!(p, c[i].x, c[i].y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, label=labels[i])
    end

    return p 

end 

#creates a shape plot of multiple cells from a MutableLinkedList{Cell}
function multiShapePlot(c:: MutableLinkedList{Cell}, labels::MutableLinkedList{String}=nothing)
       
    N = 1000
    l = length(c)
    z = zeros(2,N)
    
    if(l == 0) 
        return nothing
    end
    if(length(c) > length(labels))
        println("Use more labels in multiShapePlot")
        return nothing
    end
        
    for i = 1:N
        angle = 2*π*i / N
        z[1, i] = c[1].r(angle) * cos(angle) + c[1].centre[1]
        z[2, i] = c[1].r(angle) * sin(angle) + c[1].centre[2]
    end
        
    p = plot(z[1, :], z[2, :], seriestype=:shape, aspect_ratio=:equal, opacity=.25, label=labels[1], ticks = false, dpi=500, legendfont=font(14))
        
    for j = 2:l
        
        for i = 1:N
            angle = 2*π*i / N
            z[1, i] = c[j].r(angle) * cos(angle) + c[j].centre[1]
            z[2, i] = c[j].r(angle) * sin(angle) + c[j].centre[2]
        end
        plot!(p, z[1, :], z[2, :], seriestype=:shape, aspect_ratio=:equal, opacity=.25, label=labels[j], ticks = false, dpi=500)
         
    end 
       
    return p 
    
end 

# used for plot of rectangles
function xRect(r:: Rectangle) :: Vector{Float64}
    return [r.x0, r.x1, r.x1, r.x0]
end

# used for plot of rectangles
function yRect(r:: Rectangle) :: Vector{Float64}
    return [r.y0, r.y0, r.y1, r.y1]
end

# returns Rectangle that overlaps the cell (CRF cells)
function rectangle(c1:: Cell) :: Rectangle
   
    x0 = c1.centre[1]
    x1 = c1.centre[1]
    y0 = c1.centre[2]
    y1 = c1.centre[2]
    
    for i = 0:100
        
        ϕ = i * π * 0.02
        x = c1.r(ϕ) * cos(ϕ) + c1.centre[1]
        y = c1.r(ϕ) * sin(ϕ) + c1.centre[2]
        
        if x < x0
            
            x0 = x
            
        elseif x > x1
            
            x1 = x
            
        end 
        
        if y < y0 
            
            y0 = y 
            
        elseif y > y1
            
            y1 = y
        
        end
      
    end
    
    #return Float64[x_l, x_r, y_d, y_u]
    return Rectangle(x0, x1, y0, y1)
    
end


# returns Rectangle that overlaps the cell (DF cell)
function rectangle(c1:: DiscreteCell) :: Rectangle
   
    x0 = c1.x[1]
    x1 = c1.x[1]
    y0 = c1.y[1]
    y1 = c1.y[1]
    
    for i = 2 : length(c1.x)        
        
        x = c1.x[i]
        y = c1.y[i]
        
        if x < x0
            
            x0 = x
            
        elseif x > x1
            
            x1 = x
            
        end 
        
        if y < y0 
            
            y0 = y 
            
        elseif y > y1
            
            y1 = y
        
        end
      
    end
    
    #return Float64[x_l, x_r, y_d, y_u]
    return Rectangle(x0, x1, y0, y1)
    
end

# returns a Rectangle that overlaps both CRF cells 
function rectangleUnion(c1:: Cell, c2:: Cell) :: Rectangle
            
    r1 = rectangle(c1)
    r2 = rectangle(c2)
    
    return Rectangle(min(r1.x0, r2.x0), max(r1.x1, r2.x1), min(r1.y0, r2.y0), max(r1.y1, r2.y1))
        
end

# returns true iff the intersection of 2 rectangles isnt empty
function overlap(r1:: Rectangle, r2:: Rectangle) ::Bool

    return max(r1.x0, r2.x0) <= min(r1.x1, r2.x1) && max(r1.y0, r2.y0) <= min(r1.y1, r2.y1)

end 

# returns the intersecting rectangle of 2 CRF cells overlapping rectangle 
function rectangleOverlap(c1:: Cell, c2:: Cell) :: Rectangle
    
    r1 = rectangle(c1)
    r2 = rectangle(c2)
    
    if(r1.x0 > r2.x0) 
        return rectangleOverlap(c2, c1)
    end 
    
    if(overlap(r1, r2))
        
        #case2aII
        if(r1.y1 >= r2.y0)
            #println("2aII")
                    #println(r2.x0, min(r1.x1, r2.x1), r2.y0             , min(r1.y1, r2.y1))
            return Rectangle(r2.x0, min(r1.x1, r2.x1), max(r1.y0, r2.y0), min(r1.y1, r2.y1))
        end 
        
        #case2bII
        if(r1.y0 <= r2.y1)
            #println("2bII")
            return Rectangle(r2.x0, min(r1.x1, r2.x1), r1.y0, min(r1.y1, r2.y1))
        end 
        
    else 
        
        println("no overlap")
        
    end 

    return nothing 

end 

# returns the intersecting rectangle of 2 DF cells overlapping rectangle 
function rectangleOverlap(c1:: DiscreteCell, c2:: DiscreteCell)
    
    r1 = rectangle(c1)
    r2 = rectangle(c2)
    
    if(r1.x0 > r2.x0) 
        return rectangleOverlap(c2, c1)
    end 
    
    if(overlap(r1, r2))
        
        #case2aII
        if(r1.y1 >= r2.y0)
            #println("2aII")
                    #println(r2.x0, min(r1.x1, r2.x1), r2.y0             , min(r1.y1, r2.y1))
            return Rectangle(r2.x0, min(r1.x1, r2.x1), max(r1.y0, r2.y0), min(r1.y1, r2.y1))
        end 
        
        #case2bII
        if(r1.y0 <= r2.y1)
            #println("2bII")
            return Rectangle(r2.x0, min(r1.x1, r2.x1), r1.y0, min(r1.y1, r2.y1))
        end 
        
    else 
        
        println("no overlap")
        
    end 

    return nothing

end

# computes the integral of psi1*psi2 over the overlapping rectangle of 2 CRF cells using summed up midpoint rule 
function computeOverlap(c1::Cell, c2::Cell)::Float64
    
    if(! overlap(rectangle(c1), rectangle(c2)))
        return 0
    end  
    
    r = rectangleOverlap(c1, c2)
    
    nx = 30
    ny = 30
    
    hx = (r.x1 - r.x0) / nx
    hy = (r.y1 - r.y0) / ny
    sum = 0
    
    for i = 0:nx-1
        for j = 0:ny-1
            
            xM = [r.x0, r.y0] + 0.5*[hx, hy] + [i*hx, j*hy]
            sum += psiFunction(c1, xM[1], xM[2]) * psiFunction(c2, xM[1], xM[2])
            
        end
    end

    return hx*hy*sum
    
end 

#computes Volume of a CRF cell (integral of psi over cells rectangle)
function computeVolume(c::Cell)::Float64
    
    r = rectangle(c)
    nx = 20
    ny = 20
    
    hx = (r.x1 - r.x0) / nx
    hy = (r.y1 - r.y0) / ny
    sum = 0
    
    for i = 0:nx-1
        for j = 0:ny-1
            
            xM = [r.x0, r.y0] + 0.5*[hx, hy] + [i*hx, j*hy]
            sum += psiFunction(c, xM[1], xM[2])
            
        end
    end

    return hx*hy*sum
    
end 

#computes Volume of a discrete cell with all available radiuses r by computing all triangle areas between 2 neighboring vertices and the centre
function computeVolume(r:: Vector{Float64})::Float64
    
    ϕ_amount = length(r)
    summe = r[ϕ_amount] * r[1]
    
    for i = 2: ϕ_amount
        summe += r[i-1] * r[i] 
    end
    
    return 0.5*summe*sin(2*π/ϕ_amount)

end 


r_0(ϕ) = 0.5*sin(2*ϕ) + 1
r_1(ϕ) = 2
# definition of ODE #1, just the cell reforming to its base state 
function cellStep(du, u, p, t)
    for i = 1:lastindex(p)
        du[i] = r_1(p[i]) - u[i]
    end 
end


v_0 = computeVolume(Cell([0.0, 0.0], r_0))
v_1 = computeVolume(Cell([0.0, 0.0], r_1))
# definition of ODE #2, faster reforming for different volumina 
function cellStepVolume1(du, u, p, t)
    
    v = computeVolume(u)
    for i = 1:lastindex(p)
        du[i] = (r_1(p[i]) - u[i])*abs(v_1 - computeVolume(u))
    end 
end

#=
# HOW TO USE ODE SOLVER:
ϕ_amount = 50s
ϕ_set = zeros(ϕ_amount)
for i = 1:ϕ_amount
    ϕ_set[i] = 2*π*i / ϕ_amount
end 
u0 = r_0.(ϕ_set)
tspan = (0.0,0.5)
tstep = 0.01
prob = ODEProblem(cellStepVolume1, u0, tspan, ϕ_set)
sol = solve(prob, Tsit5(), saveat=tstep)

#println(v_0)
#println(v_1)
#plot(sol)
# sol[i] gives vector of radius evolution at angle ϕ_i 
# sol[i][j] returns the radius of angle ϕ_i at time t_j



#how to make a animation out of it 
@userplot CellPlot
@recipe function f(cp::CellPlot)
    
    t, centre = cp.args 
    label --> (t-1)/100 
    
    z = zeros(2,ϕ_amount)
    for i = 1:ϕ_amount
        z[1, i] = sol[t][i] * cos(ϕ_set[i]) + centre[1];
        z[2, i] = sol[t][i] * sin(ϕ_set[i]) + centre[2];
    end

    z[1,:], z[2,:]

end


#t_amount = (tspan[2] - tspan[1]) : tstep
centre = [0.0, 0.0]
times = 1:length(sol)

anim = @animate for t ∈ times

    cellplot(t, centre, seriestype=:shape, aspect_ratio=:equal)

end

gif(anim, "anim_fps4.gif", fps = 3)

=#

# used for areaPolygon
function getInfo(x1, y1, x2, y2)
    return x1*y2 - y1*x2
end

# computes the area of a polygon using shoelace lemma
function areaPolygon(p :: Vector{Vector{Float64}})
    N = length(p)
    firstx, firsty = p[1]
    prevx, prevy = firstx, firsty
    res = 0
    
    for i = 2:N
        nextx, nexty = p[i]
        res += getInfo(prevx, prevy, nextx, nexty)
        prevx = nextx
        prevy = nexty
    end
    
    res += getInfo(prevx, prevy, firstx, firsty)
    
    return 0.5*abs(res)    
    
end 

# computes the area of a polygon using shoelace lemma
function areaPolygon(x::Vector{Float64}, y::Vector{Float64}):: Float64
    
    N = length(x)
    if(N != length(y)) 
        println("In areaPolygon x and y must have same dimensions.")
        return 0
    end
    
    firstx = x[1]
    firsty = y[1]
    prevx, prevy = firstx, firsty
    res = 0
    
    for i = 2:N
        nextx = x[i]
        nexty = y[i]
        res += getInfo(prevx, prevy, nextx, nexty)
        prevx = nextx
        prevy = nexty
    end
    
    res += getInfo(prevx, prevy, firstx, firsty)
    
    return 0.5*abs(res)    
    
end

# computes the area of a polygon using shoelace lemma
function areaPolygon(u::Vector{Float64}):: Float64
    
    if(mod(length(u), 2) != 0) 
        println("In areaPolygon the length of a single vector must be even so that u = [x;y] holds.")
        return 0
    end
    
    N = convert(Int64, 0.5 * length(u))
    x = zeros(N)
    y = zeros(N)
    for i = 1:N
        x[i] = u[i]
        y[i] = u[N+i]
    end
    
    return areaPolygon(x,y)
    
end

l(x,y) = norm(x-y)
l(x1, y1, x2, y2) = norm([x1-x2, y1-y2])



# TODO: go on 

# computes the sum for the vertex i that will be used in the ODE; x,y ∈ R^N
function edgeSumI(x::Vector{Float64}, y::Vector{Float64}, i::Int64)
    
    N = length(x)
    if(N != length(y))
        println("In edgeSum both vectors have to have the same length.")
        return 0
    end
    if( i<1 || i>N) 
        println("In edgeSum i has to be: 1 <= i <O N.")
        return 0
    end
    
    if(i == 1)
        return edgeForce(x, y, 1, N) + edgeForce(x, y, 1, 2)
    elseif (i == N)
        return edgeForce(x, y, N, N-1) + edgeForce(x, y, N, 1)
    else 
        return edgeForce(x, y, i, i-1) + edgeForce(x, y, i, i+1)
    end
        
end

# returns the interior angle at the vertex (x[i], y[i]) inside of a DF cell 
function interiorAngle(x::Vector{Float64}, y::Vector{Float64}, i::Int64)
    
    N = length(x)
    if(N != length(y))
        println("In interiorAngle x and y must have same lengths.")
        return 0
    end
    
    if(i < 1 || i > N) 
        println("In interiorAngle i must be a vertex index (be in [1,N]).")
        return 0 
    end 
    
    if(i == 1)
        
        #println("x[N]= ", x[N])
        #println("y[N]= ", y[N])
        #println("x[2]= ", x[2])
        #println("y[2]= ", y[2])
        
        a = l(x[N], y[N], x[2], y[2])        # l_N,2
        b = l(x[N], y[N], x[1], y[1])        # l_N,1
        c = l(x[1], y[1], x[2], y[2])        # l_1,2
        
        #println(a)
        #println(b)
        #println(c)
        
    elseif(i == N)
        
        a = l(x[N-1], y[N-1], x[1], y[1])        # l_N-1,1
        b = l(x[N], y[N], x[1], y[1])        # l_N,1
        c = l(x[N], y[N], x[N-1], y[N-1])    # l_N,N-1
        
    else
        
        a = l(x[i-1], y[i-1], x[i+1], y[i+1])  # l_i-1,i+1
        b = l(x[i], y[i], x[i-1], y[i-1])      # l_i,i-1
        c = l(x[i], y[i], x[i+1], y[i+1])      # l_i,i+1
        
    end 
    
    d = (b^2 + c^2 - a^2)/(2*b*c)
    
    while(d < -1)
        d = -1   
        #println("-")
    end
    
    while(d>1)
        d = 1
        #println("+")
    end
    
    return acos(d)
    
end

# returns the interior angle at the vertex (u[i], u[i+ lu/2]) inside of a DF cell 
function interiorAngle(u::Vector{Float64}, i::Int64)
    
    M = length(u)
    if(mod(M,2) != 0 )
        println("In interiorAngle(u,i) length(u) must be even.")
        return 0
    end
    
    N = convert(Int64, M/2)
    
    #println(N)
    
    x = zeros(N)
    y = zeros(N)

    for i = 1:N
        x[i] = u[i]
        y[i] = u[i+N]
    end
    
    return interiorAngle(x,y,i)
    
end

# returns the force that is applied on the x coordinate caused by the interiorAngle energy 
function d_xi_interiorAngle(x::Vector{Float64}, y::Vector{Float64}, i::Int64)
    
    h = 0.001 
    ei = zeros(length(x))
    ei[i] = 1
    return (interiorAngle(x+ei*h,y,i) - interiorAngle(x,y,i)) / h 
    
end 

# returns the force that is applied on the y coordinate caused by the interiorAngle energy 
function d_yi_interiorAngle(x::Vector{Float64}, y::Vector{Float64}, i::Int64)
    
    h = 0.001
    ei = zeros(length(x))
    ei[i] = 1
    return (interiorAngle(x,y+ei*h,i) - interiorAngle(x,y,i)) / h 
    
end 

# returns the force that is applied on the vertex (x[i], y[i]) caused by the interiorAngle energy 
function di_interiorAngle(x::Vector{Float64}, y::Vector{Float64}, i::Int64)
    return [d_xi_interiorAngle(x, y, i), d_yi_interiorAngle(x, y, i)] 
end

function solutionToXY(sol::Vector{Float64}) 
    """
    Returns for the solution of A SINGLE TIME STEP
        sol = [x1,x2,...,xM,y1,y2,...,yM],

    vectors res = X, Y such that for all cells i = 1,...,M 
        X[i] = x_i, Y[i] = y_i.      
    """
    if NumberOfCellWallPoints!=0
        X = Vector{Vector{Float64}}(undef, M)
        Y = Vector{Vector{Float64}}(undef, M)
        for i = 1:M 
            X[i] = sol[1+(i-1)*N : i*N]
            Y[i] = sol[1+(i-1)*N + N*M: i*N + N*M]
        end 
    else
        X = Vector{Float64}(undef, M)
        Y = Vector{Float64}(undef, M)
        for i = 1:M 
            X[i] = sol[i]
            Y[i] = sol[i+M]
        end 
    end 
    return X,Y
end 

function solutionToCells(sol::Vector{Float64}) 
    """
    Returns for the solution of A SINGLE TIME STEP
        sol = [x1,x2,...,xM,y1,y2,...,yM],

    a vector C = [c1, ..., cM] of all DiscreteCells,    
    such that for all cells i = 1,...,M 
       c[i].x = x_i, c[i].y = y_i.       
    """
    C = DiscreteCell[]
    if NumberOfCellWallPoints!=0
        for i = 1:M 
            x = sol[1+(i-1)*N : i*N]
            y = sol[1+(i-1)*N + N*M: i*N + N*M]
            c = DiscreteCell(x,y)
            push!(C,c)
        end 
    else
        for i = 1:M 
            x = sol[i]
            y = sol[i+M]
            c = DiscreteCell(x,y)
            push!(C,c)
        end 
    end 
    return C
end 

function getCentre(c::DiscreteCell)
    x = sum(c.x) / N 
    y = sum(c.y) / N 
    return [x,y]
end 

function getCellsFromU(u)
    """
    Returns vector [c1, c2, ..., cM] of type DiscreteCell for u from DifferentialEquations.jl 
    """

    res = DiscreteCell[]
    for i = 1:M
        c_i = DiscreteCell(u[N*(i-1)+1:N*i], u[N*(i-1+M)+1:N*(i+M)])
        push!(res, c_i)
    end 
    return res 
    
end 

function circleArea(radius, N)
    """
    Returns area of a regular polygon with N vertices and each vertex having a distance of radius to centre. 
    """
    return 0.5*N*radius^2*sin(2*pi/N) 
end 

function circleEdgeLengths(radius, N)
    """
    Returns edge lengths of a regular polygon with N vertices and each vertex having a distance of radius to centre. 
    """
    return 2*radius*sin(pi/N)*ones(N)
end 

function circleInteriorAngles(N)
    """
    Returns interior angles of a regular polygon with N vertices and each vertex having a distance of radius to centre. 
    """
    return (N-2)/N * pi *ones(N)
end 
