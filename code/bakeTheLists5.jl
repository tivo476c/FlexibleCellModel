include("cell_functionalities.jl")

# structure to save an edge with its parametrization
mutable struct CriticalEdge   
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    critical::Bool      # tells whether the edge is in the critical area 
    f::Bool             # tells if the function exists (true), or c.x[j] = c.x[i] (false)
    m::Float64          # m = (c.y[j] - c.y[i]) / (c.x[j] - c.x[i])
    a::Float64          # a = (c.x[j]*c.y[i] - c.x[i]*c.y[j]) / (c.x[j] - c.x[i])
                        # such that f(t) = m*t + a (t ∈ [min(x_i, x_{i+1}), max(x_i, x_{i+1})]) is the parametrization of the path that is exactly the line segment
end

# structure to holt a DF cell with all its edges and the covering rectangle
mutable struct CriticalCell 
    x::Vector{Float64}
    y::Vector{Float64}
    edge::Vector{CriticalEdge}
    R::Rectangle
end 

# structure to hold an intersections coordinates and the indices on which it lays in cell 1 and cell 2 (DF cells)
mutable struct Intersection
    x::Float64  # x coordinate of intersection point
    y::Float64  # y coordinate of intersection point
    i::Int64    # indice of edge of cell 1
    j::Int64    # indice of edge of cell 2
end

import Base: ==
function ==(i1::Intersection, i2::Intersection)
    return i1.x == i2.x && i1.y == i2.y && i1.i == i2.i && i1.j == i2.j
end
function ==(e1::CriticalEdge, e2::CriticalEdge)
    return e1.xmin == e2.xmin && e1.xmax == e2.xmax && e1.ymin == e2.ymin && e1.ymax == e2.ymax
end


import Base: copy
function copy(c::CriticalCell)::CriticalCell
    return CriticalCell(deepcopy(c.x), deepcopy(c.y), deepcopy(c.edge), c.R)
end 

# returns Rectangle that covers the cell (discrete cell)
function rectangle(c1::CriticalCell)::Rectangle
   
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

# returns the overlapping rectangle for 2 DF cells 
function rectangleOverlap(c1::CriticalCell, c2::CriticalCell)
    a = max(c1.R.x0, c2.R.x0)
    b = min(c1.R.x1, c2.R.x1)
    c = max(c1.R.y0, c2.R.y0)
    d = min(c1.R.y1, c2.R.y1)
    if(a<=b && c <= d)
        return Rectangle(a,b,c,d)
    else 
        return Rectangle(0.0,0.0,0.0,0.0)
    end 

end 

# moves the DF C in direction (x,y)
function moveC(cell::CriticalCell, x::Float64, y::Float64, C2::CriticalCell)
    c = copy(cell)
    #x = 0.001875
    #y = 0.001233
    for i = 1:lastindex(c.x)
        c.x[i] += x 
        c.y[i] += y
    end 
    c.edge = computeEdges(DiscreteCell(c.x, c.y), rectangleOverlap(c, C2))
    return c 

    
end 

# moves the DF C a bit to the up right
function moveC(c::CriticalCell, C2::CriticalCell)

    x = 0.001875
    y = 0.001233
    
    xRes = c.x 
    yRes = c.y 
    for i = 1:length(c.x)
        xRes[i] += x 
        yRes[i] += y 
    end 
    rRes = Rectangle(c.R.x0+x, c.R.x1+x, c.R.y0+y, c.R.y1+y)
    eRes = computeEdges(DiscreteCell(c.x, c.y), rRes)

    return CriticalCell(xRes, yRes, eRes, rRes) 

end 

# computes the CriticalCell structure for a given DiscreteCell 
function buildCritical(c::DiscreteCell)::CriticalCell
    R = rectangle(c)
    return CriticalCell(c.x, c.y, computeEdges(c, R), R)
end 

# returns true iff p lays on the cell wall of c 
function pointOnCellWall(p::Vector{Float64}, c::CriticalCell)::Bool

    for e ∈ filter(x -> x.critical, c.edge) 
        if(e.f)
            if(e.xmin <= p[1] <= e.xmax && p[2] == param(e.m, p[1], e.a) )
                return true
            end 
        else 
            if(p[1] == e.xmin && e.ymin <= p[2] <= e.ymax)
                return true
            end 
        end 
    end 
    return false 

end 

# returns true iff e1 runs above an edge in c 
function badEdge(e1::CriticalEdge, c::CriticalCell)::Bool 

    for e2 ∈ filter(x -> x.critical, c.edge) 
        
        if(max(e1.xmin, e2.xmin) <= min(e1.xmax, e2.xmax))
            if( !e1.f && !e2.f)
                if(max(e1.ymin, e2.ymin) <= min(e1.ymax, e2.ymax))
                    return true 
                end 
            elseif( e1.f && e2.f )
                if(e1.m == e2.m && e1.a == e2.a)
                    return true
                end 
            end 
        end 


    end 

    return false 

end 

# returns true iff badEdge or pointOnCellWall is true for any point or edge 
function badCells(c1::CriticalCell, c2::CriticalCell)::Bool 
    # returns true if any wallpoint of one cells locates on the other cell wall 
    # bad if two edges run above each other or a wallpoint is located on the other cell wall 
    M = length(c1.x)
    checkedP = Set()
    for i = 1:M
        if(c1.edge[i].critical)
            if(i==M)
                next = 1
            else 
                next = i+1
            end 
            p1 = [c1.x[i], c1.y[i]]
            p2 = [c1.x[next], c1.y[next]]
            if(!in(p1, checkedP))
                push!(checkedP, p1)
                if(pointOnCellWall(p1, c2))
                    return true
                end 
            end 
            if(!in(p2, checkedP))
                push!(checkedP, p2)
                if(pointOnCellWall(p2, c2))
                    return true
                end 
            end
            if(badEdge(c1.edge[i], c2))
                return true
            end 
        end 
    end 

    N = length(c2.x)
    for j = 1:N
        if(c2.edge[j].critical)
            if(j==N)
                next = 1
            else 
                next = j+1
            end 
            p1 = [c2.x[j], c2.y[j]]
            p2 = [c2.x[next], c2.y[next]]
            if(!in(p1, checkedP))
                push!(checkedP, p1)
                if(pointOnCellWall(p1, c1))
                    return true
                end 
            end 
            if(!in(p2, checkedP))
                push!(checkedP, p2)
                if(pointOnCellWall(p2, c1))
                    return true
                end 
            end
        end 
    end 

    return false 

end 

# computes the criticalEdges (what will be every edge) for a DF cell 
function computeEdges(c::DiscreteCell, R::Rectangle)

    N = length(c.x)
    res = Vector{CriticalEdge}(undef, N)
    
    for i = 1:N

        if(i==N)          # next entry 
            j = 1
        else 
            j = i+1
        end

        if(c.x[i] <= c.x[j])
            x0 = c.x[i]
            y0 = c.y[i]
            x1 = c.x[j]
            y1 = c.y[j]
        else 
            x0 = c.x[j]
            y0 = c.y[j]
            x1 = c.x[i]
            y1 = c.y[i]
        end 
        critical = overlap(R, Rectangle(x0, x1, min(y0,y1), max(y0,y1)))
        f = c.x[i] != c.x[j]
        if(f)
            if(y0 == y1)
                m = 0.0
                a = y0
            else 
                m = (c.y[j] - c.y[i]) / (c.x[j] - c.x[i]) 
                a = (c.x[j]*c.y[i] - c.x[i]*c.y[j]) / (c.x[j] - c.x[i])
            end 
        else 
            m = 0
            a = 0 
        end 

        res[i] = CriticalEdge(x0, x1, min(y0,y1), max(y0,y1), critical, f, m, a)

    end 

    return res 

end 


param(m,t,a) = m*t + a
# returns the index where e can be found in c.edge 
function getEdgeIndex(e::CriticalEdge, c::CriticalCell)::Int64
    for i = 1:length(c.x) 
        if(e == c.edge[i])
            return i 
        end 
    end 
    return 0 
end 

# returns true iff x is located in R 
function pointInRectangle(x::Vector{Float64}, R::Rectangle)::Bool 
    
    if(length(x) == 2)
       
        return (R.x0 <= x[1] <= R.x1) && (R.y0 <= x[2] <= R.y1)
    else 
        return nothing
    end
    
end

# returns an intersection of e1 and e2, or nothing if there is no intersection
function findIntersection(e1::CriticalEdge, e2::CriticalEdge, c1::CriticalCell, c2::CriticalCell)

    a = max(e1.xmin, e2.xmin)
    b = min(e1.xmax, e2.xmax)
    c = max(e1.ymin, e2.ymin)
    d = min(e1.ymax, e2.ymax)
    if(a > b || c > d) 
        return nothing 
    end 


    if(!e1.f)
        if(!e2.f)
            i = getEdgeIndex(e1,c1)
            j = getEdgeIndex(e2,c2)
            return Intersection(e1.xmin, c, i, j)
        else 
            s = [e1.xmin, param(e2.m, e1.xmin, e2.a)]
            if(pointInRectangle(s, Rectangle(a,b,c,d)))
                i = getEdgeIndex(e1,c1)
                j = getEdgeIndex(e2,c2)
                return Intersection(s[1], s[2], i, j)
            else 
                return nothing
            end 
        end 
    end 

    if(!e2.f)
         
        s = [e2.xmin, param(e1.m, e2.xmin, e1.a)]
        if(pointInRectangle(s, Rectangle(a,b,c,d)))
            i = getEdgeIndex(e1,c1)
            j = getEdgeIndex(e2,c2)
            return Intersection(s[1], s[2], i, j)
        else 
            return nothing
        end 
    end 

    if(e1.m == 0)
        if(e2.m == 0)
            i = getEdgeIndex(e1,c1)
            j = getEdgeIndex(e2,c2)
            return Intersection(a, e1.ymin, i, j)
        else 

            t = (e1.ymin - e2.a)/e2.m
            if( a <= t <= b)
                i = getEdgeIndex(e1,c1)
                j = getEdgeIndex(e2,c2)
                return Intersection(t , e1.ymin, i, j)
            else 
                return nothing
            end 
        end 
    end 
    if(e2.m == 0)
         
        t = (e2.ymin - e1.a)/e1.m
        if( a <= t <= b)
            i = getEdgeIndex(e1,c1)
            j = getEdgeIndex(e2,c2)
            return Intersection(t , e2.ymin, i, j)
        else 
            return nothing
        end 

    end 

    if(e1.m == e2.m)
        if(e1.a == e2.a)
            i = getEdgeIndex(e1,c1)
            j = getEdgeIndex(e2,c2)
            return Intersection(a, param(e1.m, a, e1.a), i, j)
        else
            return nothing 
        end 
    end 

    t = (e2.a - e1.a) /(e1.m - e2.m) 
    if(a <= t <= b)
        i = getEdgeIndex(e1,c1)
        j = getEdgeIndex(e2,c2)
        return Intersection(t, param(e1.m, t, e1.a), i, j)
    else 
        return nothing
    end 

end 

# returns true iff e1 and e2 intersect 
function hasIntersection(e1::CriticalEdge, e2::CriticalEdge)::Bool 

    a = max(e1.xmin, e2.xmin)
    b = min(e1.xmax, e2.xmax)
    c = max(e1.ymin, e2.ymin)
    d = min(e1.ymax, e2.ymax)
    if(a > b || c > d) 
        return false  
    end 


    if(!e1.f)
        if(!e2.f)
            return true 
        else 
            s = [e1.xmin, param(e2.m, e1.xmin, e2.a)]
            if(pointInRectangle(s, Rectangle(a,b,c,d)))
                return true 
            else 
                return false
            end 
        end 
    end 

    if(!e2.f)
         
        s = [e2.xmin, param(e1.m, e2.xmin, e1.a)]
        if(pointInRectangle(s, Rectangle(a,b,c,d)))
            return true 
        else 
            return false 
        end 
    end 

    if(e1.m == 0)
        if(e2.m == 0)
            return true 
        else 

            t = (e1.ymin - e2.a)/e2.m
            if( a <= t <= b)
                return true 
            else 
                return false
            end 
        end 
    end 
    if(e2.m == 0)
         
        t = (e2.ymin - e1.a)/e1.m
        if( a <= t <= b)
            return true 
        else 
            return false 
        end 

    end 

    if(e1.m == e2.m)
        if(e1.a == e2.a)
            return true 
        else
            return false  
        end 
    end 

    t = (e2.a - e1.a) /(e1.m - e2.m) 
    if(a <= t <= b)
        return true 
    else 
        return false 
    end 

end

# collects all intersections of c1 and c2
function findAllIntersections(C1::CriticalCell, C2::CriticalCell)::MutableLinkedList{Intersection}

    res = MutableLinkedList{Intersection}()

    for e1 ∈ filter(x -> x.critical, C1.edge)
        for e2 ∈ filter(x -> x.critical, C2.edge)
            
            i = findIntersection(e1, e2, C1, C2)
            if(i !== nothing)
                push!(res, i)
            end 

        end 
    end 

    return res 

end 

# checks if c1 has an intersection on the same edge where i is located 
function hasIntersectionFirst1(i::Intersection, c1::CriticalCell, direction::Int64, I::MutableLinkedList{Intersection})
    # if there is another intersection on same edge in given direction, it will get returned (the closest)
    # returns 'nothing' otherwise  
    # if direction == 1 search on line from i to x_i+1 
    # if direction == -1 search on line from i to x_i 
    M = length(c1.x)
    Inew = MutableLinkedList{Intersection}()
    iP = [i.x, i.y]
    res = nothing 
    if(direction == 1 )
    
        if(i.i == M)
            next = 1 
        else 
            next = i.i + 1
        end 
        x = [c1.x[next], c1.y[next]] 
        d1 = norm( x - iP , 2 )
        dCurrent = d1 
        for inters ∈ filter(x -> x.i == i.i, I)
            d = norm( [inters.x, inters.y] - iP , 2 ) 
            if(0 < d < dCurrent && norm( x - [inters.x, inters.y] , 2) < d1)       
                dCurrent = d   
                res = inters 
            end 
        end 

    else 
    
        x = [c1.x[i.i], c1.y[i.i]]
        d1 = norm(iP - x, 2) 
        dCurrent = d1 
        for inters ∈ filter(x -> x.i == i.i, I)
            d = norm( [inters.x, inters.y] - iP , 2 ) 
            if(0< d < dCurrent && norm( x - [inters.x, inters.y] , 2) < d1)       
                dCurrent = d   
                res = inters 
            end 
        end 
    
    end 
    
    return res 
    
end

# checks if c1 has an intersection on the same edge where i is located 
function hasIntersectionFirst2(i::Intersection, c1::CriticalCell, direction::Int64, I::MutableLinkedList{Intersection})
    # if there is another intersection on same edge in given direction, it will get returned (the closest)
    # returns 'nothing' otherwise  
    # if direction == 1 search on line from i to x_i+1 
    # if direction == -1 search on line from i to x_i 
    M = length(c1.x)
    Inew = MutableLinkedList{Intersection}()
    iP = [i.x, i.y]
    res = nothing 
    if(direction == 1 )
    
        if(i.j == M)
            next = 1 
        else 
            next = i.j + 1
        end 
        x = [c1.x[next], c1.y[next]] 
        d1 = norm( x - iP , 2 )
        dCurrent = d1 
        for inters ∈ filter(x -> x.j == i.j, I)
            d = norm( [inters.x, inters.y] - iP , 2 ) 
            
            if(0< d < dCurrent && norm( x - [inters.x, inters.y] , 2) < d1 )        
                dCurrent = d   
                res = inters 
            end 
        end 

    else 
    
        x = [c1.x[i.j], c1.y[i.j]]
        d1 = norm(iP - x, 2) 
        dCurrent = d1 
        for inters ∈ filter(x -> x.j == i.j, I)
            d = norm( [inters.x, inters.y] - iP , 2 ) 
            if(0< d < dCurrent && norm( x - [inters.x, inters.y] , 2) < d1)       
                d1 = d   
                res = inters 
            end 
        end 
    
    end 
    
    return res 
    
end

# if there is an intersection on same edge in given direction, it will get returned (the closest)
# returns 'nothing' otherwise  
# if direction == 1 return closest to x_i
# if direction == -1 return closest to x_i+1 
function hasIntersection1(c1::CriticalCell, i::Int64, direction::Int64, I::MutableLinkedList{Intersection})

    if(direction == 1)
        x = [c1.x[i], c1.y[i]]
    else 
        if(i == length(c1.x))
            next = 1
        else 
            next = i + 1
        end 
        x = [c1.x[next], c1.y[next]]
    end 
    
    Inew = filter(x -> x.i == i, I)
    if(isempty(Inew))
        return nothing 
    else 

        res = Inew[1]
        d1 = norm( x - [res.x, res.y], 2)
        for ind = 2:lastindex(Inew)
            inters = Inew[ind]
            d = norm( x - [inters.x, inters.y], 2) 
            if(0 < d < d1)
                d1 = d
                res = inters 
            end 
        end 
        return res 

    end 

end

# same as hasIntersection1, but must be applied on other cell 
function hasIntersection2(c1::CriticalCell, i::Int64, direction::Int64, I::MutableLinkedList{Intersection})
    # if there is an intersection on same edge in given direction, it will get returned (the closest)
    # returns 'nothing' otherwise  
    # if direction == 1 return closest to x_i
    # if direction == -1 return closest to x_i+1 

    if(direction == 1)
        x = [c1.x[i], c1.y[i]]
    else 
        if(i == length(c1.x))
            next = 1
        else 
            next = i + 1
        end 
        x = [c1.x[next], c1.y[next]]
    end 
    
    Inew = filter(x -> x.j == i, I)
    if(isempty(Inew))
        return nothing 
    else 

        res = Inew[1]
        d1 = norm( x - [res.x, res.y], 2)
        for ind = 2:lastindex(Inew)
            inters = Inew[ind]
            d = norm( x - [inters.x, inters.y], 2) 
            if(0 < d < d1)
                d1 = d
                res = inters 
            end 
        end 
        return res 

    end 

end

# returns true iff p is located in c 
function isInsideDiscreteCell(p::Vector{Float64}, c::CriticalCell)::Bool
    # p point that gets tested
    # c cell that gets tested 

    Rc = rectangle(c)
    xMax = Rc.x1
    if(xMax < p[1])
        return false
    end 
    testEdge = CriticalEdge(p[1], Rc.x1, p[2], p[2], true, true, 0.0, p[2])
    
    # TODO: GO ON implement this function only using testEdge 
    numOfinters = 0     
    for edge ∈ filter(x -> x.ymin <= p[2] <= x.ymax, c.edge)
        if(hasIntersection(edge, testEdge))
            numOfinters += 1
        end
    end 


    if(iseven(numOfinters))
        return false 
    else 
        return true 
    end 
end 

# tells in which direction one has to traverse through c1, in order to continue on the edges that are part of the overlap
function findDirection(i::Intersection, c1::CriticalCell, ind::Int64, c2::CriticalCell)::Int64
    # returns 1 if you should go in positive list direction 
    # returns -1 if you should go in negative list direction
    
    # i is current intersection; c1 is the current cell; e is the current critical edge on the current cell, c2 is the OTHER cell 

    if(ind == length(c1.x))
        next = 1 
    else 
        next = ind + 1
    end 
    test = [i.x, i.y] 
    # now go slightly in positive list direction 

    p1 = [c1.x[ind], c1.y[ind]]
    p2 = [c1.x[next], c1.y[next]]
    direction = 0.0001 * (p2 - p1) / norm(p2 - p1, 2) 
    test += direction
    if(isInsideDiscreteCell(test, c2))
        return 1 
    else 
        return -1 
    end 
end 

# returns the path from currentInters to the next intersection that gets found 
function findPath1(c1::CriticalCell, c2::CriticalCell, currentInters::Intersection, I::MutableLinkedList{Intersection})
    # ich geh davon aus dass immer noch eine intersection kommen muss 


    # TODO: Go on fix this shit 
    
    M = length(c1.x)
    xRes = MutableLinkedList{Float64}()
    yRes = MutableLinkedList{Float64}()

    e = c1.edge[currentInters.i]
    listInd = currentInters.i

    if(listInd == 0)
        println("Given currentInters is not on the given edge list.")
        return xRes, yRes, nothing 
    end 
    
    direction = findDirection(currentInters, c1, listInd, c2)  
    if( direction == 1)
        
        #println("going in direction 1 on cell 1")

        iNew = hasIntersectionFirst1(currentInters, c1, 1, I)
        if(iNew !== nothing && iNew != currentInters)
            return xRes, yRes, iNew 
        end 
        
        listInd += 1
        if(listInd > M)
            listInd = 1 
        end 

        e = c1.edge[listInd] 
        #println(" c1.x[listInd] = ", c1.x[listInd])
        #println(" c1.y[listInd] = ", c1.y[listInd])
        push!(xRes, c1.x[listInd])
        push!(yRes, c1.y[listInd])

        counter = 0
        while(counter <= M)

            iNew = hasIntersection1(c1, listInd, 1, I)
            if(iNew === nothing)
                # no intersection on that edge 
                # move on listInd and push res 
                listInd += 1
                if(listInd > M)
                    listInd = 1 
                end 
                e = c1.edge[listInd] 
                #println(" c1.x[listInd] = ", c1.x[listInd])
                #println(" c1.y[listInd] = ", c1.y[listInd])
                push!(xRes, c1.x[listInd])
                push!(yRes, c1.y[listInd])
            else 
                # intersection found 
                # return current path + nearest intersection 
                return xRes, yRes, iNew
            end 

            counter += 1

        end         

    else 
        # now iterate through l1 in reversed order 
        #println("going in direction -1 on cell 1")

        iNew = hasIntersectionFirst1(currentInters, c1, -1, I)
        if(iNew !== nothing)
            return xRes, yRes, iNew 
        end 
        push!(xRes, c1.x[listInd])
        push!(yRes, c1.y[listInd])
        listInd -= 1
        if(listInd == 0)
            listInd = M
        end 
        
        counter = 0
        while(counter <= M)

            e = c1.edge[listInd]
            iNew = hasIntersection1(c1, listInd, -1, I)
            if(iNew === nothing)     
                push!(xRes, c1.x[listInd])
                push!(yRes, c1.y[listInd])
                listInd -= 1 
                if(listInd == 0)
                    listInd = M
                end 
                e = c1.edge[listInd]
            else  
                return xRes, yRes, iNew
            end

        end 
    end 

    println("In findPath1 konnte keine intersection auf dem verbleibenden weg gefunden werden")
    return xRes, yRes, nothing 

end 

# same as findPath1, but must be applied on other cell 
function findPath2(c2::CriticalCell, c1::CriticalCell, currentInters::Intersection, I::MutableLinkedList{Intersection})
    # ich geh davon aus dass immer noch eine intersection kommen muss 
    # actually muss die 1. zelle als 1. argument übergeben werden 

    M = length(c1.x)
    xRes = MutableLinkedList{Float64}()
    yRes = MutableLinkedList{Float64}()

    e = c1.edge[currentInters.j]
    listInd = currentInters.j

    if(listInd == 0)
        println("Given currentInters is not on the given edge list.")
        return xRes, yRes, nothing 
    end 
    
    direction = findDirection(currentInters, c1, listInd, c2)  
    if( direction == 1)
        
        #println("going in direction 1 on cell 2")
        iNew = hasIntersectionFirst2(currentInters, c1, 1, I)
        if(iNew !== nothing && iNew != currentInters)
            return xRes, yRes, iNew 
        end 
        
        listInd += 1
        if(listInd > M)
            listInd = 1 
        end 
        e = c1.edge[listInd] 
        push!(xRes, c1.x[listInd])
        push!(yRes, c1.y[listInd])

        counter = 0
        while(counter <= M)

            iNew = hasIntersection2(c1, listInd, 1, I)
            if(iNew === nothing)
                # no intersection on that edge 
                # move on listInd and push res 
                listInd += 1
                if(listInd > M)
                    listInd = 1 
                end 
                e = c1.edge[listInd] 
                push!(xRes, c1.x[listInd])
                push!(yRes, c1.y[listInd])
            else 
                # intersection found 
                # return current path + nearest intersection 
                return xRes, yRes, iNew
            end 

            counter += 1

        end         

    else 
        # now iterate through l1 in reversed order 
        # println("going in direction -1 on cell 2")

        iNew = hasIntersectionFirst2(currentInters, c1, -1, I)
        if(iNew !== nothing)
            return xRes, yRes, iNew 
        end 
        push!(xRes, c1.x[listInd])
        push!(yRes, c1.y[listInd])
        listInd -= 1
        if(listInd == 0)
            listInd = M
        end 
        
        counter = 0
        while(counter <= M)

            e = c1.edge[listInd]
            iNew = hasIntersection2(c1, listInd, -1, I)
            if(iNew === nothing)     
                push!(xRes, c1.x[listInd])
                push!(yRes, c1.y[listInd])
                listInd -= 1 
                if(listInd == 0)
                    listInd = M
                end 
                e = c1.edge[listInd]
            else  
                return xRes, yRes, iNew
            end

        end 
    end 

    println("In findPath2 konnte keine intersection auf dem verbleibenden weg gefunden werden")
    return xRes, yRes, nothing 

end 

# appends l2 at the end of l1 
function anhängen(l1, l2)
    # fick append!
    # hängt l2 hinten an l1 an 
   
    for l ∈ l2
        push!(l1, l)
    end 
  
end 

function forceOrientation(c)

    N = length(c.x)
    sum = 0
    for i = 1:N
        
        x1 = c.x[i]
        y1 = c.y[i]

        if(i==N)
            next = 1
        else
            next = i+1
        end 
        x2 = c.x[next]
        y2 = c.y[next] 

        sum += (x2 - x1)*(y2 + y1)
            
    end

    if(sum > 0) 
        return DiscreteCell( reverse(c.x), reverse(c.y) )
    else 
        return c
    end 
    
end

# constructs an overlap of c1 and c2, returns it as a DF cell as well as all used intersections 
function constructOverlap(c1::CriticalCell, c2::CriticalCell, I::MutableLinkedList{Intersection})

    currentI = first(I) 
    usedI = MutableLinkedList{Intersection}(currentI)
    xVal = MutableLinkedList{Float64}(currentI.x)
    yVal = MutableLinkedList{Float64}(currentI.y)

    #println("firstI =[", currentI.x, ", ", currentI.y, "]")

    for counter = 1:length(I)

        if(iseven(counter))
            xPath, yPath, iNew = findPath2(c1, c2, currentI, I)
        else 
            xPath, yPath, iNew = findPath1(c1, c2, currentI, I)
        end 

        #println("counter = ", counter)
        #println("iNew = [", iNew.x, ", ", iNew.y, "]")
        #println("xPath = ", xPath)
        #println("yPath = ", yPath)
        #println()
    
        anhängen(xVal, xPath)
        anhängen(yVal, yPath)


        if(iNew == first(I))
            # the cell is then closed 
            return forceOrientation(DiscreteCell(collect(xVal), collect(yVal))), usedI

        else 


            push!(usedI, iNew)

            currentI = iNew 
            push!(xVal, currentI.x)

            push!(yVal, currentI.y)

        end 
    end 

    #println(" yolo ")
    return forceOrientation(DiscreteCell(collect(xVal), collect(yVal))), usedI

end 

# returns the index where  inters can be found in I 
function giveIntersIndex(I::MutableLinkedList{Intersection}, inters::Intersection)
    for j = 1:length(I)
        if(I[j] == inters)
            return j
        end 
    end 
    return 0
end 

# deletes all intersections from Iused out of I 
function deleteIntersection(I::MutableLinkedList{Intersection}, Iused::MutableLinkedList{Intersection})

    for i ∈ Iused 
        j = giveIntersIndex(I, i)
        if(j != 0)
            delete!(I, j)
        end 
    end 

end

# returns a list of all overlaps of C11 and c2 
function getOverlap(C11::DiscreteCell, c2::DiscreteCell)
    
    overlaps = MutableLinkedList{DiscreteCell}()
    M = length(C11.x)
    N = length(c2.x) 
    R1 = rectangle(C11)
    R2 = rectangle(c2)
    if(!overlap(R1, R2))
        return overlaps
    end 

    R = Rectangle( max(R1.x0, R2.x0), min(R1.x1, R2.x1), max(R1.y0, R2.y0), min(R1.y1, R2.y1))
    c1 = CriticalCell(C11.x, C11.y, computeEdges(C11, R), R1)
    C2 = CriticalCell(c2.x, c2.y, computeEdges(c2, R), R2)

   

    C1 = copy(c1)
    while(badCells(C1, C2))
        C1 = moveC(C1, C2)
    end

    

    intersections = findAllIntersections(C1, C2)

    #println("length(intersections) = ", length(intersections))

    while(!isempty(intersections))

        

        o, Iused = constructOverlap(C1, C2, copy(intersections))
        if(o !== nothing )
            push!(overlaps, o)
            deleteIntersection(intersections, Iused)
        end 

    end 

    #=
    res = 0
    for o ∈ overlaps
        res += areaPolygon(o.x, o.y)
    end 
    =#

    return overlaps
end

# ___________________________________________________________________________________________________________________________

function multiShapePlot(c, labels)
    if(length(c) == 0) 
        return nothing
    end
    if(length(c) > length(labels))
        println("Use more labels in multiShapePlot")
        return nothing
    end

    plot!(c[1].x, c[1].y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, label=labels[1])

    for i = 2:lastindex(c)
        plot!(c[i].x, c[i].y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, label=labels[i])
    end
end 

#=
c = DiscreteCell([ 5.0, 4.0, 1.0, 1.0, 2.0, 3.0, 2.0, 0.0, 0.0, 7.0, 8.0, 7.0, 1.0, 1.0, 4.0 ], [0.0, 1.0, 1.0, 3.0, 3.0, 4.0, 5.0, 5.0, -5.0, -5.0, -4.0, -3.0, -3.0,-1.0,  -1.0])
d = DiscreteCell([0.5, 2.0, 4.0], [4.1, -2.0, 4.0])
e = cellToDiscreteCell(peanutCell([2,5, 2.8], 2.0), 10)
f = DiscreteCell([0.1, 2.0, 4.0, 2.0], [3.0, 0.0, 3.0, 6.0])
g = DiscreteCell([0.0, 5.0, 5.0, 0.0], [-3.0, -3.0, 5.0, 5.0])

# ___________________________________________________________________________________________________________________________


c1_0 = cellToDiscreteCell(circleCell([2.0 , 2.0], 2.0), 30)
c1_1 = cellToDiscreteCell(circleCell([10.0,2.0], 2.0), 30)

c2_0 = cellToDiscreteCell(circleCell([10.0, 4.0], 2.0), 30)
c2_1 = cellToDiscreteCell(circleCell([2.0 , 4.0], 2.0), 30)  


abstand = norm([0.001875, 0.001233], 2)

function findDEindices(vertex, u, M::Int64, N::Int64)
    
    for i = 1:M
        point = [u[i], u[i+M]] 
        if(point == vertex || norm(point - vertex, 2) <= abstand)
            return [i, i+M]
        end
    end
    
    for i = 1:N
        point = [u[2*M + i], u[2*M + i + N]]
        if(point == vertex || norm(point - vertex, 2) <= abstand)
            return [2*M + i, 2*M + i + N]
        end
    end
    
    return [0,0]

end 


function addOverlap(du, u, p, t)
    # u is the vector of [c1.x; c1.y; c2.x; c2.y](t)
    # p must tell how long c1 and c2 are such that 2*(p[1]+p[2]) = length(u) 
    M, N = p
    
    x = zeros(M)
    y = zeros(M)
    for i = 1:M
        x[i] = u[i]
        y[i] = u[N+i]
    end
    a = zeros(N)
    b = zeros(N)
    for i = 1:N
        a[i] = u[2*M + i]
        b[i] = u[2*M + i + N]
    end
    
    for i=1:M
        scale = norm( [c1_1.x[i], c1_1.y[i]] - [u[i], u[i+N]], 2)
        if(scale > 0.5)
            du[i] = (c1_1.x[i] - u[i]) / scale 
            du[i+N] = (c1_1.y[i] - u[i+N]) / scale 
        else 
            du[i] = (c1_1.x[i] - u[i])
            du[i+N] = (c1_1.y[i] - u[i+N])
        end 
    end 
    for i=1:N
        scale = norm( [c2_1.x[i], c2_1.y[i]] - [u[2*M + i], u[2*M + i+N]], 2)
        if(scale > 0.5)
            du[2*M + i] = (c2_1.x[i] - u[2*M + i]) / scale 
            du[2*M + i+N] = (c2_1.y[i] - u[2*M + i+N]) / scale 
        else 
            du[2*M + i] = (c2_1.x[i] - u[2*M + i])
            du[2*M + i+N] = (c2_1.y[i] - u[2*M + i+N])
        end 
    end 
    
    overlaps = getOverlap(DiscreteCell(x,y), DiscreteCell(a,b))
    for o ∈ overlaps
        
        #println("length x = ", length(o.x), "; length y = ", length(o.y))
        K = length(o.x)
        area = areaPolygon(o.x, o.y)
        for i = 1:K
            
            vertex = [o.x[i], o.y[i]]
            DEindices = findDEindices(vertex, u, M, N)
            if DEindices != [0,0] 
                
                if(i == 1)
                    prev = K
                    next = 2
                elseif(i == K)
                    prev = K-1
                    next = 1
                else 
                    prev = i-1
                    next = i+1
                end 
                
                #x indice
                du[DEindices[1]] += -0.5 * area * ( o.y[next] - o.y[prev] )
                du[DEindices[2]] += -0.5 * area * ( o.x[prev] - o.x[next] )         
            
            end 
        end 
    end
end 

u0 = [c1_0.x; c1_0.y; c2_0.x; c2_0.y]
p = [length(c1_0.x), length(c2_0.x)]
tspan = (0.0, 10.0)
tstep = 0.5
prob = ODEProblem(addOverlap,u0,tspan,p)
sol = solve(prob, Tsit5(), saveat=tstep)

println(":-)")


multiShapePlot(MutableLinkedList{DiscreteCell}(C1, C2), MutableLinkedList{String}("C1", "C2"))

=#




