include("energies.jl")

# One has to set: M, N;  u0, A1, E1, I1, 


# -------------- animations for area Force 

M = 1
N = 6

C = cellToDiscreteCell( circleCell([0.,0.],3.), N )
u0 = [C.x; C.y]

#for i = 1:length(C.x)
#    println(vertex(C, i))
#end

A1 = [10.]

tspan = (0.0, 20.0)
Δt = 1 / 2^(8)
D = 1.5
p = [Δt, D]
prob_cell1 = SDEProblem( energies, nomotion, u0, tspan, p, noise=WienerProcess(0., 0.))        # how to correctly pass Δt and D in p ?
sol = solve(prob_cell1, EM(), dt=Δt)
domain = (-3.5,3.5)

animSDE = @animate for t ∈ 0:300

    time = t[1]*4+1
    x = sol[time][1:N]
    y = sol[time][ M*N+1 : N*(M+1)]
    #lab = string("time: ", time)
    area = areaPolygon(x,y)
    xlab = string("a(C(", time ,")) = ", round(area, digits=2))
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=domain, ylims = domain, xguidefontsize=16, xlabel = xlab )
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[time][ (i-1)*N+1  : i*N]
        y = sol[time][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 


    # plot arrows: 
    forces = energies(sol[time][:], sol[time][:], 0., 0.)
    c1 = DiscreteCell(sol[time][ 1:N], sol[time][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M

        for i = 1:N 

            c1 = DiscreteCell( sol[time][(j-1)*N+1 : j*N],  sol[time][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 
    
end

gif(animSDE, fps = 5)



# t ∈ {0, 50, 100, 1000}
times = 1001
x = sol[times][1:N]
y = sol[times][ M*N+1 : N*(M+1)]
#lab = string("times: ", times)
area = areaPolygon(x,y)
xlab = string("a(C(", times-1 ,")) = ", round(area, digits=2))
plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=500, 
         label = false, xlims=domain, ylims = domain, xguidefontsize=16, xlabel = xlab )
forces = energies(sol[times][:], sol[times][:], 0., 0.)
c1 = DiscreteCell(sol[times][ 1:N], sol[times][ M*N+1 : N*(M+1)])
x1 = MutableLinkedList{Float64}()
y1 = MutableLinkedList{Float64}()
GR.setarrowsize(1)
for j = 1:M

    for i = 1:N 

        c1 = DiscreteCell( sol[times][(j-1)*N+1 : j*N],  sol[times][(j-1+M)*N+1 : (j+M)*N] )
        cellForceX = forces[(j-1)*N+1 : j*N ]
        cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
        v = vertex(c1, i)
        push!(x1, v[1])
        push!(x1, cellForceX[i] + v[1])
        push!(y1, v[2])
        push!(y1, cellForceY[i] + v[2])
        push!(x1, NaN)
        push!(y1, NaN)

    end 
    plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

end 
p1000 = plot!()


savefig("areaEnergy1000_2")


# -------------- animations for edge Force 


M = 1
N = 4

C = rectangleCell(Rectangle(-3.,3.,-1.,1.), N)
u0 = [C.x; C.y]
E1 = [2., 6., 2., 6.]

tspan = (0.0, 20.0)
Δt = 1 / 2^(8)
D = 1.5
p = [Δt, D]
prob_cell1 = SDEProblem( energies, nomotion, u0, tspan, p, noise=WienerProcess(0., 0.))        # how to correctly pass Δt and D in p ?
sol = solve(prob_cell1, EM(), dt=Δt)
domain = (-3.5,3.5)

animSDE = @animate for t ∈ 0:300

    time = t[1]*4+1
    x = sol[time][1:N]
    y = sol[time][ M*N+1 : N*(M+1)]
    #lab = string("time: ", time)
    area = areaPolygon(x,y)
    x1 = [x[1], y[1]]
    x2 = [x[2], y[2]]
    x3 = [x[3], y[3]]
    xlab = string( "e1(C(", time ,")) = ", round(norm( x1-x2, 2), digits=1), "; e2(C(", time ,")) = ", round(norm( x3-x2, 1),digits=2)  )
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=domain, ylims = (-5.,5.), xguidefontsize=13, xlabel = xlab )
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[time][ (i-1)*N+1  : i*N]
        y = sol[time][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 


    # plot arrows: 
    forces = energies(sol[time][:], sol[time][:], 0., 0.)
    c1 = DiscreteCell(sol[time][ 1:N], sol[time][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M

        for i = 1:N 

            c1 = DiscreteCell( sol[time][(j-1)*N+1 : j*N],  sol[time][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 
    
end

gif(animSDE, fps = 5)

# t ∈ {0, 100, 200, 1000}
begin 
    times = 1001
    x = sol[times][1:N]
    y = sol[times][ M*N+1 : N*(M+1)]
    #lab = string("times: ", times)
    area = areaPolygon(x,y)
    x1 = [x[1], y[1]]
    x2 = [x[2], y[2]]
    x3 = [x[3], y[3]]
    xlab = string( "e1(C(", times-1 ,")) = ", round(norm( x1-x2, 2), digits=1), "; e2(C(", times-1 ,")) = ", round(norm( x3-x2, 1),digits=2)  )
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=500, 
            label = false, xlims=domain, ylims = (-5.,5.), xguidefontsize=13, xlabel = xlab )
    forces = energies(sol[times][:], sol[times][:], 0., 0.)
    c1 = DiscreteCell(sol[times][ 1:N], sol[times][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M

        for i = 1:N 

            c1 = DiscreteCell( sol[times][(j-1)*N+1 : j*N],  sol[times][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 
    p0 = plot!()
end 
savefig("edgeEnergy1000")


# ---------------- Interior angle energy 
N = 6
M = 1

C = DiscreteCell([0.0, 2.0, 1.0, 2.0, -0.0, -1.0], [-1.0, -1.0, 0.0, 1.0, 1.0, 0.0])
shapePlot(C)

C1 = DiscreteCell([0.0, 2.0, 3.0, 2.0, -0.0, 1.0], [-1.0, -1.0, 0.0, 1.0, 1.0, 0.0])
shapePlot(C1)

I1 = zeros(6)
I1[1] = intAngle2( vertex(C1, 6),vertex(C1, 1), vertex(C1, 2) )
I1[6] = intAngle2( vertex(C1, 5),vertex(C1, 6), vertex(C1, 1) )
for i = 2:5
    I1[i] = intAngle2( vertex(C1, i-1),vertex(C1, i), vertex(C1, i+1) )
end 

u0 = [C.x; C.y]

tspan = (0.0, 15.0)
Δt = 1 / 2^(8)
D = 1.5
p = [Δt, D]
prob_cell1 = SDEProblem( energies, nomotion, u0, tspan, p, noise=WienerProcess(0., 0.))        # how to correctly pass Δt and D in p ?
sol = solve(prob_cell1, EM(), dt=Δt)
domain = (-3.5,3.5)

animSDE = @animate for t ∈ 0: round(Int64, length(sol)/4 - 1 )  

    time = t[1]*4+1
    x = sol[time][1:N]
    y = sol[time][ M*N+1 : N*(M+1)]
    ct = DiscreteCell(x,y)
    #lab = string("time: ", time)
    area = areaPolygon(x,y)
    a1 = round(  intAngle2( vertex(ct, 2), vertex(ct, 3),vertex(ct, 4) ) / π * 180 , digits=1)
    a2 = round(  intAngle2( vertex(ct, 5), vertex(ct, 6),vertex(ct, 1) ) / π * 180 , digits=1)

    xlab = string( "ι3(C(", time ,")) = ", a1, "°; ι6(C(", time ,")) = ", a2, "°"  )
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=domain, ylims = domain, xguidefontsize=13, xlabel = xlab )
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[time][ (i-1)*N+1  : i*N]
        y = sol[time][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 


    # plot arrows: 
    forces = energies(sol[time][:], sol[time][:], 0., 0.)
    c1 = DiscreteCell(sol[time][ 1:N], sol[time][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M

        for i = 1:N 

            c1 = DiscreteCell( sol[time][(j-1)*N+1 : j*N],  sol[time][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 
    
end

gif(animSDE, fps = 10)


# t ∈ {0, 50, 100, 3000}
begin 
    times = 1
    x = sol[times][1:N]
    y = sol[times][ M*N+1 : N*(M+1)]

    ct = DiscreteCell(x,y)

    a1 = round(  intAngle2( vertex(ct, 2), vertex(ct, 3),vertex(ct, 4) ) / π * 180 , digits=1)
    a2 = round(  intAngle2( vertex(ct, 5), vertex(ct, 6),vertex(ct, 1) ) / π * 180 , digits=1)

    xlab = string( "ι3(C(", times-1 ,")) = ", a1, "°; ι6(C(", times-1 ,")) = ", a2, "°"  )
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
            label = false, xlims=(-1.5, 5.5), ylims = (-3., 3.), xguidefontsize=13, xlabel = xlab )
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[times][ (i-1)*N+1  : i*N]
        y = sol[times][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 


    # plot arrows: 
    forces = energies(sol[times][:], sol[times][:], 0., 0.)
    c1 = DiscreteCell(sol[times][ 1:N], sol[times][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M

        for i = 1:N 

            c1 = DiscreteCell( sol[times][(j-1)*N+1 : j*N],  sol[times][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 
    p0 = plot!()
end 

#savefig("angleEnergy3000")


# ------------------ OVERLAP PLOT 
C = Cell([3.75, 5.0], x -> 3.0*(1 - 0.3*cos(2*x)))

M = 2   # amount of cells 
N = 20  # amount of vertices per cell 

cDF = cellToDiscreteCell(C, N) 
shapePlot(cDF)
cDF2 = moveC(cDF, 2.5, 0.0)
u0 = [cDF.x; cDF2.x; cDF.y; cDF2.y]


tspan = (0.0, 15.0)
Δt = 1 / 2^(8)
D = 1.5
p = [Δt, D]
prob_cell1 = SDEProblem( energies, nomotion, u0, tspan, p, noise=WienerProcess(0., 0.))        # how to correctly pass Δt and D in p ?
sol = solve(prob_cell1, EM(), dt=Δt)
domain = (1.0, 9.0)

animSDE = @animate for t ∈ 0: round(Int64, length(sol)/4 - 1 )  

    time = t[1]*4+1
    x = sol[time][1:N]
    y = sol[time][ M*N+1 : N*(M+1)]
    ct = DiscreteCell(x,y)
    #lab = string("time: ", time)
    area = areaPolygon(x,y)
    a1 = round(  intAngle2( vertex(ct, 2), vertex(ct, 3),vertex(ct, 4) ) / π * 180 , digits=1)
    a2 = round(  intAngle2( vertex(ct, 5), vertex(ct, 6),vertex(ct, 1) ) / π * 180 , digits=1)

    xlab = string( "ι3(C(", time ,")) = ", a1, "°; ι6(C(", time ,")) = ", a2, "°"  )
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=domain, ylims = domain, xguidefontsize=13, xlabel = xlab )
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[time][ (i-1)*N+1  : i*N]
        y = sol[time][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 


    # plot arrows: 
    forces = energies(sol[time][:], sol[time][:], 0., 0.)
    c1 = DiscreteCell(sol[time][ 1:N], sol[time][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M

        for i = 1:N 

            c1 = DiscreteCell( sol[time][(j-1)*N+1 : j*N],  sol[time][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 
    
end

gif(animSDE, fps = 10)

# times in {0,20, 50,500}
begin  
    times = 1
    domain = (0., 10.)
    xdomain = (-2.,12.)
    x = sol[times][1:N]
    y = sol[times][ M*N+1 : N*(M+1)]

    c1 = DiscreteCell( sol[times][ 1  : N], sol[times][ (M)*N+1  : (1+M)*N] )
    c2 = DiscreteCell( sol[times][ (2-1)*N+1  : 2*N], sol[times][ (2-1+M)*N+1  : (2+M)*N] )
    o = getOverlap(c1, c2)[1]

    xlab = string( "a(D(", times-1, ")) = ", round(areaPolygon(o.x, o.y), digits=2) )
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
        label = false, xlims=xdomain, ylims = domain, xguidefontsize=15, xlabel = xlab )
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[times][ (i-1)*N+1  : i*N]
        y = sol[times][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 


    # plot arrows: 
    forces = energies(sol[times][:], sol[times][:], 0., 0.)
    c1 = DiscreteCell(sol[times][ 1:N], sol[times][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M
        x1 = MutableLinkedList{Float64}()
        y1 = MutableLinkedList{Float64}()
        for i = 1:N 

            c1 = DiscreteCell( sol[times][(j-1)*N+1 : j*N],  sol[times][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 
    plot!()
end 
#savefig("overlap0")


# ------------------ Combined Forces #1: THE SHAPE SAVER 

#C = Cell([0., 0.0], x -> 2.0*(1 - 0.3*cos(2*x)))
C_D = Cell([0., 0.0], x -> 4.0*(0.5*sin(2*x) + 0.7))


M = 1   # amount of cells 
N = 40  # amount of vertices per cell 

#C = cellToDiscreteCell(C, N)
C = rectangleCell(Rectangle(-3.,3.,-1.,1.), N)
C_D = cellToDiscreteCell(C_D, N)
areaPolygon(C_D.x, C_D.y)
multiShapePlot(MutableLinkedList{DiscreteCell}(C, C_D), MutableLinkedList{String}("cell0", "cellD"))
u0 = [C.x; C.y]
A1 = ones(M) * areaPolygon(C_D.x, C_D.y) # ∈ R^M
E1 = ones(N*M)              # ∈ (R^N)^M
I1 = ones(N*M)              # ∈ (R^N)^M
e = computeEdgeLengths(C_D)
ia = computeInteriorAngles(C_D)
for i = 1:M 
    E1[(i-1)*N+1 : i*N] = e
    I1[(i-1)*N+1 : i*N] = ia
end 

tspan = (0.0, 15.0)
Δt = 1 / 2^(8)
D = 1.5
p = [Δt, D]
prob_cell1 = SDEProblem( energies, nomotion, u0, tspan, p, noise=WienerProcess(0., 0.))        # how to correctly pass Δt and D in p ?
sol = solve(prob_cell1, EM(), dt=Δt)
domain = (-5.0, 5.0)

animSDE = @animate for t ∈ 0: round(Int64, 150)  

    time = t[1]*4+1
    x = sol[time][1:N]
    y = sol[time][ M*N+1 : N*(M+1)]
    ct = DiscreteCell(x,y)

    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=domain, ylims = domain)
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[time][ (i-1)*N+1  : i*N]
        y = sol[time][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 


    # plot arrows: 
    forces = energies(sol[time][:], sol[time][:], 0., 0.)
    c1 = DiscreteCell(sol[time][ 1:N], sol[time][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M

        for i = 1:N 

            c1 = DiscreteCell( sol[time][(j-1)*N+1 : j*N],  sol[time][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 
    
end

gif(animSDE, fps = 10)

t = length(sol)

begin
    # times ∈ {0, 20, 200, 3000}
    times = 3001
    x = sol[times][1:N]
    y = sol[times][ M*N+1 : N*(M+1)]
    ct = DiscreteCell(x,y)

    xlab = string("t = ", times-1)
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=domain, ylims = domain, xguidefontsize=15, xlabel = xlab )
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[times][ (i-1)*N+1  : i*N]
        y = sol[times][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false)

    end 


    # plot arrows: 
    forces = energies(sol[times][:], sol[times][:], 0., 0.)
    c1 = DiscreteCell(sol[times][ 1:N], sol[times][ M*N+1 : N*(M+1)])
    x1 = MutableLinkedList{Float64}()
    y1 = MutableLinkedList{Float64}()
    GR.setarrowsize(1)
    for j = 1:M

        for i = 1:N 

            c1 = DiscreteCell( sol[times][(j-1)*N+1 : j*N],  sol[times][(j-1+M)*N+1 : (j+M)*N] )
            cellForceX = forces[(j-1)*N+1 : j*N ]
            cellForceY = forces[(j-1+M)*N+1 : (j+M)*N ]
            v = vertex(c1, i)
            push!(x1, v[1])
            push!(x1, cellForceX[i] + v[1])
            push!(y1, v[2])
            push!(y1, cellForceY[i] + v[2])
            push!(x1, NaN)
            push!(y1, NaN)

        end 
        plot!(collect(x1),collect(y1), arrow=(:closed, 2.0), color = j, label=false)

    end 

    plot!()
end
#savefig("shape3000")

plot(C.x, C.y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=500, 
         label = false, xlims=domain, ylims = domain, xguidefontsize=15, xlabel = "C(0)" )
#savefig("shape_initial")

# ------------------ Combined Forces #2: THE GROßE AND GANZE 

C = Cell([1.0, 1.75], x -> .80*(1 - 0.3*cos(2*x)))

M = 9   # amount of cells 
N = 20    # amount of vertices per cell 

cDF = cellToDiscreteCell(C, N) 
u0 = zeros(2*M*N) 

for i = 0:2
    for j = 0:2

        c = moveC(cDF, j*1.0, i*1.5)
        u0[ i*N*3 + N*j + 1 : i*N*3 + N*j + N] = c.x
        u0[ i*N*3 + N*j + 1 + M*N : i*N*3 + N*j + N + M*N] = c.y

    end 
end 

A1 = ones(M) * areaPolygon(cDF.x, cDF.y) # ∈ R^M
E1 = ones(N*M)              # ∈ (R^N)^M
I1 = ones(N*M)              # ∈ (R^N)^M
e = computeEdgeLengths(cDF)
ia = computeInteriorAngles(cDF)
for i = 1:M 
    E1[(i-1)*N+1 : i*N] = e
    I1[(i-1)*N+1 : i*N] = ia
end 



tspan = (0.0, 30.0)
Δt = 1 / 2^(6)
D = 2.0
p = [Δt, D]
prob_cell1 = SDEProblem( energies, brownian, u0, tspan, p, noise=WienerProcess(0., 0.))        
steps = 4*Δt
sol = solve(prob_cell1, EM(), dt=Δt, saveat=steps)
domain = (0.0, 10.0)

xdomain = (-1.,5.0)
ydomain = (-1., 7.0)
titte = "Areax50; Edgex1.5; Anglex1; Overlapx80"
animSDE = @animate for t ∈ 1:round(Int64, (length(sol)-4)/2)

    time = t[1]*2
    x = sol[time][1:N]
    y = sol[time][ M*N+1 : N*(M+1)]

    xlab = string("t = ", time)
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=xdomain, ylims = ydomain, xguidefontsize=15, xlabel = xlab, title = titte)
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[time][ (i-1)*N+1  : i*N]
        y = sol[time][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false )

    end 

end


# animSDE: scaling factor 10 for overlap 

gif(animSDE, "firstInteraction_Areax50_Edgex1.5_Anglex1_Overlapx80.mp4", fps = 15)

begin
    # t ∈ {1, ,2001, 3501 }
    times = 501
    x = sol[times][1:N]
    y = sol[times][ M*N+1 : N*(M+1)]

    xlab = string("t = ", 2*(times-1))
    plot(x, y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, 
         label = false, xlims=xdomain, ylims = ydomain, xguidefontsize=15, xlabel = xlab)
    #xlabel!(xlab, fontsize=1)

    for i = 2:M 

        x = sol[times][ (i-1)*N+1  : i*N]
        y = sol[times][ (i-1+M)*N+1  : (i+M)*N]
        plot!(x,y, seriestype=:shape, aspect_ratio=:equal, opacity=.25, dpi=300, label = false )

    end
    plot!() 
end

savefig("overlap_rescaled_1000")