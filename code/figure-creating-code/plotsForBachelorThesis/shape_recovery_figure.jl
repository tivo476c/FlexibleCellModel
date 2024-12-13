include("cell_functionalities.jl")

function radiusfunc(initialCell::Cell, desiredCell::Cell, t::Float64, ϕ::Float64)
    return (initialCell.r(ϕ) - desiredCell.r(ϕ))*exp(-t) + desiredCell.r(ϕ) 
end 

function centrefunc(initialCell::Cell, desiredCell::Cell, t::Float64)
    return (initialCell.centre - desiredCell.centre)*exp(-t) + desiredCell.centre
end 



@userplot CellPlotShape
@recipe function fShape(cp::CellPlotShape)
    
    t = cp.args
    label --> "t = " * string(t[1])
    
    z = zeros(2, N)
    for i = 1:N
        ϕ = (i-1)*2*π / N 
        c = centrefunc(initial, desired, convert(Float64, t[1]))
        r = radiusfunc(initial, desired, convert(Float64, t[1]), ϕ)
        z[1,i] = c[1] + r*cos(ϕ)
        z[2,i] = c[2] + r*sin(ϕ)
    end
    
    z[1,:], z[2,:]   
 
end
times = 0:0.1:10
animShape = @animate for t ∈ times

    cellplotshape(t, seriestype=:shape, aspect_ratio=:equal)
    xlims!(-2,8)
    ylims!(-1.5,1.5) 
end
gif(animShape, "Form1_u0_anim.gif", fps = 5)


C1 = circleCell([0.0, 0.0], 1.0)
C2 = peanutCell([5.0, 0.0], 0.3)
C3 = Cell([0.0, 0.0], x -> 1 - 0.2*cos(2*x))

initial = C1 
desired = C2 
N = 500

t = [0.0, 0.5, 1.0, 10.0]
z = zeros(2, N)
time = t[4]
for i = 1:N
    ϕ = (i-1)*2*π / N 
    c = centrefunc(initial, desired, time)
    r = radiusfunc(initial, desired, time, ϕ)
    z[1,i] = c[1] + r*cos(ϕ)
    z[2,i] = c[2] + r*sin(ϕ)
end

ticker = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
p11 = plot(z[1,:], z[2,:], seriestype=:shape, label= false, guidefontsize=14, xaxis=("t = " * string(t[1])), dpi=500, aspect_ratio=:equal, xlims=(-1.5,6.5), ylims=(-1.5,1.5), ticks=ticker)
p12 = plot(z[1,:], z[2,:], seriestype=:shape, label= false, guidefontsize=14, xaxis=("t = " * string(t[2])), dpi=500, aspect_ratio=:equal, xlims=(-1.5,6.5), ylims=(-1.5,1.5), ticks=ticker)
p13 = plot(z[1,:], z[2,:], seriestype=:shape, label= false, guidefontsize=14, xaxis=("t = " * string(t[3])), dpi=500, aspect_ratio=:equal, xlims=(-1.5,6.5), ylims=(-1.5,1.5), ticks=ticker)
p14 = plot(z[1,:], z[2,:], seriestype=:shape, label= false, guidefontsize=14, xaxis=("t = " * string(t[4])), dpi=500, aspect_ratio=:equal, xlims=(-1.5,6.5), ylims=(-1.5,1.5), ticks=ticker)

p21 = plot(z[1,:], z[2,:], seriestype=:shape, label= false, guidefontsize=14, xaxis=("t = " * string(t[1])), dpi=500, aspect_ratio=:equal, xlims=(-1.5,6.5), ylims=(-1.5,1.5), ticks=ticker)
p22 = plot(z[1,:], z[2,:], seriestype=:shape, label= false, guidefontsize=14, xaxis=("t = " * string(t[2])), dpi=500, aspect_ratio=:equal, xlims=(-1.5,6.5), ylims=(-1.5,1.5), ticks=ticker)
p23 = plot(z[1,:], z[2,:], seriestype=:shape, label= false, guidefontsize=14, xaxis=("t = " * string(t[3])), dpi=500, aspect_ratio=:equal, xlims=(-1.5,6.5), ylims=(-1.5,1.5), ticks=ticker)
p24 = plot(z[1,:], z[2,:], seriestype=:shape, label= false, guidefontsize=14, xaxis=("t = " * string(t[4])), dpi=500, aspect_ratio=:equal, xlims=(-1.5,6.5), ylims=(-1.5,1.5), ticks=ticker)



savefig("shapePlotC1_C2_04.png") 
