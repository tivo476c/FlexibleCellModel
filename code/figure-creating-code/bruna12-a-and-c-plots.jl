using Plots
using Printf
using ColorSchemes

using DifferentialEquations
using LinearAlgebra


# 2D Gaussian with μ = 0 and σ^2 as argument: 
function gaussian2d(x, y, sigma_square)
    return (1 / (2 * pi * sigma_square)) * exp(-(x^2 + y^2) / (2 * sigma_square))
end

# PDE dynamic of Figure 2a 
function heat_equ(du, u, p, t)
    # du and u must be vectors 
    uMatrix = reshape(u, Nx, Nx)
    _, D = p

    # interior 
    res = (HeatMatrix * uMatrix + uMatrix * HeatMatrix)

    # reflective boundary condition -> u_0 = u_1 
    for i = 1:Nx
        res[1, i] += uMatrix[1, i]
        res[Nx, i] += uMatrix[Nx, i]
        res[i, 1] += uMatrix[i, 1]
        res[i, Nx] += uMatrix[i, Nx]
    end

    res = D * res / (dx^2)

    res = vec(res)
    for i = 1:length(du)
        du[i] = res[i]
    end

    return du
end

function area_equ(du, u, p, t)

end 

### Simulation parameters  

## spatial discretisation  
Nx = Ny = 300
x = y = range(-0.5, 0.5, length=Nx)  # grid in both x and y
dx = dy = x[2] - x[1]
# x = y = range(-5.0, 5.0, length=2000)  # grid in both x and y
X, Y = [x[i] for i in 1:length(x), j in 1:length(y)], [y[j] for i in 1:length(x), j in 1:length(y)]

## time discretisation
delta_t = 10^-5
# D = 100 * dx^2 / delta_t
D = 100
p = [delta_t, D]
T = 0.05
tspan = (0, T)

## inital condition 
sigma_square = 0.09^2
u0 = gaussian2d.(X, Y, sigma_square)

sum(u0 .* dx^2)

## Problem and solution 
HeatMatrix = zeros(Nx, Nx)
for i = 1:Nx
    for j = 1:Nx
        if i == j
            HeatMatrix[i, j] = -2
        elseif abs(i - j) == 1
            HeatMatrix[i, j] = 1
        end
    end
end

HeatProblem = ODEProblem(heat_equ, vec(u0), tspan, p)
@time sol = solve(HeatProblem, Tsit5(); saveat=[0, 0.01, 0.02, 0.03, 0.04, T])
# @time sol = solve(HeatProblem, Tsit5())


u_t05 = reshape(sol.u[6], Nx, Nx)


tit = string("N(0,0.09)^2 performing Heat equ \n with reflective BC \n\n")

heatmap(x, y, u_t05,
    xlimits=(-0.5, 0.5),
    ylimits=(-0.5, 0.5),
    # xlimits=(-5.0, 5.0),
    # ylimits=(-5.0, 5.0),
    # title="N_2(0,0.09) \n\n",
    # title=tit,
    c=palette(reverse(cgrad(:hot)), 60),
    # cgrad=cgrad(reverse(palette(:hot, 8))),
    clim=(minimum(u_t05), maximum(u_t05)),
    ratio=:equal,
    dpi=300
)

savefig("figures/sanity-check/dt10^-5-small.png")
