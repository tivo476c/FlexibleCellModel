using Plots
using Printf
using ColorSchemes

using DifferentialEquations
using LinearAlgebra

# 2D Gaussian with μ = 0 and σ^2 as argument: 
function gaussian2d(x, y, sigma_square)
    return (1 / (2*pi*sigma_square)) * exp(- (x^2 + y^2) / (2*sigma_square))
end

function laplacian_1d(N, dx)
    A = spzeros(N, N)
    for i in 2:N-1
        A[i, i-1] = -1
        A[i, i]   =  2
        A[i, i+1] = -1
    end
    # Neumann BCs (zero derivative at boundaries)
    A[1,1] = 1
    A[1,2] = -1
    A[N,N-1] = -1
    A[N,N] = 1

    return A / dx^2
end

# PDE dynamic of Figure 2a 
function heat_equ!(du, u, p, t)
    laplace_x = laplacian_1d(Nx, dx)
    Identity = spdiagm(0 => ones(Nx))
    L = LinearAlgebra.kron(Identity, laplace_x) + LinearAlgebra.kron(laplace_x, Identity)
    
    du .= D * (L * u)
end 

function heat2(du, u, p, t)
    laplace_u = zeros(size(u))
    for i in 2:size(u, 1)-1
        for j in 2:size(u, 2)-1
            laplace_u[i, j] = (u[i+1, j] - 2*u[i, j] + u[i-1, j]) / dx^2 + 
                        (u[i, j+1] - 2*u[i, j] + u[i, j-1]) / dy^2
        end
    end
    du .= u .* laplace_u

    #reflective boundary condition: 
    du[1, :] .= 0
    du[end, :] .= 0
    du[:, 1] .= 0
    du[:, end] .= 0
end 


### Simulation parameters  

## spatial discretisation  
Nx = Ny = 100
x = y = range(-0.5, 0.5, length=Nx)  # grid in both x and y
dx = dy = x[2] - x[1] 
# x = y = range(-5.0, 5.0, length=2000)  # grid in both x and y
X, Y = [x[i] for i in 1:length(x), j in 1:length(y)], [y[j] for i in 1:length(x), j in 1:length(y)]

## time discretisation
delta_t = dx^2
# D = 100 * dx^2 / delta_t
D = 1
p = [delta_t, D]
T = 0.05
tspan = (0, T) 

## inital condition 
sigma_square = 0.09
u0 = gaussian2d.(X, Y, sigma_square)

## Problem and solution 
HeatProblem = ODEProblem(heat2, u0, tspan, p)
@time sol = solve(HeatProblem, Tsit5(); saveat=T)

u_t05 = sol.u[1]


heatmap(x, y, u_t05,
    xlimits = (-0.5,0.5), 
    ylimits = (-0.5,0.5), 
    # xlimits = (-5.0,5.0), 
    # ylimits = (-5.0,5.0), 
    # title="N(0,0.09)^2 \n\n",
    title="N(0,0.09)^2 performing Heat equ \n at t=0.05 with reflective BC \n\n",
    c=reverse(cgrad(:hot)),
    # clim = (0.55, 1.55), 
    ratio=:equal,
    dpi=500
)

savefig("figures/sanity-check/normal-distribution-T0_005-.png")