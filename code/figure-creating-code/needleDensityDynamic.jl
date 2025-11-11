using Plots
using Printf
using ColorSchemes

using DifferentialEquations
using LinearAlgebra

# 2D Gaussian with μ = 0 and σ^2 as argument: 
function gaussian2d(x, y, sigma_square)
    return (1 / (2 * pi * sigma_square)) * exp(-(x^2 + y^2) / (2 * sigma_square))
end

# PDE dynamic of needle density 
function density_equ(du, u, p, t)
    # du and u must be vectors 
    uMatrix = reshape(u, Nx, Nx)
    dt, D, ld = p
    dx = 1.0 / Nx

    res = zeros(Nx,Nx)
    # outer values where uMatrix[:,Nx+1] is not defined 
    # we apply reflective BC 
    for i = 1:Nx-1 
        alpha = -(abs(dx*i- dx*Nx) - ld) 
        # first we do the right boundary: 
        res[i,Nx] = uMatrix[i,Nx] + dt * (2*uMatrix[i,Nx] - alpha / dx * (uMatrix[i,Nx] - uMatrix[i+1,Nx]))
        # second we do the lower boundary: 
        res[Nx,i] = uMatrix[Nx,i] + dt * (2*uMatrix[Nx,i] - alpha / dx * (uMatrix[Nx,i+1] - uMatrix[Nx,i]))
    end 
    res[Nx,Nx] = uMatrix[Nx,Nx] + dt * (2*uMatrix[Nx,Nx])
    # inner values with no problems 
    for i = 1:Nx-1
        for j = 1:Nx-1  

            if i >= j
                alpha = abs(dx*i- dx*j) - ld 
            else 
                alpha = -(abs(dx*i- dx*j) - ld) 
            end 
            res[i,j] += uMatrix[i,j] + dt * (2*uMatrix[i,j] - alpha / dx * (uMatrix[i,j+1] - uMatrix[i+1,j]))

        end 
    end 

    res = vec(res)
    for i = 1:length(du)
        du[i] = res[i]
    end

    return du
end


## spatial discretisation  
Nx = Ny = 300
x = y = range(-0.5, 0.5, length=Nx)  # grid in both x and y
dx = dy = x[2] - x[1]
X, Y = [x[i] for i in 1:length(x), j in 1:length(y)], [y[j] for i in 1:length(x), j in 1:length(y)]


## time discretisation
dt = 1e-2
NumberOfTimeSteps = 2e2
T = NumberOfTimeSteps * dt
tspan = (0, T)

## inital condition 
sigma_square = 0.09^2
u0 = gaussian2d.(X, Y, sigma_square)
# desired cell length 
ld = 0.2
p = [dt, 1.0, ld]

println("initial mass = $(sum(u0)* dx^2)")

DensProblem = ODEProblem(density_equ, vec(u0), tspan, p)
@time sol = solve(DensProblem, Tsit5())



println("last mass = $(sum(sol[end])* dx^2)")

tit = string("N_2(0,0.09^2) performing needle density equation \n with no BC \n\n")
savePath = joinpath(homedir(), "simulations", "density", "needles", "density-evo")
climits = (0, maximum([maximum(sol.u[t]) for t=1:7]))
for t = 1:7
    t = round(Int64, t)
    u_T = reshape(sol.u[t], Nx, Nx)
    heatmap(x, y, u_T,
        xlimits=(-0.5, 0.5),
        ylimits=(-0.5, 0.5),
        # c=palette(reverse(cgrad(:hot)), 60),
        cgrad=cgrad(reverse(palette(:hot, 8))),
        clim=climits,
        aspect_ratio=:equal,
        dpi=500
    )

    name = "density-t$t.png"
    savefig(joinpath(savePath, name))
end 
