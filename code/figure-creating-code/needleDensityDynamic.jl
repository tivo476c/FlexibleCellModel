using Plots
using Printf
using ColorSchemes
using LaTeXStrings
using DifferentialEquations
using LinearAlgebra

# 2D Gaussian with μ = 0 and σ^2 as argument: 
function gaussian2d(x, y, mu, sigma_square)
    return (1 / (2 * π * sigma_square)) * exp(-((x - mu)^2 + (y - mu)^2) / (2 * sigma_square))
end

# PDE dynamic of needle density 
function density_equ_alphamu(du, u, p, t)
    # du and u must be vectors 
    dt, D, ld = p
    dx = 1.0 / Nx
    uMatrix = reshape(u, Nx, Nx)
    ghostM = zeros(Nx + 2, Nx + 2)
    ghostM[2:Nx+1, 2:Nx+1] = uMatrix
    Alpha = zeros(Nx + 2, Nx + 2)
    for i = 2:Nx+1 
        for j = 2:Nx+1 
            Alpha[i, j] = sign(j-i)*(dx*abs(i - j) - ld) 
        end 
    end 
    # copy upper and lower rows and columns
    # upper row 
    ghostM[1, 2:Nx+1] = uMatrix[1, :]
    Alpha[1, 2:Nx+1] = Alpha[2, 2:Nx+1]
    # left column
    ghostM[2:Nx+1, 1] = uMatrix[:, 1]
    Alpha[2:Nx+1, 1] = Alpha[2:Nx+1, 2]
    # lower row 
    ghostM[Nx+2, 2:Nx+1] = uMatrix[Nx, :]
    Alpha[Nx+2, 2:Nx+1] = Alpha[Nx+1, 2:Nx+1]
    # right column
    ghostM[2:Nx+1, Nx+2] = uMatrix[:, Nx]
    Alpha[2:Nx+1, Nx+2] = Alpha[2:Nx+1, Nx+1]
    
    res = zeros(Nx, Nx)
    for i = 1:Nx
        for j = 1:Nx
            
            res[i, j] = 1/(2*dx) * (Alpha[i+1,j+2]*ghostM[i+1,j+2] - Alpha[i+1,j]*ghostM[i+1,j] - Alpha[i+2,j+1]*ghostM[i+2,j+1] + Alpha[i,j+1]*ghostM[i,j+1])
            
        end
    end

    res = vec(res) 
    for i = 1:length(du)
        du[i] = res[i]
    end

    return du
end

## spatial discretisation  
L = 0.5         # usual initial vertex position 
Nx = Ny = 10
x = y = range(0, 2*L, length=Nx)  # grid in both x and y
dx = dy = x[2] - x[1]
X, Y = [x[i] for i in 1:length(x), j in 1:length(y)], [y[j] for i in 1:length(x), j in 1:length(y)]


## time discretisation
dt = 1e-3
NumberOfTimeSteps = 1e3
NumberOfTimeSteps = 1
T = NumberOfTimeSteps * dt
tspan = (0, T)

## inital condition 
sigma_square = 0.09^2

# u0 = [v^1_1, ..., v^400_1, v^1_2, ..., v^400_2] cell vertex positions
u0 = gaussian2d.(X, Y, L, sigma_square)
u0 ./= (sum(u0)*dx^2)

# desired cell length 
ld = 0.2
p = [dt, 1.0, ld]

println("initial mass = $(sum(u0)* dx^2)")
save_times = range(0, T, length=6)
DensProblem = ODEProblem(density_equ_alphamu, vec(u0), tspan, p)
@time sol = solve(DensProblem, Tsit5(); saveat=save_times)



println("last mass = $(sum(sol[end])* dx^2)")

tit = string("N_2(0,0.09^2) performing needle density equation \n with no BC \n\n")
savePath = joinpath(homedir(), "simulations", "density", "needles", "density-evo")
climits = (0, maximum([maximum(sol.u[t]) for t = 1:length(sol.u)]))
for t = 1:length(sol.u)
    t = round(Int64, t)
    u_T = reshape(sol.u[t], Nx, Nx)
    v1string = L"v_1"
    timestring = "t = $(sol.t[t])"
    xlab = "$v1string\n$timestring"
    heatmap(x, y, u_T,
        xlimits=(0.0, 2*L),
        ylimits=(0.0, 2*L),
        xlab=xlab,
        ylab=L"v_2",
        # c=palette(reverse(cgrad(:hot)), 60),
        cgrad=cgrad(reverse(palette(:hot, 8))),
        clim=climits,
        aspect_ratio=:equal,
        dpi=500
    )
    name = "equal-scale-density-t$t.png"
    savefig(joinpath(savePath, name))
    
    heatmap(x, y, u_T,
        xlimits=(0.0, 2*L),
        ylimits=(0.0, 2*L),
        xlab=xlab,
        ylab=L"v_2",
        # c=palette(reverse(cgrad(:hot)), 60),
        cgrad=cgrad(reverse(palette(:hot, 8))),
        # clim=climits,
        aspect_ratio=:equal,
        dpi=500
    )

    name = "dynamic-scale-density-t$t.png"
    savefig(joinpath(savePath, name))
end
