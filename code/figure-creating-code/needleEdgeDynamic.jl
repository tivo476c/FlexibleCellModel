using Distributions
using Plots

N = 400         # number of cells 
ld = 0.2        # desired needle length 
L = 0.5         # half domain length 


# initialize all vertices 
C = rand(Truncated(Normal(0, 0.09), -L, L), 2 * N)  # full cell vector with [N*x, N*y] (first all x coords, then all y coords)

# now do the dynamic: 

dt = 1e-2
NumberOfTimeSteps = 2e2
T = NumberOfTimeSteps * dt
# numSteps = Int(T / dt)
NumberOfSampleTimes = 6
sol = [zeros(2 * N) for i in 1:NumberOfSampleTimes]

timestep = 0

sampleTimeFactor = round(Int64, (T / dt) * (1 / (NumberOfSampleTimes - 1)))
# this is the integer index, st. k*sampleTimeFactor for k=0:NumberOfSampleTimes-1 gives the index to access the according sampletimes for the solution vector. 
while timestep <= T

    for cell = 1:N

        v1 = C[cell]
        v2 = C[cell+N]
        l = abs(v1 - v2)

        C[cell] += dt * (-sign(v1 - v2) * (l - ld)^1)
        C[cell+N] += dt * (sign(v1 - v2) * (l - ld)^1)

    end

    # save at sample times 
    k = round(Int64, timestep / dt)

    if k == 0 || mod(k, sampleTimeFactor) == 0
        sol[round(Int64, k / sampleTimeFactor + 1)] = copy(C)
    end
    global timestep += dt
    global timestep = round(timestep, sigdigits=3)


end


for i = 1:length(sol)
    midCellLength = 1.0/N * sum(abs.(sol[i][1:N] - sol[i][N+1:2*N]))
    println("midCellLength = $midCellLength at t=$i")
end 

savePath = joinpath(homedir(), "simulations", "density", "needles", "histograms")
# if !isdir(savePath)
#     mkdir(savePath)
# end

for sampleTime = 1:NumberOfSampleTimes

    scatterplt = scatter(
        sol[sampleTime][1:N],
        sol[sampleTime][N+1:2*N],
        xlab="v1\nt=$(round((sampleTime-1)*sampleTimeFactor*dt, sigdigits=3))",
        ylab="v2",
        lab="cells",
        aspect_ratio=:equal,
        legend=:topleft,
        xlims=(-L, L),
        ylims=(-L, L),
        dpi=500,
    )
    name = "histogram_t$(sampleTime).png"
    savefig(joinpath(savePath, name))


end
