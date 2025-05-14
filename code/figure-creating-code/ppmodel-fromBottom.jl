include("heatmap.jl")
include("../energies.jl")
using Distributions

NumParticle = 400
NumSims = 1000
NumberOfHeatGridPoints = 50
domainL = 0.5
HeatStepSize = 2 * domainL / NumberOfHeatGridPoints
HeatGrid = [-domainL + k * HeatStepSize for k in 0:NumberOfHeatGridPoints]


T = 0.05
dt = 10^-3
D = 1


initialDistribution = MvNormal([0, 0], 0.09^2 * I(2))


matrix = zeros(Int64, NumberOfHeatGridPoints, NumberOfHeatGridPoints)
for i_sim = 1:NumSims

    # println("starting with sim ", i_sim)

    # initialize 
    allcoords = rand(initialDistribution, NumParticle) # first NumParticle coordinates for x coordinate, others for y coordinate 
    currentState = [allcoords[1, :]; allcoords[2, :]]
    # compute solution at T=0.05 
    time = 0
    while time < T
        # energies(du, u, p, t)
        du = zeros(800)
        energies!(du, currentState, (dt, D), time)
        currentState += dt * du + sqrt(2 * D * dt) * randn(2 * NumParticle)
        # for i = 1:2*NumParticle
        #     if currentState[i] < -0.5
        #         currentState[i] = -1.0 - currentState[i]  # reflect back into domain
        #     elseif currentState[i] > 0.5
        #         currentState[i] = 1.0 - currentState[i]   # reflect back into domain
        #     end
        # end
        time += dt
    end

    # add to matrix
    for i = 1:NumParticle
        coords = [currentState[i], currentState[i+NumParticle]]
        row, column = getMatrixIndex(coords)
        matrix[row, column] += 1
    end
end

# plot 
matrix = matrix ./ (NumSims * NumParticle * (HeatStepSize)^2)

mass1 = sum(matrix) * HeatStepSize^2
println("mass1 = ", mass1)


heatmap(HeatGrid, HeatGrid, matrix,
    xlimits=(-domainL, domainL),
    ylimits=(-domainL, domainL),
    c=reverse(cgrad(:hot)),
    # clim=(minVal, maxVal),
    clim=(0.55, 1.55),  # activate for bruna scaling 
    ratio=:equal,
    dpi=500
)
vline!(HeatGrid, c=:grey, linewidth=0.1, label=false)
hline!(HeatGrid, c=:grey, linewidth=0.1, label=false)

savefig("bottomTest/testPPModel-includingEnergies-newBoundaryCondition2.png")

