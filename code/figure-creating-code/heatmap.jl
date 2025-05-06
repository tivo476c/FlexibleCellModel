using Plots
using Printf
using ColorSchemes

include("../cell_functionalities.jl")
include("../parameters.jl")

function createLocationFile(sol, sim::Int64, locationsPath)

    sim = @sprintf("%04d", sim)
    fileName = string("locations-", simulationName, "-", sim, ".txt")
    filePath = joinpath(locationsPath, fileName)

    open(filePath, "w") do file

        for u in sol.u
            X, Y = solutionToCells(u)
            if NumberOfCellWallPoints != 0
                for i = 1:M
                    xCentre = round(sum(X[i]) / Float64(N), digits=3)
                    yCentre = round(sum(Y[i]) / Float64(N), digits=3)
                    println(file, "$xCentre $yCentre")
                end
            else
                for i = 1:M
                    xCentre = round(X[i], digits=3)
                    yCentre = round(Y[i], digits=3)
                    println(file, "$xCentre $yCentre")
                end
            end
            println(file)
        end
    end


end

function makeMatrices()
    # ONE MATRIX STORES THE ENTRIES FOR ONE SAMPLE TIME AT ALL SIMULATION ITERATIONS 

    matrices = [zeros(Int64, NumberOfHeatGridPoints, NumberOfHeatGridPoints) for _ in 1:NumberOfSampleTimes]
    for file in readdir(locationsPath)
        runPath = joinpath(locationsPath, file)
        i = 1
        open(runPath, "r") do file_stream
            for line in eachline(file_stream)
                coords = parse.(Float64, split(line))
                if length(coords) == 2
                    # centre coordinates obtained 
                    row, column = getMatrixIndex(coords)
                    matrices[i][row, column] += 1
                else
                    i += 1
                    if (i > NumberOfSampleTimes)
                        break
                    end
                end
            end
        end
    end

    return matrices

end

function getMatrixIndex(coords::Vector{Float64})

    x, y = coords

    i = NumberOfHeatGridPoints - (Int(floor((y + domainL) / HeatStepSize)))     # row (y-axis)
    j = Int(floor((x + domainL) / HeatStepSize)) + 1            # column (x-axis)
    
    if i > NumberOfHeatGridPoints
        i = NumberOfHeatGridPoints
    end
    if i < 1
        i = 1
    end
    if j > NumberOfHeatGridPoints
        j = NumberOfHeatGridPoints
    end
    if j < 1
        j = 1
    end
    return i, j

end

function createHeatmaps(matrices)

    matrices .= [Float64.(M) for M in matrices]
    maxVal = maximum([maximum(matrices[i]) for i = 1:NumberOfSampleTimes])
    minVal = minimum([minimum(matrices[i]) for i = 1:NumberOfSampleTimes])
    println("old minVal = ", minVal, "; old maxVal = ", maxVal) 
    matrices = [matrices[i]./(NumberOfSimulations * NumberOfCells * (HeatStepSize)^2) for i=1:NumberOfSampleTimes]

    maxVal = maximum([maximum(matrices[i]) for i = 1:NumberOfSampleTimes])
    minVal = minimum([minimum(matrices[i]) for i = 1:NumberOfSampleTimes])
    println("new minVal = ", minVal, "; new maxVal = ", maxVal) 
    mass1 = sum(matrices[1]) * HeatStepSize^2
    massN = sum(matrices[NumberOfSampleTimes]) * HeatStepSize^2
    println("mass1 = ", mass1, "; massN = ", massN)

    for i = 1:NumberOfSampleTimes

        sampleTime = sampleTimes[i]
        # heatMapName = string("heatmap-", simulationName, "-sampleTime", @sprintf("%.4f", sampleTime), ".png")
        heatMapName = string("heatmap-", simulationName, "-sampleTime", sampleTime, "bruna12scale.png")
        title = string("Heatmap of simulation '", simulationName, "'")
        caption = string("Number of simulations: ", NumberOfSimulations, ", sample time t = ", @sprintf("%.4f", sampleTime))

        heatmap(HeatGrid, HeatGrid, matrices[i],
            xlimits=domain,
            ylimits=domain,
            xlabel=caption,
            c=reverse(cgrad(:hot)),
            # clim=(minVal, maxVal),
            clim=(0.55, 1.55),  # activate for bruna scaling 
            ratio=:equal,
            dpi=500
        )
        vline!(HeatGrid, c=:grey, linewidth=0.1, label=false)
        hline!(HeatGrid, c=:grey, linewidth=0.1, label=false)

        savefig(joinpath(heatMapsPath, heatMapName))

    end
end
