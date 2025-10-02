using Plots
using Printf
using ColorSchemes

include("../cell_functionalities.jl")
# include("../parameters.jl")

function createLocationFile(sol, sim::Int64, locationsPath)

    sim = @sprintf("%07d", sim)
    fileName = string("locations-", simulationName, "-", sim, ".txt")
    filePath = joinpath(locationsPath, fileName)

    open(filePath, "w") do file

        for u in sol.u
            X, Y = solutionToXY(u)
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


function mirrowhorizontal!(matrices)

    for i in eachindex(matrices)
        secondMatrix = reverse(matrices[i], dims=2)
        matrices[i] += secondMatrix
        matrices[i] .*= 0.5
    end

end

function mirrowvertical!(matrices)

    for i in eachindex(matrices)
        secondMatrix = reverse(matrices[i], dims=1)
        matrices[i] += secondMatrix
        matrices[i] .*= 0.5
    end

end

function addTransposed!(matrices)

    for i in eachindex(matrices)
        secondMatrix = transpose(matrices[i])
        matrices[i] += secondMatrix
        matrices[i] .*= 0.5
    end

end



function smoothenMatrix!(matrices)
    """
    this function enhances the number of sims by copying matrices in different ways. 
    """

    mirrowhorizontal!(matrices)
    mirrowvertical!(matrices)
    addTransposed!(matrices)

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

function addSolToMatrices(solution, matrices)

    for counter = 1:NumberOfSampleTimes

        # extract centres from solution 
        X, Y = solutionToXY(solution[counter])
        centresX = zeros(NumberOfCells)
        centresY = zeros(NumberOfCells)
        if NumberOfCellWallPoints == 0
            centresX = X
            centresY = Y
        else
            for i = 1:NumberOfCells
                centresX[i] = sum(X[i]) / Float64(NumberOfCellWallPoints)
                centresY[i] = sum(Y[i]) / Float64(NumberOfCellWallPoints)
            end
        end

        # add centres to matrices
        for i = 1:NumberOfCells
            coords = [centresX[i], centresY[i]]
            row, column = getMatrixIndex(coords)
            matrices[counter][row, column] += 1
        end
    end
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
    matrices = [matrices[i] ./ (NumberOfSimulations * NumberOfCells * (HeatStepSize)^2) for i = 1:NumberOfSampleTimes]
    oldMatrices = copy(matrices)
    # smoothenMatrix!(matrices)
    maxVal = maximum([maximum(matrices[i]) for i = 1:NumberOfSampleTimes])
    minVal = minimum([minimum(matrices[i]) for i = 1:NumberOfSampleTimes])
    println("NumberOfSimulations = $NumberOfSimulations \nNumberOfCells = $NumberOfCells \nHeatStepSize = $HeatStepSize")
    mass1 = sum(matrices[1]) * HeatStepSize^2
    massN = sum(matrices[NumberOfSampleTimes]) * HeatStepSize^2
    println("mass1 = $mass1, massN = $massN")

    for i = 1:NumberOfSampleTimes

        # sampleTime = sampleTimes[i]
        # heatMapName = string("heatmap-", simulationName, "-sampleTime", @sprintf("%.3f", sampleTime), "bruna12scale")
        # # heatMapName = string("heatmap-sampleTime", sampleTime, "bruna12scale.png")
        # title = string("Heatmap of simulation '", simulationName, "'")
        # # caption = string("Number of cells = $NumberOfCells")
        # caption = string("Number of simulations: ", 10000, ", sample time t = ", @sprintf("%.4f", sampleTime))

        # heatmap(HeatGrid, HeatGrid, matrices[i],
        #     xlimits=domain,
        #     ylimits=domain,
        #     xlabel=caption,
        #     c=reverse(cgrad(:hot)),
        #     # clim=(minVal, maxVal),
        #     clim=(0.55, 1.55),  # activate for bruna scaling 
        #     ratio=:equal,
        #     dpi=500
        # )
        # vline!(HeatGrid, c=:grey, linewidth=0.1, label=false)
        # hline!(HeatGrid, c=:grey, linewidth=0.1, label=false)
        # savefig(joinpath(heatMapsPath, "$(heatMapName).png"))


        # # plot without color legend
        # heatmap(HeatGrid, HeatGrid, matrices[i],
        #     xlimits=domain,
        #     ylimits=domain,
        #     xlabel=caption,
        #     c=reverse(cgrad(:hot)),
        #     # clim=(minVal, maxVal),
        #     clim=(0.55, 1.55),  # activate for bruna scaling 
        #     ratio=:equal,
        #     dpi=500,
        #     colorbar=false,
        #     margin=0Plots.mm
        # )
        # vline!(HeatGrid, c=:grey, linewidth=0.1, label=false)
        # hline!(HeatGrid, c=:grey, linewidth=0.1, label=false)
        # savefig(joinpath(heatMapsPath, "$(heatMapName)-nolegend.png"))

        heatmap(
            HeatGrid, HeatGrid, matrices[i],
            xlimits=domain,
            ylimits=domain,
            # xlabel=caption,
            c=reverse(cgrad(:hot)),
            clim=(0.55, 1.55),
            ratio=:equal,
            size=(1000, 800),          # width x height in pixels (affects aspect)
            colorbar=true,
            colorbar_position=:bottom     # horizontal
        )

        vline!(HeatGrid, c=:grey, linewidth=0.1, label=false)
        hline!(HeatGrid, c=:grey, linewidth=0.1, label=false)
        
        plot!(Plots.current(), colorbar=:horizontal)
        # Save as scalable PDF
        savefig(joinpath(heatMapsPath, "colorbar.pdf"))

    end
end
