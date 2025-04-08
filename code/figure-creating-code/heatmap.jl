using Plots
using Printf
using ColorSchemes

include("../cell_functionalities.jl")

function createLocationFile(sol, sim::Int64, locationsPath) 
    
    sim = @sprintf("%04d", sim)
    fileName = string("locations-", simulationName, "-", sim, ".txt")
    filePath = joinpath(locationsPath, fileName)

    open(filePath, "w") do file 

        for u in sol.u
            X,Y = solutionToCells(u)
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

    matrices = [zeros(Int64, 40, 40) for _ in 1:NumberOfSampleTimes]
    for file in readdir(locationsPath)
        runPath = joinpath(locationsPath,file)
        i=1
        open(runPath, "r") do file_stream
            for line in eachline(file_stream)
                coords = parse.(Float64, split(line))
                if length(coords) == 2
                    # centre coordinates obtained 
                    row, column = getMatrixIndex(coords)
                    matrices[i][row, column] += 1
                else
                    i += 1
                    if(i>NumberOfSampleTimes) 
                        break 
                    end 
                end
            end
        end 
    end 

    return matrices 

end 

function getMatrixIndex(coords::Vector{Float64})
    x,y = coords
    k = trunc(4*x) 
    column = Int(k + 21) 
    l = trunc(4*y) 
    row = Int(20 - l)          # = 41 - (4*l + 21) 
    if(row==0)
        row = 1
    end 
    if(column==41)
        column = 40
    end  

    return row, column
end

function createHeatmaps(matrices)
    
    for i = 1:NumberOfSampleTimes

        sampleTime = sampleTimes[i]
        matrix = matrices[i] ./  NumberOfSimulations
        maxValue = maximum(matrix)        
        
        heatMapName = string("heatmap-", simulationName, "-sampleTime", sampleTime, ".png") 
        title = string("Heatmap of simulation '", simulationName, "'")
        caption = string("Number of simulations: ", NumberOfSimulations, ", sample time t = ", @sprintf("%.4f", sampleTime))

        grid = -5.0:0.25:5.0
        heatmap(grid, grid, matrix,
            xlimits = (-5.0,5.0), 
            ylimits = (-5.0,5.0), 
            xlabel=caption, 
            c=reverse(cgrad(:hot)),
            colorrange = (0, maxValue), 
            ratio=:equal,
            dpi=500
            )
        vline!(-5.0:0.25:5.0, c=:grey, linewidth=0.1, label=false)
        hline!(-5.0:0.25:5.0, c=:grey, linewidth=0.1, label=false)

        savefig(joinpath(heatMapsPath, heatMapName))
    end 
end
