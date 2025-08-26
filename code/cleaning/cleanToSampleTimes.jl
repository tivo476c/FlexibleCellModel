# save as filter_timesteps.jl
# run with: julia filter_timesteps.jl

# Parameters
dt = 1e-5
sample_times = 0.00:0.01:0.05  # times we want to keep 

# Compute which step indices to keep
sample_steps = Int.(round.(sample_times ./ dt)) .+ 1

# Process every .txt file in the same directory
for file in filter(f -> endswith(f, ".txt"), readdir("."))
    println("Processing $file ...")

    # Read all lines from file
    lines = readlines(file)

    # Split into blocks separated by blank lines
    blocks = Vector{Vector{String}}()
    current = String[]
    for line in lines
        if isempty(strip(line))
            push!(blocks, current)
            current = String[]
        else
            push!(current, line)
        end
    end
    if !isempty(current)
        push!(blocks, current)
    end

    println("  Found $(length(blocks)) time steps")

    # Keep only the wanted sample steps
    keep_blocks = [blocks[i] for i in sample_steps if i ≤ length(blocks)]

    if length(keep_blocks) == 6
        # Write reduced output to new file
        outname = replace(file, ".txt" => "_filtered.txt")
        # outname = file
        open(outname, "w") do io
            for (k, block) in enumerate(keep_blocks)
                for line in block
                    println(io, line)
                end
                if k < length(keep_blocks)
                    println(io)  # blank line between timesteps
                end
            end
        end

        println("  Wrote $outname with $(length(keep_blocks)) time steps")
        rm(file)
        println("  Deleted original $file")
    end 
end
