using Plots
using LaTeXStrings

# Define points
points = [
    (0.5, 5),   # P1
    (1.5, 1), # P2
    (5, 1.5),   # P3
    (2.5, 3), # P4
    (6, 4.5)  # P5
]

# Extract x and y
x = [p[1] for p in points]
y = [p[2] for p in points]

# Polygon order (P1 -> P2 -> P4 -> P3 -> P5 -> P1)
poly_order = [1, 2, 4, 3, 5]


# Color-coded arrows
edges_red = [(1, 2), (2, 3), (4, 5)]    # Red edges
edges_green = [(3, 4), (5, 1)]  # Green edges

# Plot polygon and fill
plot(x, y, seriestype=:shape, fillalpha=0.2, linecolor=:lightblue, legend=false, aspect_ratio=1, xlims=(0, 7), ylims=(0, 6), dpi=300)
scatter!(x, y, color=:blue, markerstrokecolor=:black)

# # Draw red arrows
for (i, j) in edges_red
    dx, dy = x[j] - x[i], y[j] - y[i]
    annotate!((x[i] + dx / 2, y[i] + dy / 2), text("→", :red, 14))
    plot!([x[i], x[j]], [y[i], y[j]], arrow=true, color=:red, linewidth=2)
end

# # Draw green arrows
for (i, j) in edges_green
    dx, dy = x[j] - x[i], y[j] - y[i]
    annotate!((x[i] + dx / 2, y[i] + dy / 2), text("→", :green, 14))
    plot!([x[i], x[j]], [y[i], y[j]], arrow=true, color=:green, linewidth=2)
end

# # Label points
labels = [L"\vec{v}_1", L"\vec{v}_2", L"\vec{v}_3", L"\vec{v}_4", L"\vec{v}_5"]
for (i, label) in enumerate(labels)
    annotate!(x[i] + 0.2, y[i] + 0.4, label)
end

## Add vertical lines from each point to the x axis 
for p in points

    x = [p[1], p[1]]
    y = [p[2], 0]
    plot!(x, y, color="black")

end

# # Add green and red plus/minus signs
symbols = [
    ("+", 0.85, 1),
    ("-", 1.1, 1),
    ("+", 1.85, 0.5),
    ("-", 2.1, 0.5),
    ("++", 3.6, 0.7),
    ("- -", 4.2, 0.7),
    ("+", 4.5, 2.8),
    ("-", 4.8, 2.8),
    ("+", 1.8, 3.5),
    ("++", 3.2, 1.8),
    ("-", 3.7, 1.8),
]

for (s, px, py) in symbols
    color = occursin("-", s) ? :red : :green
    annotate!(px, py, text(s, color, 14))
end

# Show plot
savefig("shoelace_new.png")