using Plots
using LaTeXStrings

# Define the parameters
L = 5.0
eps = 2.0

# Create the plot
domain = [-L, L]
ticks = [-L, -L + eps/2, 0, L - eps/2, L]
tick_labels = [L"-L", L"-L + \frac{\epsilon}{2}", L"0", L"L - \frac{\epsilon}{2}", L"L"]

plt = plot(
    xlims = domain,
    ylims = domain,
    aspect_ratio=:equal,
    xlabel = "Position of cell 1",
    ylabel = "Position of cell 2",
    xticks = (ticks, tick_labels),
    yticks = (ticks, tick_labels),
    xguidefont = font(10),
    yguidefont = font(10),
    legend = true,
    # color = blue,  # RGB(173/255, 216/255, 230/255),
    dpi =500,
    )

myRed = RGB(216/255,48/255,48/255) 
hspan!(plt, domain, color=RGB(135/255, 206/255, 235/255), alpha=0.2, labels = "feasible area")
hspan!(plt, [-L, -L+eps/2], color=myRed, alpha=1.0, labels = "excluded area")
hspan!(plt, [L- eps/2, L ], color=myRed, alpha=1.0, labels =false)
vspan!(plt, [-L, -L+eps/2], color=myRed, alpha=1.0, labels =false)
vspan!(plt, [L- eps/2, L ], color=myRed, alpha=1.0, labels =false)
plot!(plt, domain, domain,  color=myRed, alpha=1.0, linewidth=34*eps, label=false)
vline!(ticks, c=:black, linewidth=0.2, label=false)
hline!(ticks, c=:black, linewidth=0.2, label=false)


savefig(plt, joinpath(homedir(), "simulations", "hardsphere_exclusion.png"))