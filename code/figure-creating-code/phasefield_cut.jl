using Plots
using LaTeXStrings

alpha = 12.0

phasefield1(x) = tanh(8*(x+3.2)) - tanh(5*(x-1.9)) - 1
x = -4:0.01:4
y = phasefield1.(x)

plt = plot(
    x,
    y,
    linewidth=3,
    color=:red,
    label=L"\Phi",
    fill = (-1, 0.5, RGB(0.1,0.6,0.1)),
    xlims = [-4,4],
    ylims = [-1.2,1.2],
    aspect_ratio=:equal,
    xlabel = "x",
    # ylabel = false,
    yticks = [-1,0,1], 
    xticks = (-4:1.0:4, ["","","","","","","",""]),
    title="Phase field cut",
    dpi =500,
    )
hline!([0.0], c=:black, linewidth=1, linestyle=:dot, label="Cell wall level set")


savefig(plt, joinpath(homedir(), "simulations", "phasefield_cut.png"))