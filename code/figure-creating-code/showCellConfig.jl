include("../parameters.jl") 

NumberOfCellWallPoints = 6                  # number of wall points per cell [OLD NAME: "N"]
N = NumberOfCellWallPoints
NumberOfCells = 2                           # number of cells [OLD NAME: "M"] TODO: change to 400
M = NumberOfCells

c1 = DiscreteCell([0.0018788942990846757, 1.1131757133992795e14, -0.004027864045000419, -0.0065278640450004205, -0.004027864045000423, 1.1131757133992764e14], [1.626303258728257e-19, 6.426922977864271e13, 0.004330127018922193, 6.123233995736766e-19, -0.004330127018922192, -6.42692297786425e13])
c2 = DiscreteCell([0.0058579910639226145, 2.792564706188567e12, -0.00125564214596432, -0.001255642145964321, -8.377694118565703e12, 0.005857991063922613], [0.0024999999999999996, -4.836863954542246e12, 0.0024192503572867674, -0.002419250357286765, -1.4510591863626748e13, -0.0025000000000000022])

u = [c1.x; c2.x; c1.y; c2.y]
X, Y = solutionToXY(u) 
plt = plot(X[1], Y[1],
                seriestype=:shape,
                aspect_ratio=:equal,
                opacity=0.25,
                dpi=500,
                # title=title,
                label=false,
                xlims=domain,
                ylims=domain,
                # xguidefontsize=13,
                # xlabel="t = $(@sprintf("%.1e", time))")
                # xlabel="t = $(round(Int, time/timeStepSize))e-05"
                )

# scatter!(plt, X[1], Y[1])

for i = 2:NumberOfCells
    plot!(plt, X[i], Y[i],
        seriestype=:shape,
        aspect_ratio=:equal,
        opacity=0.25,
        label=false,
    )
end

plot!()
