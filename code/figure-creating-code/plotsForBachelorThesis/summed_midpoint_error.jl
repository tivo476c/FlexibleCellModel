using Plots, LinearAlgebra

#create mesh 
h = 2.0/10
X = -1:h:1

function psi(x:: Vector{Float64}) :: Float64
    s = 50
    return 0.5 * (1 + tanh(s*(1 - (norm(x,2)))))
    #return 1
end



function errorö(h)
    return abs(approximation(h) - pi)
end

plot( 0.001:0.001:1, realArea, label="Unit Circle Area", xlabel = "h", 
)
plot!(approximation, label="Approximated Area")




