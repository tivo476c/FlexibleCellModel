### TODO 

CONTINUE WITH:
- check whether initialization of circle cells works in sanityCheck 
- find good parameters for a heatmap simulation similar to Bruna2012:
    * read Chapman12's simulation setup


    * which diffusitivity constant D, is brownian() correct in SDEProblem(energies, brownian, ...)? try D=1

    * is it ok to let the domain be (-5,5)^2 -> lets try 

    * whats with that one equation that must still hold: Ne^2 = 0.04 (N number of cells, e cell diameter)
    * NoCells M = 100
    * cell radius = 0.01 (-> diameter = 0.02) | BUT MULTIPLY *10  BECAUSE OF LARGER DOMAIN
    * ->> cell radius = 0.1

    * NoWallPoints N = 20 

    - NoSimulations = 10000 how long does 1 sim take ?! 


    * Time interval = 0.05
    * time stepsize = 10^(-5)
    * sample times = 0, 0.01, 0.03, 0.05

    * force scalings = bachelor scalings

- run a nice heatmap sim for markus for all the different overlap types and compare them 
- save the creating files 

change heatmap look:
- they should display: time 

change Parameter.jl into a tree structured structure (.yaml?) 

clean up cell_functionalities.jl and computeOverlap.jl




### Things to talk about with Markus 

Is my diffusitivity constant D chosen right?

interior Angle Problem:
- show instable simulation
- leave it be for now, fix it later?? Do we need it?? Maybe implement a simpler function that just keeps all angles in a range of (20°, 340°) or so 

sim run with all forces but not the interior angle force:
- timestepsize = 2^-10 
- time interval = (0,1)
- 8.135905 seconds (66.28 M allocations: 2.343 GiB, 10.45% gc time, 62.02% compilation time)