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

    * NoCells M = 200  
    * diameter = 1/(sqrt(2)*50)
    * cell radius = 1/(sqrt(2)*100) (-> diameter = 0.02) | BUT MULTIPLY *10  BECAUSE OF LARGER DOMAIN
    * ->> cell radius = 1/(sqrt(2)*1000)

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
- D ~ dx^2 / dt, we need to rescale dx with factor 10 --> D = 100? 

interior Angle Problem:
- show instable simulation
- leave it be for now, fix it later?? Do we need it?? Maybe implement a simpler function that just keeps all angles in a range of (20°, 340°) or so 
- let it be for now

sim run with all forces but not the interior angle force:
- timestepsize = 2^-10 
- time interval = (0,1)
- 8.135905 seconds (66.28 M allocations: 2.343 GiB, 10.45% gc time, 62.02% compilation time)

sim run with all forces but not the interior angle force:
- timestepsize = 10^-4
- time interval = (0,0.05)
- N = 20                              # number of wall points per cell 
- M = 200                             # number of cells 
- D = 100                             # diffusitivity constant 
- radius = 1/(sqrt(2)*1000)           # cell radius 
- 34.911674 seconds (712.24 M allocations: 37.026 GiB, 15.61% gc time)


sim run with all forces but not the interior angle force:
- timestepsize = 10^-5 [everything is the same as one above]
- time interval = (0,0.05)
- N = 20                              # number of wall points per cell 
- M = 200                             # number of cells 
- D = 100                             # diffusitivity constant 
- radius = 1/(sqrt(2)*1000)           # cell radius 
- 335.656799 seconds (7.12 G allocations: 370.253 GiB, 15.45% gc time)


TODOs:
- rescale heatmaps such that:  ... / (noSimulations*NoCells*deltaX*deltaY) (deltaX = deltaY = 1/4) 
- DO THE POINTPARTICLE HEATMAP 1st (just brownian motion nothing else, with D=100)
- for heatmap use color: hot, or cmap=Hot, oder so
- parallelize simulazion look at using distributed julia documentation 
- new overlap scaling: IF overlap x_j = x_j + (2*r - |x_i - x_j|)(x_j-x_i)/|x_i - x_j|
- DO SIMULATION JUST WITH THIS OVERLAP AND NO DEFORMATION 
- THEN WITH DEFORMATION 

- SIMULATION PARAMETERS: 

    * delta t = 10^-3
    * NoCells M = 200  
    * NoWallpoints = 20 (control whether its enough for deformation pde)
    * diameter = 1/(sqrt(2)*50)
    * cell radius = 1/(sqrt(2)*100) (-> diameter = 0.02) | BUT MULTIPLY *10  BECAUSE OF LARGER DOMAIN
    * ->> cell radius = 1/(sqrt(2)*1000)
    * noSimulations: look where the point particles do a nice gaussian 
    * cells get initialized via: N(0,0.9)^2 (in paper: N(0, 0.09)) 

- START WRITING LATEX PDF!!! , also add different overlap forces
