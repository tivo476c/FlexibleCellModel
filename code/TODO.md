

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


QUESTIONS: 
* use \mathrm{d} also in \int ... \mathrm{d} x or just for SDEs?