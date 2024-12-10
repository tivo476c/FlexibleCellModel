### TODO 

change heatmap look:
- they should display: time 

change Parameter.jl into a tree structured structure (.yaml?) 

clean up cell_functionalities.jl and computeOverlap.jl




### Things to talk about with Markus 

interior Angle Problem:
- show instable simulation
- leave it be for now, fix it later?? Do we need it?? Maybe implement a simpler function that just keeps all angles in a range of (20°, 340°) or so 

sim run with all forces but not the interior angle force:
- timestepsize = 2^-10 
- time interval = (0,1)
- 8.135905 seconds (66.28 M allocations: 2.343 GiB, 10.45% gc time, 62.02% compilation time)