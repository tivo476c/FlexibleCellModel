using Pkg
"""
This file is the entry point of my FlexibleCellModel code. 
It gets executed, whenever parameters.jl gets executed.
It uses its parameter configuration to start the according simulation. 
""" 
Pkg.develop(path="C:/Users/voglt/Desktop/FlexibleCellModel")

println("Started entryPoint.jl")

include("figure-creating-code/sanityCheck.jl")
