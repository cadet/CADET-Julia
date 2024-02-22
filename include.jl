
# A file that includes the relevant scripts, refers to the correct functions and loads basic packages for ODE solving


# Load paths 
rel_path = joinpath(@__DIR__) #Should refer to main overview folder where DG folder is located
push!(LOAD_PATH, joinpath(rel_path,"SRC\\DG")) #Load path to DG folder 

# Include the Julia functions
include(joinpath(rel_path,"SRC\\isotherms\\bindingBase.jl"))
include(joinpath(rel_path,"SRC\\DG\\modelBase.jl"))
include(joinpath(rel_path,"SRC\\utils\\solverUtils.jl"))
include(joinpath(rel_path,"SRC\\isotherms\\Jacobians.jl"))
include(joinpath(rel_path,"SRC\\DG\\ConvDispOperatorDGJac.jl"))
include(joinpath(rel_path,"SRC\\utils\\evaluate_convergence.jl"))

# Import necessary Julia packages and custom modules
using DGElements
using ConvDispOperatorDG
# using OrdinaryDiffEq
using DifferentialEquations
using SparseArrays
using SpecialFunctions,LinearAlgebra

using IncompleteLU
include(joinpath(rel_path,"SRC\\utils\\evaluate_ODEsolvers.jl"))
