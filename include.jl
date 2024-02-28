
# A file that includes the relevant scripts, refers to the correct functions and loads basic packages for ODE solving


# Load paths 
rel_path = joinpath(@__DIR__) #Should refer to main overview folder where DG folder is located
push!(LOAD_PATH, joinpath(rel_path,"SRC\\DG")) #Load path to DG folder 

# Import necessary Julia packages and custom modules
using DGElements
using ConvDispOperatorDG
# using OrdinaryDiffEq
using DifferentialEquations
using SparseArrays # To use sparse matrices 
using DataStructures # To get OrderedDict
using SpecialFunctions,LinearAlgebra

using IncompleteLU

# Include the Julia functions
include(joinpath(rel_path,"SRC\\isotherms\\bindingBase.jl"))
include(joinpath(rel_path,"SRC\\DG\\modelBase.jl"))
include(joinpath(rel_path,"SRC\\utils\\connections.jl"))
include(joinpath(rel_path,"SRC\\utils\\dictionaryReader.jl"))
include(joinpath(rel_path,"SRC\\utils\\repeatFunctions.jl"))
include(joinpath(rel_path,"SRC\\utils\\initialConditionSpecification.jl"))
include(joinpath(rel_path,"SRC\\utils\\solverUtils.jl"))
include(joinpath(rel_path,"SRC\\utils\\solver.jl"))
include(joinpath(rel_path,"SRC\\isotherms\\Jacobians.jl"))
include(joinpath(rel_path,"SRC\\DG\\ConvDispOperatorDGJac.jl"))
include(joinpath(rel_path,"SRC\\utils\\evaluate_convergence.jl"))
include(joinpath(rel_path,"SRC\\utils\\evaluate_ODEsolvers.jl"))



