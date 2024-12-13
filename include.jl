
# A file that includes the relevant scripts, refers to the correct functions and loads basic packages for ODE solving


# Load paths 
rel_path = joinpath(@__DIR__) #Should refer to main overview folder where DG folder is located
push!(LOAD_PATH, joinpath(rel_path, "src", "DG")) #Load path to DG folder 

# Import necessary Julia packages and custom modules
using DGElements
using ConvDispOperatorDG
# using OrdinaryDiffEq
using DifferentialEquations
using SparseArrays # To use sparse matrices 
using DataStructures # To get OrderedDict
using HDF5 # To read HDF5 files
using SpecialFunctions,LinearAlgebra
using Sundials
using Plots
using BenchmarkTools

using IncompleteLU

# Include the Julia functions
include(joinpath(rel_path,"src", "isotherms", "binding_base.jl"))
include(joinpath(rel_path,"src", "DG", "model_base.jl"))
include(joinpath(rel_path,"src", "utils", "flow_rates.jl"))
include(joinpath(rel_path,"src", "utils", "cstr.jl"))
include(joinpath(rel_path,"src", "utils", "connections.jl"))
include(joinpath(rel_path,"src", "utils", "file_reader.jl"))
include(joinpath(rel_path,"src", "utils", "repeat_functions.jl"))
include(joinpath(rel_path,"src", "utils", "initial_condition_specification.jl"))
include(joinpath(rel_path,"src", "utils", "solver_utils.jl"))
include(joinpath(rel_path,"src", "utils", "solver.jl"))
include(joinpath(rel_path,"src", "utils", "solverDAE.jl"))
include(joinpath(rel_path,"src", "isotherms", "Jacobians.jl"))
include(joinpath(rel_path,"src", "DG", "conv_disp_operator_dg_jac.jl"))



