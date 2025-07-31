module CADETJulia


# A file that includes the relevant scripts, refers to the correct functions and loads basic packages for ODE solving

# Local modules
include("DG/DGElements.jl")
using .DGElements
include("DG/ConvDispOperatorDG.jl")
using .ConvDispOperatorDG

#Packages in General registry
# using OrdinaryDiffEq
using Reexport
@reexport using DifferentialEquations
using SparseArrays # To use sparse matrices 
@reexport using DataStructures # To get OrderedDict
using HDF5 # To read HDF5 files
using SpecialFunctions
using LinearAlgebra
using Sundials
using Plots
#using BenchmarkTools
using IncompleteLU

# Include the Julia functions
include("isotherms/binding_base.jl")
include("DG/model_base.jl")
include("utils/flow_rates.jl")
include("utils/cstr.jl")
include("utils/connections.jl")
include("utils/file_reader.jl")
include("utils/repeat_functions.jl")
include("utils/initial_condition_specification.jl")
include("utils/solver_utils.jl")
include("utils/solver.jl")
include("utils/solverDAE.jl")
include("isotherms/Jacobians.jl")
include("DG/conv_disp_operator_dg_jac.jl")

export OrderedDict, solve_model, create_units

end