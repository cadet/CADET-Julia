module CADETJulia

using DifferentialEquations
using HDF5
using CSV
using DataFrames
using ComponentArrays
using Plots
using Sundials

include("solver.jl")
include("solver_dae.jl")


export solve, solver_dae

end