

# Add the include file custom to load packages and scripts. 
# the file is located on the main from which the file takes care of the rest. 
include(joinpath(@__DIR__,"..\\..\\include.jl"))



filename = joinpath(@__DIR__,"model1DG.h5")
file = h5open(filename, "r")
model = file


# read(file["input/model"]["unit_000"]["sec_001"]["CONST_COEFF"])
inlets, outlets, columns, switches, solverOptions = create_units(model)

solve_model(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, # Defaults to (0,) as output is also written to units 
			alg = QNDF(autodiff=false), # Defaults to alg = QNDF(autodiff=false)
			)

using Plots
plot(columns[1].solution_outlet[:,1])

# Compare to analytical solution 
using CSV,DataFrames 
c_analytical = CSV.read((joinpath(@__DIR__,"Analytical_LRM_Linear.csv")),DataFrame)
err = maximum(abs.(columns[1].solution_outlet[:,1]-c_analytical[:,"C0"]))
if err<1e-5
    println("Test succesful - error lower than 1e-5")
else
    println("Test unsuccesful - error larger than 1e-5")
end

