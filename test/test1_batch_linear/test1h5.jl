
using Test, CADETJulia, CSV, DataFrames



filename = joinpath(@__DIR__,"test1.h5")
file = h5open(filename, "r")
model = file


inlets, outlets, columns, switches, solverOptions = create_units(model)

solve_model(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, # Defaults to (0,) as output is also written to units 
			alg = QNDF(autodiff=false), # Defaults to alg = QNDF(autodiff=false)
			)

#plot(columns[1].solution_outlet[:,1])

# Compare to analytical solution 
using CSV,DataFrames 
c_analytical = CSV.read((joinpath(@__DIR__,"Analytical_LRM_Linear.csv")),DataFrame)
err = maximum(abs.(columns[1].solution_outlet[:,1]-c_analytical[:,"C0"]))

@test err < 1e-5
