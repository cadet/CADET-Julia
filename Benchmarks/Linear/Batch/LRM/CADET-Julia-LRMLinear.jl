#Run a convergence test where the results are compared to the semi-analytical solution.


# Add paths for the include file
include(joinpath(@__DIR__, fill("..", 4)..., "include.jl"))

# Specify number of cells, polynomial degree and number of components 
nCell =  [2,4,8,16]
polyDeg = [4,5,6]
nComp = 1

# Preload allocation vectors and matrices
runtime_e = []
maxE_e = []
runtime_i = []
maxE_i = []
DOF = []

nCellu = []
degreeu = []

# Load semi analytical solution
using CSV, DataFrames
c_analytical = CSV.read((joinpath(@__DIR__,"Analytical_LRM_Linear.csv")),DataFrame)


function model_setup(nCells, polyDeg, exactInt, analJac=false, solver=QNDF(autodiff=false))

	# Specify transport parameters for LRMP
	model = LRM(
		nComp = 1, 
		colLength = 1, 
		d_ax = [1e-4], 
		eps_c = 0.6, 
		u = 2/60, 
		switch_time = [0, 60, 130], 
		cIn_c = [[1.0] [0.0]], # Inlet concentration at each switch time
		# cIn_l = [], # Inlet concentration linear slope at each switch time
		# cIn_q = [[0] [0]], # Inlet concentration quadratic slope at each switch time
		# cIn_cube = [[0] [0]], # Inlet concentration cubic slope at each switch time
		polyDeg = polyDeg, # defaults to 4
		nCells = nCells, # defaults to 8
		exactInt = exactInt # 
		)

	bind = Linear(
				ka = [1.0],
				kd = [1.0,],
				is_kinetic = false, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
				nBound = ones(Bool,model.nComp), # Number of bound components, salt not used anyway
				bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
				# nBound =  [1,0,1,1] # Specify non-bound states by a zero, defaults to assume all bound states
				)

	solverOptions = solverCache(model, # Should refer to the model instance
								c0 = 0,
								# cp0 = 0, # only relevant for LRMP, GRM
								q0 = 0,
								# x0 = 0, # only provide if using a whole x0
								abstol = 1e-12,
								# reltol = 1e-10,
								dt = 1,
								# prototypeJacobian = true #defaults to true
								analyticalJacobian = analJac #defaults to false
								)
	output = solve_model(model, bind, solverOptions, solver)
	return output
end

# Evaluate the convergence using the evaluate_convergence function 
using BenchmarkTools
evaluate_convergence(model_setup, c_analytical, nComp, nCell, polyDeg, 1, "LRM", @__DIR__)

# Evaluate ODe solvers and save results in ODETests folder 
# using Plots
# [evaluate_ODEsolvers(model_setup, c_analytical, nComp, nCell, [polyDeg[1]], 1, "LRM", joinpath(@__DIR__,"ODETests")) for _ in 1:2]
