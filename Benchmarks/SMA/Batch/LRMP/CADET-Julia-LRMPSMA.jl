#Run a convergence test where the results are compared to the semi-analytical solution.


# Add paths for the include file
include(joinpath(@__DIR__, fill("..", 4)..., "include.jl"))

# Specify number of cells, polynomial degree and number of components 
nCell =  [2,4,8,16]
polyDeg = [4,5,6]
nComp = 4


# Load semi analytical solution
using CSV, DataFrames
c_analytical = CSV.read((joinpath(@__DIR__,"Semi-analytical_LRMP_SMA.csv")),DataFrame)

function model_setup(nCells, polyDeg, exactInt, analJac=false, solver=QNDF(autodiff=false))

	# Specify transport parameters for LRMP
	model = LRMP(
		nComp = 4, 
		colLength = 0.014, 
		d_ax = 1e-5, 
		eps_c = 0.37, 
		eps_p = 0.75,
		u = 5.75e-4, 
		kf = [3.3e-3, 3.3e-3, 3.3e-3, 3.3e-3],
		Rp = 4.5e-5,
		switch_time = [0, 10, 90, 1800], 
		cIn_c = [[50,1,1,1] [50,0,0,0] [100,0,0,0]], # Inlet concentration at each switch time
		cIn_l = [[0, 0, 0, 0] [0, 0, 0, 0] [0.2, 0, 0, 0]], # Inlet concentration linear slope at each switch time
		# cIn_q = [[0] [0]], # Inlet concentration quadratic slope at each switch time
		# cIn_cube = [[0] [0]], # Inlet concentration cubic slope at each switch time
		polyDeg = polyDeg, # defaults to 4
		nCells = nCells, # defaults to 8
		exactInt = exactInt # 
		)

	bind = SMA(
				ka = [0.0, 35.5e-3, 1.59e-3, 7.70e-3],
				kd = [1.0, 1, 1, 1],
				ionicCapacity = 1200.0,
				v = [0.0, 4.7, 5.29, 3.7 ], # [-], charge
				sigma = [0.0, 11.83, 10.6, 10.0000], # [-], shielding
				is_kinetic = true, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
				nBound = ones(Bool,model.nComp), # Number of bound components, salt not used anyway
				bindStride = model.bindStride # Not necessary for Linear model, only for Langmuir and SMA 
				# nBound =  [1,0,1,1] # Specify non-bound states by a zero, defaults to assume all bound states
				)

	solverOptions = solverCache(model, # Should refer to the model instance
								c0 = [50, 0, 0, 0],
								# cp0 = 0, # only relevant for LRMP, GRM
								q0 = [1200.0, 0, 0, 0],
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
using BenchmarkTools, Plots
evaluate_convergence(model_setup, c_analytical, nComp, nCell, polyDeg, 1, "LRMP", @__DIR__)

# Evaluate ODe solvers and save results in ODETests folder 
# evaluate_ODEsolvers(model_setup, c_analytical, nComp, nCell, polyDeg, 1, "LRMP", "ODETests//")
