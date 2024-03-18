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

function model_setup(nCells, polyDeg, exactInt, analJac=false)


	nComp = 4
	model = OrderedDict(
		"root" => OrderedDict(
			"input" => OrderedDict(
				"model" => OrderedDict()
			)
		)
	)


	# Set elements sequentially for unit_000
	model["root"]["input"]["model"]["unit_000"] = OrderedDict()
	model["root"]["input"]["model"]["unit_000"]["unit_type"] = "INLET"
	model["root"]["input"]["model"]["unit_000"]["ncomp"] = nComp
	model["root"]["input"]["model"]["unit_000"]["inlet_type"] = "PIECEWISE_CUBIC_POLY"

	model["root"]["input"]["model"]["unit_000"]["sec_000"] = OrderedDict()
	model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [50, 1, 1, 1]

	model["root"]["input"]["model"]["unit_000"]["sec_001"] = OrderedDict()
	model["root"]["input"]["model"]["unit_000"]["sec_001"]["const_coeff"] = [50, 0, 0, 0]

	model["root"]["input"]["model"]["unit_000"]["sec_002"] = OrderedDict()
	model["root"]["input"]["model"]["unit_000"]["sec_002"]["const_coeff"] = [100, 0, 0, 0]
	model["root"]["input"]["model"]["unit_000"]["sec_002"]["lin_coeff"] = [0.2, 0, 0, 0]


	# Set elements sequentially for unit_001
	model["root"]["input"]["model"]["unit_001"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["unit_type"] = "LUMPED_RATE_MODEL_WITH_PORES"
	model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
	model["root"]["input"]["model"]["unit_001"]["col_porosity"] = 0.37
	model["root"]["input"]["model"]["unit_001"]["col_dispersion"] = 1.0e-5
	model["root"]["input"]["model"]["unit_001"]["col_length"] = 0.014
	model["root"]["input"]["model"]["unit_001"]["velocity"] = 5.75e-4
	model["root"]["input"]["model"]["unit_001"]["par_porosity"] = 0.75
	model["root"]["input"]["model"]["unit_001"]["film_diffusion"] = [3.3e-3, 3.3e-3, 3.3e-3, 3.3e-3]
	model["root"]["input"]["model"]["unit_001"]["par_radius"] = 4.5e-5
	model["root"]["input"]["model"]["unit_001"]["adsorption_model"] = "STERIC_MASS_ACTION"

	model["root"]["input"]["model"]["unit_001"]["adsorption"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["is_kinetic"] = true
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["SMA_KA"] = [0.0, 35.5e-3, 1.59e-3, 7.70e-3]
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["SMA_KD"] = [1.0, 1, 1, 1]
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["SMA_LAMBDA"] = 1200.0
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["SMA_NU"] = [0.0, 4.7, 5.29, 3.7 ]
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["SMA_SIGMA"] = [0.0, 11.83, 10.6, 10.0]

	model["root"]["input"]["model"]["unit_001"]["init_c"] = [50.0, 0, 0, 0]
	model["root"]["input"]["model"]["unit_001"]["init_q"] = [1200.0, 0, 0, 0]

	model["root"]["input"]["model"]["unit_001"]["discretization"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["discretization"]["polyDeg"] = polyDeg
	model["root"]["input"]["model"]["unit_001"]["discretization"]["ncol"] = nCells
	model["root"]["input"]["model"]["unit_001"]["discretization"]["exact_integration"] = exactInt
	model["root"]["input"]["model"]["unit_001"]["discretization"]["nbound"] = ones(Bool, nComp)
	model["root"]["input"]["model"]["unit_001"]["discretization"]["use_analytic_jacobian"] = analJac

	# Set elements for unit_002
	model["root"]["input"]["model"]["unit_002"] = OrderedDict()
	model["root"]["input"]["model"]["unit_002"]["unit_type"] = "OUTLET"
	model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp


	# Set elements for solver
	model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
	model["root"]["input"]["solver"]["sections"]["nsec"] = 3
	model["root"]["input"]["solver"]["sections"]["section_times"] = [0.0, 10, 90, 1800]
	model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


	# Set elements for connections
	model["root"]["input"]["model"]["connections"] = OrderedDict()
	model["root"]["input"]["model"]["connections"]["nswitches"] = 1
	model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
	model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
	model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] = [0, 1, -1, -1, 5.75e-4, 
																				   1, 2, -1, -1, 5.75e-4]


	# Set elements for user_solution_times
	model["root"]["input"]["solver"]["user_solution_times"] = LinRange(0, 1800, 1800+1)

	# Set elements for time_integrator
	model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
	model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
	model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-10
	model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-10

	
	inlets, outlets, columns, switches, solverOptions = create_units(model)

	return inlets, outlets, columns, switches, solverOptions

end


# Evaluate the convergence using the evaluate_convergence function 
evaluate_convergence(QNDF(autodiff=false), c_analytical, nComp, nCell, polyDeg, 1, "LRMP", @__DIR__)

# Evaluate ODe solvers and save results in ODETests folder 
# evaluate_ODEsolvers(model_setup, c_analytical, nComp, nCell, polyDeg, 1, "LRMP", "ODETests//")
