#Run a convergence test where the results are compared to the semi-analytical solution.


# Add paths for the include file
include(joinpath(@__DIR__, fill("..", 4)..., "include.jl"))

# Specify number of cells, polynomial degree and number of components 
nCell =  [2,4,8,16,32,64,128]
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

function model_setup(nCells, polyDeg, exactInt, analJac=false)

	nComp = 1
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
	model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [1]
	model["root"]["input"]["model"]["unit_000"]["sec_001"] = OrderedDict()
	model["root"]["input"]["model"]["unit_000"]["sec_001"]["const_coeff"] = [0]


	# Set elements sequentially for unit_001
	model["root"]["input"]["model"]["unit_001"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["unit_type"] = "LUMPED_RATE_MODEL_WITHOUT_PORES_DG"
	model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
	model["root"]["input"]["model"]["unit_001"]["col_porosity"] = 0.6
	model["root"]["input"]["model"]["unit_001"]["col_dispersion"] = 1e-4
	model["root"]["input"]["model"]["unit_001"]["col_length"] = 1
	model["root"]["input"]["model"]["unit_001"]["velocity"] = 2/60
	model["root"]["input"]["model"]["unit_001"]["adsorption_model"] = "LINEAR"

	model["root"]["input"]["model"]["unit_001"]["adsorption"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["is_kinetic"] = false
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["LIN_KA"] = [1.0]
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["LIN_KD"] = [1.0]

	model["root"]["input"]["model"]["unit_001"]["init_c"] = [0]
	model["root"]["input"]["model"]["unit_001"]["init_q"] = [0]

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
	model["root"]["input"]["solver"]["sections"]["nsec"] = 2
	model["root"]["input"]["solver"]["sections"]["section_times"] = [0.0, 60, 130]
	model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


	# Set elements for connections
	model["root"]["input"]["model"]["connections"] = OrderedDict()
	model["root"]["input"]["model"]["connections"]["nswitches"] = 1
	model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
	model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
	model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] = [0, 1, -1, -1, 2/60, 
																				   1, 2, -1, -1, 2/60]


	# Set elements for user_solution_times
	model["root"]["input"]["solver"]["user_solution_times"] = LinRange(0, 130, 130+1)

	# Set elements for time_integrator
	model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
	model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
	model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-10
	model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-10



	inlets, outlets, columns, switches, solverOptions = create_units(model)
	

	return inlets, outlets, columns, switches, solverOptions
end

# Evaluate the convergence using the evaluate_convergence function 
using BenchmarkTools
evaluate_convergence(model_setup, QNDF(autodiff=false), c_analytical, nComp, nCell, polyDeg, 1, "LRM", @__DIR__)

# Evaluate ODe solvers and save results in ODETests folder 
using Plots, Sundials
[evaluate_ODEsolvers(model_setup, c_analytical, nComp, nCell, [polyDeg[1]], 1, "LRM", joinpath(@__DIR__,"ODETests")) for _ in 1:2]