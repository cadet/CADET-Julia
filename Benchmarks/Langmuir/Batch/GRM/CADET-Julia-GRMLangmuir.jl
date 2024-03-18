#Run a convergence test where the results are compared to the semi-analytical solution.

# Add paths for the include file
include(joinpath(@__DIR__, fill("..", 4)..., "include.jl"))

# Specify number of cells, polynomial degree and number of components 
nCell =  [1,2,4,8]
polyDeg = [4]
polyDegPore = [4,5,6,8,10]
nComp = 2


# Load semi analytical solution
using CSV, DataFrames
c_analytical = CSV.read((joinpath(@__DIR__,"Semi-analytical_GRM_Langmuir.csv")),DataFrame)

function model_setup(nCells, polyDeg, polyDegPore, exactInt, analJac=false)

	nComp = 2
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
	model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [10, 10]
	model["root"]["input"]["model"]["unit_000"]["sec_001"] = OrderedDict()
	model["root"]["input"]["model"]["unit_000"]["sec_001"]["const_coeff"] = [0, 0]


	# Set elements sequentially for unit_001
	model["root"]["input"]["model"]["unit_001"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["unit_type"] = "GENERAL_RATE_MODEL"
	model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
	model["root"]["input"]["model"]["unit_001"]["col_porosity"] = 0.37
	model["root"]["input"]["model"]["unit_001"]["col_dispersion"] = 5.75e-8
	model["root"]["input"]["model"]["unit_001"]["col_length"] = 0.014
	model["root"]["input"]["model"]["unit_001"]["velocity"] = 5.75e-4
	model["root"]["input"]["model"]["unit_001"]["par_porosity"] = 0.75
	model["root"]["input"]["model"]["unit_001"]["par_radius"] = 4.5e-5
	model["root"]["input"]["model"]["unit_001"]["par_coreradius"] = 0
	model["root"]["input"]["model"]["unit_001"]["par_diffusion"] = [6.07e-11, 6.07e-11]
	model["root"]["input"]["model"]["unit_001"]["film_diffusion"] = [6.9e-6, 6.9e-6]
	model["root"]["input"]["model"]["unit_001"]["adsorption_model"] = "MULTI_COMPONENT_LANGMUIR"

	model["root"]["input"]["model"]["unit_001"]["adsorption"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["is_kinetic"] = false
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_KA"] = [0.1, 0.05] 
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_KD"] = [1.0, 1.0]
	model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_QMAX"] = [10.0, 10.0]

	model["root"]["input"]["model"]["unit_001"]["init_c"] = [0, 0]
	model["root"]["input"]["model"]["unit_001"]["init_q"] = [0, 0]

	model["root"]["input"]["model"]["unit_001"]["discretization"] = OrderedDict()
	model["root"]["input"]["model"]["unit_001"]["discretization"]["polyDeg"] = polyDeg
	model["root"]["input"]["model"]["unit_001"]["discretization"]["polyDegPore"] = polyDegPore
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
	model["root"]["input"]["solver"]["sections"]["section_times"] = [0.0, 10, 1500]
	model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


	# Set elements for connections
	model["root"]["input"]["model"]["connections"] = OrderedDict()
	model["root"]["input"]["model"]["connections"]["nswitches"] = 1
	model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
	model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
	model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] = [0, 1, -1, -1, 2/60, 
																				   1, 2, -1, -1, 2/60]


	# Set elements for user_solution_times
	model["root"]["input"]["solver"]["user_solution_times"] = LinRange(0, 1500, 1500+1)

	# Set elements for time_integrator
	model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
	model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
	model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-10
	model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-10


	inlets, outlets, columns, switches, solverOptions = create_units(model)
	return inlets, outlets, columns, switches, solverOptions 

end

# Evaluate the convergence using the evaluate_convergence function 
evaluate_convergence(QNDF(autodiff=false), c_analytical, nComp, nCell, polyDeg, polyDegPore, "GRM", @__DIR__,true)

# Evaluate ODe solvers and save results in ODETests folder 
[evaluate_ODEsolvers(c_analytical, nComp, nCell, [polyDeg[1]], 4, "GRM", joinpath(@__DIR__,"ODETests"), true) for _ in 1:2]