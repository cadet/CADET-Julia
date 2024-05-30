

# Add the include file custom to load packages and scripts. 
# the file is located on the main from which the file takes care of the rest. 
include(joinpath(@__DIR__,"..\\..\\include.jl"))



# Define the dictionary representing the model structure
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
model["root"]["input"]["model"]["unit_001"]["unit_type"] = "CSTR"
model["root"]["input"]["model"]["unit_001"]["volume"] = 0.25
model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
model["root"]["input"]["model"]["unit_001"]["init_c"] = 0



# Set elements sequentially for unit_002
model["root"]["input"]["model"]["unit_002"] = OrderedDict()
model["root"]["input"]["model"]["unit_002"]["unit_type"] = "LUMPED_RATE_MODEL_WITHOUT_PORES"
model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp
model["root"]["input"]["model"]["unit_002"]["col_porosity"] = 0.6
model["root"]["input"]["model"]["unit_002"]["col_dispersion"] = 1e-4
model["root"]["input"]["model"]["unit_002"]["col_length"] = 1
model["root"]["input"]["model"]["unit_002"]["velocity"] = 2/60
model["root"]["input"]["model"]["unit_002"]["adsorption_model"] = "LINEAR"

model["root"]["input"]["model"]["unit_002"]["adsorption"] = OrderedDict()
model["root"]["input"]["model"]["unit_002"]["adsorption"]["is_kinetic"] = false
model["root"]["input"]["model"]["unit_002"]["adsorption"]["LIN_KA"] = [1.0]
model["root"]["input"]["model"]["unit_002"]["adsorption"]["LIN_KD"] = [1.0]

model["root"]["input"]["model"]["unit_002"]["init_c"] = [0]
model["root"]["input"]["model"]["unit_002"]["init_q"] = [0]

model["root"]["input"]["model"]["unit_002"]["discretization"] = OrderedDict()
model["root"]["input"]["model"]["unit_002"]["discretization"]["polyDeg"] = 4
model["root"]["input"]["model"]["unit_002"]["discretization"]["ncol"] = 16
model["root"]["input"]["model"]["unit_002"]["discretization"]["exact_integration"] = 1
model["root"]["input"]["model"]["unit_002"]["discretization"]["nbound"] = ones(Bool, nComp)


# Set elements sequentially for unit_003
model["root"]["input"]["model"]["unit_003"] = OrderedDict()
model["root"]["input"]["model"]["unit_003"]["unit_type"] = "CSTR"
model["root"]["input"]["model"]["unit_003"]["volume"] = 0.25
model["root"]["input"]["model"]["unit_003"]["ncomp"] = nComp
model["root"]["input"]["model"]["unit_003"]["init_c"] = 0

# Set elements for unit_004
model["root"]["input"]["model"]["unit_004"] = OrderedDict()
model["root"]["input"]["model"]["unit_004"]["unit_type"] = "OUTLET"
model["root"]["input"]["model"]["unit_004"]["ncomp"] = nComp


# Set elements for solver
model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
model["root"]["input"]["solver"]["sections"]["nsec"] = 2
model["root"]["input"]["solver"]["sections"]["section_times"] = [0.0, 60, 180]
model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


# Set elements for connections
model["root"]["input"]["model"]["connections"] = OrderedDict()
model["root"]["input"]["model"]["connections"]["nswitches"] = 1
model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] = [0, 1, -1, -1, 2/60, 
                                                                               1, 2, -1, -1, 2/60,
                                                                               2, 3, -1, -1, 2/60,
                                                                               3, 4, -1, -1, 2/60,]


# Set elements for user_solution_times
model["root"]["input"]["solver"]["user_solution_times"] = LinRange(0, 180, 180+1)

# Set elements for time_integrator
model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-10
model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-10



inlets, outlets, columns, switches, solverOptions = create_units(model)


solve_model(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, # Defaults to (0,) as output is also written to units 
			alg = QNDF(autodiff=false), # Defaults to alg = QNDF(autodiff=false)
			)


# Compare to semi-analytical solution 
using CSV,DataFrames 
c_analytical = CSV.read((joinpath(@__DIR__,"test9.csv")),DataFrame)
err = maximum(abs.(outlets[1].solution_outlet[:,1]-c_analytical[:,"C0"]))
if err<1e-5
    println("Test succesful - error lower than 1e-5")
else
    println("Test unsuccesful - error larger than 1e-5")
end

