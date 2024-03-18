

# Add the include file custom to load packages and scripts. 
# the file is located on the main from which the file takes care of the rest. 
include(joinpath(@__DIR__,"..\\..\\include.jl"))



# Define the dictionary representing the model structure
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
model["root"]["input"]["model"]["unit_001"]["unit_type"] = "LUMPED_RATE_MODEL_WITH_PORES"
model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
model["root"]["input"]["model"]["unit_001"]["col_porosity"] = 0.6
model["root"]["input"]["model"]["unit_001"]["par_porosity"] = 0.2
model["root"]["input"]["model"]["unit_001"]["col_dispersion"] = 1e-5
model["root"]["input"]["model"]["unit_001"]["col_length"] = 0.25
model["root"]["input"]["model"]["unit_001"]["velocity"] = 2/60
model["root"]["input"]["model"]["unit_001"]["film_diffusion"] = [3.3e-3, 3.3e-3]
model["root"]["input"]["model"]["unit_001"]["par_radius"] = 1.0e-4
model["root"]["input"]["model"]["unit_001"]["adsorption_model"] = "MULTI_COMPONENT_LANGMUIR"

model["root"]["input"]["model"]["unit_001"]["adsorption"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["adsorption"]["is_kinetic"] = true
model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_KA"] = [0.1, 0.05] 
model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_KD"] = [1.0, 1.0]
model["root"]["input"]["model"]["unit_001"]["adsorption"]["MCL_QMAX"] = [10.0, 10.0]

model["root"]["input"]["model"]["unit_001"]["init_c"] = [0.0, 0.0]
model["root"]["input"]["model"]["unit_001"]["init_q"] = [0.0, 0.0]

model["root"]["input"]["model"]["unit_001"]["discretization"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["discretization"]["polyDeg"] = 4
model["root"]["input"]["model"]["unit_001"]["discretization"]["ncol"] = 32
model["root"]["input"]["model"]["unit_001"]["discretization"]["exact_integration"] = 1
model["root"]["input"]["model"]["unit_001"]["discretization"]["nbound"] = ones(Bool, nComp)

# Set elements for unit_002
model["root"]["input"]["model"]["unit_002"] = OrderedDict()
model["root"]["input"]["model"]["unit_002"]["unit_type"] = "OUTLET"
model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp


# Set elements for solver
model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
model["root"]["input"]["solver"]["sections"]["nsec"] = 2
model["root"]["input"]["solver"]["sections"]["section_times"] = [0.0, 12, 40]
model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


# Set elements for connections
model["root"]["input"]["model"]["connections"] = OrderedDict()
model["root"]["input"]["model"]["connections"]["nswitches"] = 1
model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] = [0, 1, -1, -1, 2/60, 
                                                                               1, 2, -1, -1, 2/60]


# Set elements for user_solution_times
model["root"]["input"]["solver"]["user_solution_times"] = LinRange(0, 40, 40*10+1)

# Set elements for time_integrator
model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-10
model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-10



inlets, outlets, columns, switches, solverOptions = create_units(model)

using Sundials
solve_model_dae(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, # Defaults to (0,) as output is also written to units 
			)

            
# Compare to analytical solution 
using CSV,DataFrames 
c_analytical = CSV.read((joinpath(@__DIR__,"test5_semi_analytical.csv")),DataFrame)
err = [0.0]
for i =1:nComp
    err[1] = maximum([err[1], maximum(abs.(columns[1].solution_outlet[:,i]-c_analytical[:,"C$(i-1)"]))])
end
if err[1]<1e-3
    println("Test succesful - error lower than 1e-3")
else
    println("Test unsuccesful - error larger than 1e-3")
end

