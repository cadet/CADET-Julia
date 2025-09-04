using Test, CADETJulia, CSV, DataFrames

"""
This test is testing an acyclic model with four columns which gatheres the output from the two last columns in the outlet. 
The output from the columns are crossed such that the columns have multiple units as input, see flowsheet in folder. 

The test and its parameters are taken from:

Leweke, S., 2021. 
Unified Modeling and Efficient Simulation of Chromatographic Separation Processes (Ph.D. thesis). 
RWTH Aachen University, Aachen 
http://dx.doi.org/10.18154/RWTH-2022-01184


The test is also shown in the publication: 
Frandsen, J., Breuer, J. M., Schmölder, J., Abildskov, J., Huusom, J. K., Gernaey, K. V., & von Lieres, E. (2025). 
High-Performance C++ and Julia solvers in CADET for weakly and strongly coupled continuous chromatography problems. 
Computers and Chemical Engineering, 202. 
https://doi.org/10.1016/j.compchemeng.2025.109295

"""

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
model["root"]["input"]["model"]["unit_000"]["sec_002"] = OrderedDict()
model["root"]["input"]["model"]["unit_000"]["sec_002"]["const_coeff"] = [0] 


model["root"]["input"]["model"]["unit_001"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["unit_type"] = "INLET"
model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
model["root"]["input"]["model"]["unit_001"]["inlet_type"] = "PIECEWISE_CUBIC_POLY"

model["root"]["input"]["model"]["unit_001"]["sec_000"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["sec_000"]["const_coeff"] = [1]
model["root"]["input"]["model"]["unit_001"]["sec_001"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["sec_001"]["const_coeff"] = [5]
model["root"]["input"]["model"]["unit_001"]["sec_002"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["sec_002"]["const_coeff"] = [0]

# Set elements for outlet 
model["root"]["input"]["model"]["unit_006"] = OrderedDict()
model["root"]["input"]["model"]["unit_006"]["unit_type"] = "OUTLET"
model["root"]["input"]["model"]["unit_006"]["ncomp"] = nComp



# Set elements sequentially for unit_001
model["root"]["input"]["model"]["unit_002"] = OrderedDict()
model["root"]["input"]["model"]["unit_002"]["unit_type"] = "LUMPED_RATE_MODEL_WITH_PORES"
model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp
model["root"]["input"]["model"]["unit_002"]["col_porosity"] = 0.37
model["root"]["input"]["model"]["unit_002"]["col_dispersion"] = 2e-7
model["root"]["input"]["model"]["unit_002"]["col_length"] = 1.4e-2
model["root"]["input"]["model"]["unit_002"]["cross_section_area"] = 1
model["root"]["input"]["model"]["unit_002"]["par_porosity"] = 0.75
model["root"]["input"]["model"]["unit_002"]["film_diffusion"] = [6.9e-6]
model["root"]["input"]["model"]["unit_002"]["par_radius"] = 45e-6
LRMP_Q3 = 3.45*1e-2 / 60 * 0.37
model["root"]["input"]["model"]["unit_002"]["adsorption_model"] = "LINEAR"

model["root"]["input"]["model"]["unit_002"]["adsorption"] = OrderedDict()
model["root"]["input"]["model"]["unit_002"]["adsorption"]["is_kinetic"] = true
model["root"]["input"]["model"]["unit_002"]["adsorption"]["LIN_KA"] = [3.55] 
model["root"]["input"]["model"]["unit_002"]["adsorption"]["LIN_KD"] = [0.1]


model["root"]["input"]["model"]["unit_002"]["init_c"] = [0]
model["root"]["input"]["model"]["unit_002"]["init_q"] = [0]

model["root"]["input"]["model"]["unit_002"]["discretization"] = OrderedDict()
model["root"]["input"]["model"]["unit_002"]["discretization"]["polyDeg"] = 3
model["root"]["input"]["model"]["unit_002"]["discretization"]["ncol"] = 8
model["root"]["input"]["model"]["unit_002"]["discretization"]["exact_integration"] = true
model["root"]["input"]["model"]["unit_002"]["discretization"]["nbound"] = ones(Bool, nComp)


# Copy to remaining units
model["root"]["input"]["model"]["unit_003"] = deepcopy(model["root"]["input"]["model"]["unit_002"])
model["root"]["input"]["model"]["unit_004"] = deepcopy(model["root"]["input"]["model"]["unit_002"])
model["root"]["input"]["model"]["unit_005"] = deepcopy(model["root"]["input"]["model"]["unit_002"])


# Unit LRMP4 
model["root"]["input"]["model"]["unit_003"]["col_length"] = 4.2e-2 
model["root"]["input"]["model"]["unit_003"]["adsorption"]["is_kinetic"] = false
model["root"]["input"]["model"]["unit_003"]["adsorption"]["LIN_KA"] = [3.55] 
model["root"]["input"]["model"]["unit_003"]["adsorption"]["LIN_KD"] = [0.1]


# Unit LRMP5 
model["root"]["input"]["model"]["unit_004"]["adsorption"]["is_kinetic"] = false
model["root"]["input"]["model"]["unit_004"]["adsorption"]["LIN_KA"] = [21.4286] 
model["root"]["input"]["model"]["unit_004"]["adsorption"]["LIN_KD"] = [1.]


# Unit LRMP6
model["root"]["input"]["model"]["unit_005"]["adsorption"]["LIN_KA"] = [4.55] 
model["root"]["input"]["model"]["unit_005"]["adsorption"]["LIN_KD"] = [0.12]


#solution times
model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
model["root"]["input"]["solver"]["sections"]["nsec"] = 3
model["root"]["input"]["solver"]["sections"]["section_times"] = [0., 250, 300, 3000]  


# Set elements for connections
model["root"]["input"]["model"]["connections"] = OrderedDict()
model["root"]["input"]["model"]["connections"]["nswitches"] = 1

model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] =[
    0, 2, -1, -1, LRMP_Q3, 
    2, 4, -1, -1, LRMP_Q3/2,
    2, 5, -1, -1, LRMP_Q3/2,
    1, 3, -1, -1, LRMP_Q3,
    3, 4, -1, -1, LRMP_Q3/2,
    3, 5, -1, -1, LRMP_Q3/2,
    4, 6, -1, -1, LRMP_Q3,
    5, 6, -1, -1, LRMP_Q3,
]

# Set elements for user_solution_times
model["root"]["input"]["solver"]["user_solution_times"] = collect(0: 1.0: 3000) 

# Set elements for time_integrator
model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-12
model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-12

# Create units 
inlets, outlets, columns, switches, solverOptions = create_units(model)

# Solve model 
solve_model(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, 
			)

# Compare to analytical solution 
test11 = CSV.read((joinpath(@__DIR__,"test11.csv")),DataFrame)
err = maximum(abs.(outlets[1].solution_outlet[:,1]-test11[:,"C0"]))

@test err < 1e-5