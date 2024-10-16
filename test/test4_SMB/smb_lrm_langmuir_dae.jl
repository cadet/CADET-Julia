
using Test

# Add the include file custom to load packages and scripts. 
# the file is located on the main from which the file takes care of the rest. 
include(joinpath(@__DIR__,"..\\..\\include.jl"))



# Define the dictionary representing the model structure

function SMB_model(Q2,Q4,QF,QD,QE,QR,ts)
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
    model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [10,10]


    model["root"]["input"]["model"]["unit_001"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["unit_type"] = "INLET"
    model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
    model["root"]["input"]["model"]["unit_001"]["inlet_type"] = "PIECEWISE_CUBIC_POLY"

    model["root"]["input"]["model"]["unit_001"]["sec_000"] = OrderedDict()
    model["root"]["input"]["model"]["unit_001"]["sec_000"]["const_coeff"] = [0,0]

    # Set elements for outlet 
    model["root"]["input"]["model"]["unit_002"] = OrderedDict()
    model["root"]["input"]["model"]["unit_002"]["unit_type"] = "OUTLET"
    model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp

    # Set elements for outlet
    model["root"]["input"]["model"]["unit_003"] = OrderedDict()
    model["root"]["input"]["model"]["unit_003"]["unit_type"] = "OUTLET"
    model["root"]["input"]["model"]["unit_003"]["ncomp"] = nComp

    # Set elements sequentially for column - unit_004
    model["root"]["input"]["model"]["unit_004"] = OrderedDict()
    model["root"]["input"]["model"]["unit_004"]["unit_type"] = "LUMPED_RATE_MODEL_WITHOUT_PORES"
    model["root"]["input"]["model"]["unit_004"]["ncomp"] = nComp
    model["root"]["input"]["model"]["unit_004"]["col_porosity"] = 0.4
    model["root"]["input"]["model"]["unit_004"]["col_dispersion"] = 1e-4
    model["root"]["input"]["model"]["unit_004"]["col_length"] = 1.0
    model["root"]["input"]["model"]["unit_004"]["cross_section_area"] = 1
    model["root"]["input"]["model"]["unit_004"]["adsorption_model"] = "MULTI_COMPONENT_LANGMUIR"

    model["root"]["input"]["model"]["unit_004"]["adsorption"] = OrderedDict()
    model["root"]["input"]["model"]["unit_004"]["adsorption"]["is_kinetic"] = true
    model["root"]["input"]["model"]["unit_004"]["adsorption"]["MCL_KA"] = [0.1, 0.05]
    model["root"]["input"]["model"]["unit_004"]["adsorption"]["MCL_KD"] = [1.0, 1.0]
    model["root"]["input"]["model"]["unit_004"]["adsorption"]["MCL_QMAX"] = [10.0, 10.0]

    model["root"]["input"]["model"]["unit_004"]["init_c"] = [0,0]
    model["root"]["input"]["model"]["unit_004"]["init_q"] = [0,0]

    model["root"]["input"]["model"]["unit_004"]["discretization"] = OrderedDict()
    model["root"]["input"]["model"]["unit_004"]["discretization"]["polyDeg"] = 4
    model["root"]["input"]["model"]["unit_004"]["discretization"]["ncol"] = 64
    model["root"]["input"]["model"]["unit_004"]["discretization"]["exact_integration"] = 1
    model["root"]["input"]["model"]["unit_004"]["discretization"]["nbound"] = ones(Bool, nComp)

    
    # Copy to remaining units
    model["root"]["input"]["model"]["unit_005"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
    model["root"]["input"]["model"]["unit_006"] = deepcopy(model["root"]["input"]["model"]["unit_004"])
    model["root"]["input"]["model"]["unit_007"] = deepcopy(model["root"]["input"]["model"]["unit_004"])

    # Set elements for solver
    n_cycles = 2
    switch_time = ts #s
    model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
    model["root"]["input"]["solver"]["sections"]["nsec"] = 4*n_cycles
    # Define switch times where cycle change
    section_times =  Float64[]
    push!(section_times, 0)
    for i in 0:n_cycles-1
        push!(section_times, Int64((4*i+1)*ts))
        push!(section_times, Int64((4*i+2)*ts))
        push!(section_times, Int64((4*i+3)*ts))
        push!(section_times, Int64((4*i+4)*ts))
    end

    model["root"]["input"]["solver"]["sections"]["section_times"] = section_times
    model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]


    # Set elements for connections
    model["root"]["input"]["model"]["connections"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["nswitches"] = 4
    model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
    model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] =[
        4, 5, -1, -1, Q4,#flowrates, Q, m3/s
        5, 6, -1, -1, Q4,
        6, 7, -1, -1, Q2,
        7, 4, -1, -1, Q2,
        0, 4, -1, -1, QF,
        1, 6, -1, -1, QD,
        4, 3, -1, -1, QR,
        6, 2, -1, -1, QE
    ]

    model["root"]["input"]["model"]["connections"]["switch_001"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_001"]["section"] = 1
    model["root"]["input"]["model"]["connections"]["switch_001"]["connections"] =[
        4,	5,	-1,	-1,	Q2,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q4,
        6,	7,	-1,	-1,	Q4,
        7,	4,	-1,	-1,	Q2,
        0,	5,	-1,	-1,	QF,
        1,	7,	-1,	-1,	QD,
        5,	3,	-1,	-1,	QR,
        7,	2,	-1,	-1,	QE
    ]

    model["root"]["input"]["model"]["connections"]["switch_002"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_002"]["section"] = 2
    model["root"]["input"]["model"]["connections"]["switch_002"]["connections"] =[
        4,	5,	-1,	-1,	Q2,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q2,
        6,	7,	-1,	-1,	Q4,
        7,	4,	-1,	-1,	Q4,
        0,	6,	-1,	-1,	QF,
        1,	4,	-1,	-1,	QD,
        6,	3,	-1,	-1,	QR,
        4,	2,	-1,	-1,	QE
    ]

    model["root"]["input"]["model"]["connections"]["switch_003"] = OrderedDict()
    model["root"]["input"]["model"]["connections"]["switch_003"]["section"] = 3
    model["root"]["input"]["model"]["connections"]["switch_003"]["connections"] =[
        4,	5,	-1,	-1,	Q4,#flowrates, Q, m3/s
        5,	6,	-1,	-1,	Q2,
        6,	7,	-1,	-1,	Q2,
        7,	4,	-1,	-1,	Q4,
        0,	7,	-1,	-1,	QF,
        1,	5,	-1,	-1,	QD,
        7,	3,	-1,	-1,	QR,
        5,	2,	-1,	-1,	QE
    ]


    # Set elements for user_solution_times
    model["root"]["input"]["solver"]["user_solution_times"] = LinRange(0, n_cycles*4*ts, n_cycles*4*ts*10+1)

    # Set elements for time_integrator
    model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
    model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-12
    model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-10
    model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-10

    return model 
end

# Determining flows 
CfA = 10                             #mol/m3
CfB = 10                             #mol/m3

ka = [0.1, 0.05]
kd = [1.0,1.0]
qmax = [10.0,10.0]


m1 = 1.2*ka[1]*qmax[1]
m2 = 1.2*ka[2]*qmax[2]
m3 = 0.9*ka[1]*qmax[1]
m4 = 0.8*ka[2]*qmax[2]
ColLength = 1.0
ts = 20 #s

eps = 0.4
Fc = (1-eps)/eps

Q1 = -(m1*eps - m1 - eps)*ColLength/ts
Q2 = -(m2*eps - m2 - eps)*ColLength/ts
Q3 = -(m3*eps - m3 - eps)*ColLength/ts
Q4 = -(m4*eps - m4 - eps)*ColLength/ts  
QD = -ColLength*(m1*eps - m4*eps - m1 + m4)/ts
QE = -ColLength*(m1*eps - m2*eps - m1 + m2)/ts
QF = ColLength*(m2*eps - m3*eps - m2 + m3)/ts
QR = -ColLength*(m3*eps - m4*eps - m3 + m4)/ts

# Generate dictionary 
model = SMB_model(Q2,Q4,QF,QD,QE,QR,ts)


# Create units for CADET-Julia 
inlets, outlets, columns, switches, solverOptions = create_units(model)


# Solve model 
using Sundials
solve_model_dae(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, # Defaults to (0,) as output is also written to units 
			)

# Compare to analytical solution 
using CSV,DataFrames 
c_analytical = CSV.read((joinpath(@__DIR__,"SMB_LRM_Langmuir_semi_analytical.csv")),DataFrame)

using Plots
plot(c_analytical[1:1601,"C0_E"])
plot!(c_analytical[1:1601,"C1_E"])
plot!(outlets[1].solution_outlet[:,1])
plot!(outlets[1].solution_outlet[:,2])

plot(c_analytical[1:1601,"C0_R"])
plot!(c_analytical[1:1601,"C1_R"])
plot!(outlets[2].solution_outlet[:,1])
plot!(outlets[2].solution_outlet[:,2])

err = [0.0]
for i =1:columns[1].nComp 
    err[1] = maximum([err[1], maximum(abs.(outlets[1].solution_outlet[1:1601,i]-c_analytical[1:1601,"C$(i-1)_E"]))])
    err[1] = maximum([err[1], maximum(abs.(outlets[2].solution_outlet[1:1601,i]-c_analytical[1:1601,"C$(i-1)_R"]))])
end

@test err[1] < 1e-3
