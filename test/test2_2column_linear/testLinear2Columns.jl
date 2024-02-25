

# Add the include file custom to load packages and scripts. 
# the file is located on the main from which the file takes care of the rest. 
include(joinpath(@__DIR__,"..\\..\\include.jl"))



# Define the dictionary representing the model structure
nComp = 1
model = Dict(
    "root" => Dict(
        "input" => Dict(
            "model" => Dict(
				# Specify number of units 
                "nunits" => 4,

				# Inlet 
                "unit_000" => Dict(
                    "unit_type" => "INLET",
                    "ncomp" => nComp,
                    "inlet_type" => "PIECEWISE_CUBIC_POLY",
					# Inlet concentrations in the sections 
					"sec_000" => Dict(
						"const_coeff" => [1]
					),
					"sec_001" => Dict(
						"const_coeff" => [0]
					)
                ),

				# Column 1
                "unit_001" => Dict(
                    "unit_type" => "LUMPED_RATE_MODEL_WITHOUT_PORES_DG",
                    "ncomp" => nComp,
                    "col_porosity" => 0.6,
                    "col_dispersion" => 1e-4,
                    "col_length" => 1,
                    "velocity" => 2/60,
                    "adsorption_model" => "LINEAR",
                    "adsorption" => Dict(
                        "is_kinetic" => true,
                        "LIN_KA" => [1e8],
                        "LIN_KD" => [1e8]
                    ),
					# Initial conditions for column 1
                    "init_c" => [0],
                    "init_q" => [0],

					# Discretization 
					"discretization" => Dict(
						"polydeg" => 4,
						"ncol" => 16,
						"exact_integration" => 1,
						"nbound" => ones(Bool,nComp)
					)
                ),
				
				# Copy to column 2
				"unit_002" => Dict(
                    "unit_type" => "LUMPED_RATE_MODEL_WITHOUT_PORES_DG",
                    "ncomp" => nComp,
                    "col_porosity" => 0.6,
                    "col_dispersion" => 1e-4,
                    "col_length" => 1,
                    "velocity" => 2/60,
                    "adsorption_model" => "LINEAR",
                    "adsorption" => Dict(
                        "is_kinetic" => true,
                        "LIN_KA" => [1e8],
                        "LIN_KD" => [1e8]
                    ),
					# Initial conditions for column 1
                    "init_c" => [0],
                    "init_q" => [0],

					# Discretization 
					"discretization" => Dict(
						"polydeg" => 4,
						"ncol" => 16,
						"exact_integration" => 1,
						"nbound" => ones(Bool,nComp)
					)
                ),

				# Outlet 
				"unit_003" => Dict(
					"unit_type" => "OUTLET",
                    "ncomp" => nComp,
						),

				# Solver options for switches 	
				"solver" => Dict(
					"sections" => Dict(
						"nsec" => 4,
						"section_times" => [0.0, 60, 130, 190, 270],
						"section_continuity" => [0]
						)
				),

			# Sections and switches 
            "connections" => Dict(
                "nswitches" => 1,
                "switch_000" => Dict(
                    "section" => 0,
                    "connections" => [
                        0, 1, -1, -1, 2/60,  # [unit_000, unit_001, all components, all components, Q/ m^3*s^-1]
                        1, 2, -1, -1, 2/60   # [unit_001, unit_002, all components, all components, Q/ m^3*s^-1]
                    	]
					),
				"switch_001" => Dict(
                    "section" => 1,
                    "connections" => [
                        0, 1, -1, -1, 2/60,  # [unit_000, unit_001, all components, all components, Q/ m^3*s^-1]
                        1, 2, -1, -1, 2/60   # [unit_001, unit_002, all components, all components, Q/ m^3*s^-1]
                    	]
					)
				)
        	),
			"user_solution_times" => LinRange(0, 270, 270+1),
					"time_integrator" => Dict(
						"abstol" => 1e-12,
						"algtol" => 1e-10,
						"reltol" => 1e-10
			)
		)
	)
)
inlets, outlets, columns, switches, solverOptions = create_units(model)



solve_model(
			columns = columns,
			switches = switches,
			solverOptions = solverOptions, 
			outlets = outlets, # Defaults to (0,) as output is also written to units 
			alg = QNDF(autodiff=false), # Defaults to alg = QNDF(autodiff=false)
			)

using Plots
plot(columns[1].solution_outlet[:,1])
plot!(columns[2].solution_outlet[:,1])


# Compare to analytical solution 
using HDF5
filename = joinpath(@__DIR__,"model1.h5")
file = h5open(filename, "r")

C0 = vec(read(file["output/solution/unit_001"]["SOLUTION_OUTLET"]))
C1 = vec(read(file["output/solution/unit_002"]["SOLUTION_OUTLET"]))

plot(C0)
err = maximum([maximum(C0 - columns[2].solution_outlet[:,1]), maximum(maximum(C1 - columns[1].solution_outlet[:,1]))])
if err<1e-3
    println("Test succesful - error lower than 1e-3")
else
    println("Test unsuccesful - error larger than 1e-3")
end
