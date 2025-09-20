"""
    solve_model_dae(; columns, switches::Switches, solverOptions, outlets=(0,), alg=IDA(init_all=false))

Solves the system of differential-algebraic equations (DAEs) for the specified model using the provided DAE solver algorithm.

# Arguments
- `columns`: Tuple or array of column unit objects to be simulated.
- `switches::Switches`: Switches object containing section times and connection information.
- `solverOptions`: SolverCache object with initial conditions, tolerances, and solver settings.
- `outlets`: Tuple of outlet unit objects (optional).
- `alg`: DAE solver algorithm to use (defaults to IDA).

# Details
- Iterates over all section times, updating parameters and solving the DAE system for each section.
- Handles both static and analytical Jacobians, as well as prototype Jacobians for efficiency.
- Updates solution matrices for columns and outlets after each section.
- Supports dynamic inlets and multiple unit types.

# Returns
Nothing. Results are stored in the `solution_outlet` and `solution_times` fields of the column and outlet objects.
"""
function solve_model_dae(; columns, switches::Switches, solverOptions, outlets=(0,), alg=IDA(init_all=false))
	# To have outlets as a tuple
	if typeof(columns)<:ModelBase
		columns = (columns,)
	end
	
	# To have outlets as a tuple
	if typeof(outlets)==CreateOutlet
		outlets = (outlets,)
	end

	x0 = solverOptions.x0
	differential_vars = ones(Int64,length(x0))
	dx0 = zeros(Float64,length(x0))
	
	#running simulations
	for i = 1: length(switches.section_times) - 1 # corresponds to sections i=1
	
		# Set up parameter vector and empty elements 
		jacProto = nothing
		p_jac = nothing
		analytical_jac = nothing
		fill!(dx0, 0.0)
		p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)

		# Update inlet concentrations at for new section for each unit. If they are static, they are not updated in RHS 
		@inbounds for sink = 1:solverOptions.nColumns 
			@inbounds for j=1:columns[sink].nComp
				inlet_concentrations!(columns[sink].cIn, switches, j, i, sink, solverOptions.x0, 0.0, DynamicInlets()) 
			end 
		end
		
		# If Analytical Jacobian == yes, set analytical Jacobian
		if solverOptions.analyticalJacobian == true
			# determine static jacobian and allocation matrices that are stored in p_jac
			p_jac = jac_static(columns[1],switches.ConnectionInstance.u_tot[switches.switchSetup[i], 1], p)
			p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
			analytical_jac = analytical_jac_dae! #refers to the function
		end

		# If jacobian prototype, compute at every section time as switches might change Jacobian 
		if solverOptions.prototypeJacobian == true
			
			# determine jacobian prototype
			jacProto = sparse(jac_finite_diff_dae(problemDAE!,p,solverOptions.x0 .+ 1e-6, 1e-8))
			#analytical_jac(jacProto, solverOptions.x0 .+ 1e-6, solverOptions.x0 .+ 1e-6, p, 1e-7, 0.0) #XXX: seems wrong; replaced by line above (adapted from solver.jl)
			
			# set dq0dq to zero when computing SMA w. formulation 1 as dq/dt is not needed to be computed
			# This makes it slightly faster
			# if typeof(bind)==SMA 
			# 	@. @views jacProto[1 +columns[h].adsStride +columns[h].nComp*units[h].bindStride :columns[h].bindStride +columns[h].adsStride +columns[h].nComp*units[h].bindStride] = 0
			# end
		else
			# make jacProto empty such that it is not used
			jacProto = nothing
		end
		
		# Update dx0 for the DAE for inlet 
		initialize_dae!(dx0, x0, p)

		# update the tspan and the inlets through i to the system
		tspan = (switches.section_times[i], switches.section_times[i+1])
		fun = DAEFunction(problemDAE!; jac = analytical_jac, jac_prototype = jacProto)
		prob = DAEProblem(fun, dx0, x0, (0.0, tspan[2]-tspan[1]), p, differential_vars = differential_vars)
		idx_1 = findfirst(==(switches.section_times[i]), solverOptions.solution_times)
		idx_2 = findfirst(==(switches.section_times[i+1]), solverOptions.solution_times)
		sol_times = solverOptions.solution_times[idx_1 : idx_2] .- switches.section_times[i]
		sol = solve(prob, alg, saveat=sol_times, abstol=solverOptions.abstol, reltol=solverOptions.reltol) 

		
		#New initial conditions
		x0 = sol.u[end]
		
		#Extract solution in solution_outlet in each unit 
		for j = 1: solverOptions.nColumns
			for k = 1:columns[j].nComp 
				if typeof(columns[j]) == cstr
					columns[j].solution_outlet[length(columns[j].solution_times) + 1 : length(columns[j].solution_times) + length(sol.t[2:end]), k] = sol(sol.t[2:end], idxs=k + solverOptions.idx_units[j]).u
				else
					columns[j].solution_outlet[length(columns[j].solution_times) + 1 : length(columns[j].solution_times) + length(sol.t[2:end]), k] = sol(sol.t[2:end], idxs=k*columns[j].ConvDispOpInstance.nPoints + solverOptions.idx_units[j]).u
				end				
			end 
			append!(columns[j].solution_times, sol.t[2:end] .+ tspan[1])
		end

		# Write outlets - if specified 
		if outlets != (0,) 
			for j in eachindex(outlets)
				if outlets[j].idx_outlet != [[-1]]
					for k in eachindex(outlets[j].idx_outlet[1])
						for l = 1:columns[1].nComp
							if typeof(columns[outlets[j].idx_unit[switches.switchSetup[i]][k]]) == cstr
								outlets[j].solution_outlet[length(outlets[j].solution_times) + 1 : length(outlets[j].solution_times) + length(sol.t[2:end]), l] += sol(sol.t[2:end], idxs=l + outlets[j].idx_outlet[switches.switchSetup[i]][k]).u .* outlets[j].u_outlet[switches.switchSetup[i]][k]
							else
								outlets[j].solution_outlet[length(outlets[j].solution_times) + 1 : length(outlets[j].solution_times) + length(sol.t[2:end]), l] += sol(sol.t[2:end], idxs=l*columns[outlets[j].idx_unit[switches.switchSetup[i]][k]].ConvDispOpInstance.nPoints + outlets[j].idx_outlet[switches.switchSetup[i]][k]).u .* outlets[j].u_outlet[switches.switchSetup[i]][k]
							end
						end 
					end
					# Divide outlet with total velocity to account for contribution for potential multiple sources for the outlet
					outlets[j].solution_outlet[length(outlets[j].solution_times) + 1 : length(outlets[j].solution_times) + length(sol.t[2:end]), :] ./= outlets[j].u_tot[switches.switchSetup[i]]

					append!(outlets[j].solution_times,sol.t[2:end] .+ tspan[1])
				end
			end
		end
		
		# Write to HDF5 using a function if relevant 
		
	end
    return nothing
end

"""
    problemDAE!(out, RHS, x, p, t)

Defines the system of DAEs for the model, assembling the residual vector for the DAE solver.

# Arguments
- `out`: Output vector for the DAE residuals.
- `RHS`: Vector for the right-hand side of the ODE part.
- `x`: Current state vector.
- `p`: Tuple of parameters, including columns, allocation vectors, section index, and switches.
- `t`: Current simulation time.

# Details
- Loops over all units (columns) and calls the appropriate transport and binding functions for each to fill the residual vector.
- Intended for use as the DAE function in the solver.

# Returns
Nothing. Modifies `out` in place.
"""
function problemDAE!(out, RHS, x, p, t)
    columns, RHS_q, cpp, qq, i, nColumns, idx_units, switches = p
	# i corresponds to section 
	
	@inbounds for h = 1:nColumns 

		# Compute binding term. 
		# The cpp, qq and rhs_q are set up to ease code reading
		cpp = @view x[1 + columns[h].adsStride + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h]]
		qq = @view x[1 + columns[h].adsStride + columns[h].bindStride*columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		compute_binding!(columns[h].RHS_q, cpp, qq, columns[h].bind, columns[h].nComp, columns[h].bindStride, t)

		# Compute transport term
		RHS_q = @view RHS[1 + columns[h].adsStride + columns[h].bindStride*columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		compute_transport!(out, RHS_q, cpp, x, columns[h], t, i, h, switches, idx_units)
		
		# Subtract mobile phase RHS
		columns[1].idx = 1 + idx_units[h] : idx_units[h] + columns[h].adsStride + columns[h].bindStride * columns[h].nComp 
		@. @views out[columns[1].idx] -= RHS[columns[1].idx]

		
		#Add the isotherm part to the DAE - now steady
		@. @views out[1 + columns[1].adsStride + columns[1].bindStride*columns[1].nComp + idx_units[h] : columns[1].adsStride + 2 * columns[1].bindStride * columns[1].nComp + idx_units[h]] = columns[h].RHS_q

	end
	nothing
end


"""
    initialize_dae!(dx0, x0, p)

Initializes the DAE system by computing consistent initial derivatives for the state variables.

# Arguments
- `dx0`: Output vector for the initial derivatives.
- `x0`: Initial state vector.
- `p`: Tuple of parameters for the model.

# Returns
Nothing. Modifies `dx0` in place.
"""
function initialize_dae!(dx0, x0, p)
	columns, RHS_q, cpp, qq, i, nColumns, idx_units, switches = p

    # Update dx0 for the DAE for inlet 
	problemDAE!(dx0, zeros(Float64,length(x0)), x0,p,0.0)

	# Provide consistent initialization using the ODE formulation 
	RHS_ODE = zeros(Float64,length(x0))
	problem!(RHS_ODE, x0,p,0.0)
	
	# insert into dx0 vector to ensure consistent initialization 
	# Hence insert RHS_ODE into pore and stationary phase of dx0
	for j = 1:nColumns 
		columns[1].idx = 1 + idx_units[j] + columns[j].ConvDispOpInstance.nPoints * columns[j].nComp : idx_units[j] + columns[j].adsStride + 2 *columns[j].bindStride * columns[j].nComp 
		dx0[columns[1].idx] = RHS_ODE[columns[1].idx]
	end
	nothing
end
