





# Solve the differential equations using the ODE solver
function solve_model_dae(; columns, switches::Switches, solverOptions, outlets=(0,), alg=IDA(init_all=false))

	# To have outlets as a tuple
	if typeof(columns)<:ModelBase
		columns = (columns,)
	end
	
	# To have outlets as a tuple
	if typeof(outlets)==CreateOutlet
		outlets = (outlets,)
	end

	# If Analytical Jacobian == yes, set analytical Jacobian
	if solverOptions.analyticalJacobian == true
		# determine static jacobian and allocation matrices that are stored in p_jac
		p_jac = jac_static(columns[1])
		analytical_jac = analytical_jac! #refers to the function
	else
		# make p_jac and analytical_jac empty such that they are not used
		p_jac = nothing
		analytical_jac = nothing
	end
	
	x0 = solverOptions.x0
	differential_vars = ones(Int64,length(x0))
	dx0 = zeros(Float64,length(x0))
	
	#running simulations
	for i = 1: length(switches.section_times) - 1 # corresponds to sections i=1

		# If jacobian prototype, compute at every section time as switches might change Jacobian 
		if solverOptions.prototypeJacobian == true
			# set up a parameter vector 
			p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
			# determine jacobian prototype using finite differences - saves compilation time but perhaps change
			jacProto = sparse(jac_finite_diff_dae(problemDAE!,p,solverOptions.x0 .+ 1e-6, 1e-8))
			
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
		fill!(dx0,0.0)
		p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
		initialize_dae!(dx0, x0, p)
		
		# update the tspan and the inlets through i to the system
		tspan = (switches.section_times[i], switches.section_times[i+1])
		fun = DAEFunction(problemDAE!; jac_prototype = jacProto ) #jac_prototype = jacProto, jac = analytical_jac
		prob = DAEProblem(fun, dx0, x0, (0.0, tspan[2]-tspan[1]), p, differential_vars = differential_vars)
		sol = solve(prob, alg, saveat=solverOptions.solution_times, abstol=solverOptions.abstol, reltol=solverOptions.reltol) 
		
		#New initial conditions
		x0 = sol.u[end]
		
		#Extract solution in solution_outlet in each unit 
		for j = 1: solverOptions.nColumns
			for k = 1:columns[j].nComp 
				
				columns[j].solution_outlet[length(columns[j].solution_times) + 1 : length(columns[j].solution_times) + length(sol.t[2:end]), k] = sol(sol.t[2:end], idxs=k*columns[j].ConvDispOpInstance.nPoints + solverOptions.idx_units[j]).u
			end 
			append!(columns[j].solution_times, sol.t[2:end] .+ tspan[1])
		end

		# Write outlets - if specified 
		if outlets != (0,) 
			for j in eachindex(outlets)
				if outlets[j].idx_outlet != [-1]
					for k = 1:columns[1].nComp
						
						outlets[j].solution_outlet[length(outlets[j].solution_times) + 1 : length(outlets[j].solution_times) + length(sol.t[2:end]), k] = sol(sol.t[2:end], idxs=k*columns[outlets[j].idx_unit[switches.switchSetup[i]]].ConvDispOpInstance.nPoints + outlets[j].idx_outlet[switches.switchSetup[i]]).u
					end 
					append!(outlets[j].solution_times,sol.t[2:end] .+ tspan[1])
				end
			end
		end
		xx = outlets[1].solution_outlet
		# Write to HDF5 using a function if relevant 
		
	end
    return nothing
end


# Define the function representing the differential equations for transport and binding
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


# Define the function representing the differential equations for transport and binding
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



