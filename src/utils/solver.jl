





# Solve the differential equations using the ODE solver
function solve_model(; columns, switches::Switches, solverOptions, outlets=(0,), alg=QNDF(autodiff=false))

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
		p_jac = jac_static(model)
		analytical_jac = analytical_jac! #refers to the function
	else
		# make p_jac and analytical_jac empty such that they are not used
		p_jac = nothing
		analytical_jac = nothing
	end

	x0 = solverOptions.x0
	#running simulations
	for i = 1: length(switches.section_times) - 1 # corresponds to sections i=1

		# If jacobian prototype, compute at every section time as switches might change Jacobian 
		if solverOptions.prototypeJacobian == true
			# set up a parameter vector 
			p = (columns, columns[1].RHS_q, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
			# determine jacobian prototype using finite differences - saves compilation time but perhaps change
			jacProto = sparse(jac_finite_diff(problem!,p,solverOptions.x0 .+ 1e-6, 1e-8))
			
			# set dq0dq to zero when computing SMA w. formulation 1 as dq/dt is not needed to be computed
			# This makes it slightly faster
			# if typeof(bind)==SMA 
			# 	@. @views jacProto[1 +columns[h].adsStride +columns[h].nComp*units[h].bindStride :columns[h].bindStride +columns[h].adsStride +columns[h].nComp*units[h].bindStride] = 0
			# end
		else
			# make jacProto empty such that it is not used
			jacProto = nothing
		end

		# update the tspan and the inlets through i to the system
		tspan = (switches.section_times[i], switches.section_times[i+1])
		p = (columns, columns[1].RHS_q, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
		fun = ODEFunction(problem!; jac_prototype = jacProto, jac = analytical_jac)
		prob = ODEProblem(fun, x0, (0, tspan[2]-tspan[1]), p)
		sol = solve(prob, alg, saveat=solverOptions.solution_times, abstol=solverOptions.abstol, reltol=solverOptions.reltol) #perhaps replace dt with solution_times
		
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

		# Write to HDF5 using a function if relevant 
		tt = columns[1].solution_outlet
		
	end
    return nothing
end


# Define the function representing the differential equations for transport and binding
function problem!(RHS, x, p, t)
    columns, RHS_q, i, nColumns, idx_units, switches = p
	# i corresponds to section 
	
	@inbounds for h = 1:nColumns 
		# Compute binding term. 
		# The cpp, qq and rhs_q are set up to ease code reading
		columns[h].cpp = @view x[1 + columns[h].adsStride + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h]]
		columns[h].qq = @view x[1 + columns[h].adsStride + columns[h].bindStride*columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		RHS_q = @view RHS[1 + columns[h].adsStride + columns[h].bindStride * columns[h].nComp + idx_units[h] : columns[h].adsStride + columns[h].bindStride * columns[h].nComp * 2 + idx_units[h]]
		compute_binding!(RHS_q, columns[h].cpp, columns[h].qq, columns[h].bind, columns[h].nComp, columns[h].bindStride, t)

		# Compute transport term
		compute_transport!(RHS, RHS_q, x, columns[h], t, i, h, switches, idx_units)

	end
	nothing
end



