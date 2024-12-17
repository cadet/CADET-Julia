





# Solve the differential equations using the ODE solver
function solve_model_hybrid(; columns, switches::Switches, solverOptions, hybrid_model_setup, p_NN, outlets=(0,), alg=QNDF(autodiff=false), sensealg = nothing)
	"""
	Solves a hybrid model for all sections and stores the solution. 
	For now, having a hybrid model is treated as a special case but will be merged into the solve_model function later 
	
	"""

	# To have outlets as a tuple
	if typeof(columns)<:ModelBase
		columns = (columns,)
		if length(columns)>1
			throw("Multiple columns not supported for hybrid models. Yet.")
		end
	end
	
	# To have outlets as a tuple
	if typeof(outlets)==CreateOutlet
		outlets = (outlets,)
	end

	

	x0 = solverOptions.x0
	#running simulations
	for i = 1: length(switches.section_times) - 1 # corresponds to sections i=1

		# Set up parameter vector and empty elements 
		jacProto = nothing
		p_jac = nothing
		analytical_jac = nothing
		p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)

		# Update inlet concentrations at for new section for each unit. If they are static, they are not updated in RHS 
		@inbounds for sink = 1:solverOptions.nColumns 
			@inbounds for j=1:columns[sink].nComp
				inlet_concentrations!(columns[sink].cIn, switches, j, i, sink, solverOptions.x0, 0.0, DynamicInlets()) 
			end 
		end


		# If Analytical Jacobian == yes, set analytical Jacobian
		# Is only supported for batch operation! 
		if solverOptions.analyticalJacobian == true
		
			# determine static jacobian and allocation matrices that are stored in p_jac
			p_jac = jac_static(columns[1], switches.ConnectionInstance.u_tot[switches.switchSetup[i], 1], p) 
			p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
			analytical_jac = analytical_jac! #refers to the function
		end

		# If jacobian prototype, compute at every section time as switches might change Jacobian 
		if solverOptions.prototypeJacobian == true

			# determine jacobian prototype using finite differences - saves compilation time but perhaps change
			jacProto = sparse(jac_finite_diff(problem!,p,solverOptions.x0 .+ 1e-6, 1e-8))
			
			# set dq0dq to zero when computing SMA w. formulation 1 as dq/dt is not needed to be computed
			# This makes it slightly faster
			# if typeof(bind)==SMA 
			# 	@. @views jacProto[1 +columns[h].adsStride +columns[h].nComp*units[h].bindStride :columns[h].bindStride +columns[h].adsStride +columns[h].nComp*units[h].bindStride] = 0
			# end
		end

		# update the tspan and the inlets through i to the system
		tspan = (switches.section_times[i], switches.section_times[i+1])
		hybrid_model_setup.i = i

		fun = ODEFunction(hybrid_model_setup; jac_prototype = jacProto)
		prob = ODEProblem(fun, x0, (0, tspan[2]-tspan[1]), p_NN)
		idx_1 = findfirst(==(switches.section_times[i]), solverOptions.solution_times)
		idx_2 = findfirst(==(switches.section_times[i+1]), solverOptions.solution_times)
		sol_times = solverOptions.solution_times[idx_1 : idx_2] .- switches.section_times[i]
		sol = solve(prob, alg, saveat=sol_times, abstol=solverOptions.abstol, reltol=solverOptions.reltol, sensealg = sensealg) 
		
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
				if outlets[j].idx_outlet != [-1]
					for k = 1:columns[1].nComp
						if typeof(columns[outlets[j].idx_unit[switches.switchSetup[i]]]) == cstr
							outlets[j].solution_outlet[length(outlets[j].solution_times) + 1 : length(outlets[j].solution_times) + length(sol.t[2:end]), k] = sol(sol.t[2:end], idxs=k + outlets[j].idx_outlet[switches.switchSetup[i]]).u
						else
							outlets[j].solution_outlet[length(outlets[j].solution_times) + 1 : length(outlets[j].solution_times) + length(sol.t[2:end]), k] = sol(sol.t[2:end], idxs=k*columns[outlets[j].idx_unit[switches.switchSetup[i]]].ConvDispOpInstance.nPoints + outlets[j].idx_outlet[switches.switchSetup[i]]).u
						end
					end 
					append!(outlets[j].solution_times,sol.t[2:end] .+ tspan[1])
				end
			end
		end

		# Write to HDF5 using a function if relevant 
		
		
	end
    return nothing
end


# Define a function to compute the transport term for the LRM
function compute_transport(RHS, RHS_q, x, m::LRM, t, section, sink, switches, idx_units) 
	"""
	Transport for the LRM. 
	Allocating form. 
	
	"""

	# section = i from call 
	# sink is the unit i.e., h from previous call
    
	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	@inbounds for j = 1:m.nComp
		
		# Indices
		# For the indicies regarding mobile phase, + idx_units[sink] must be added to get the right column 
		# For the stationary phase, RHS_q is already a slice of the stationary phase of the right column
		m.idx =  1 + (j-1) * m.ConvDispOpInstance.nPoints : m.ConvDispOpInstance.nPoints + (j-1) * m.ConvDispOpInstance.nPoints

		# Determining inlet concentration 
		# inletConcentrations!(m.cIn, switches, j, switch, sink, x, t, idx_units) 
        cIn = ((switches.ConnectionInstance.cIn_c[section,sink, j] + 
				switches.ConnectionInstance.cIn_l[section,sink, j]*t +
				switches.ConnectionInstance.cIn_q[section,sink, j]*t^2 +
				switches.ConnectionInstance.cIn_cube[section,sink, j]*t^3) * switches.ConnectionInstance.u_inlet[switches.switchSetup[section], sink] +
				switches.ConnectionInstance.u_unit[switches.switchSetup[section], sink] * switches.ConnectionInstance.c_connect[switches.switchSetup[section], sink, j] * x[switches.ConnectionInstance.idx_connect[switches.switchSetup[section], sink, j]]) / switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink]


		# Convection Dispersion term
		cpp = (@view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp]) # mobile phase #
        Dh = ConvDispOperatorDGAlloc.residualImpl(cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], cIn, m.exactInt)


		# Mobile phase RHS 
		@. @views RHS[m.idx .+ idx_units[sink]] = Dh - m.Fc * RHS_q[m.idx]
	end
	
    nothing
end






