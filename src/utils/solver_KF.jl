include("update_step_solver.jl")


function systemJacobian!(J,x,p,t)
	""" 
	Calculates the Jacobian of the entire system by combining the Jacobians of each unit
	"""
	columns, RHS_q, cpp, qq, i, nColumns, idx_units, switches = p

	#initialize Jacobian
	fill!(J, 0.0)

	#Calculate Jacobian for each unit 

	@inbounds for h = 1:nColumns
		unit = columns[h]

		if typeof(unit) == cstr

			J_unit = jacobian(columns[h],t, 0.0, i, h, switches)
			start_idx = idx_units[h]
			
			@view(J[(start_idx+1):(start_idx+unit.nComp), (start_idx+1):(start_idx+unit.nComp)]) .= J_unit

		else
			println("Only CSTRs are supported")
		end 
	end

	return nothing
	
end


# Solve the differential equations using the ODE solver
function solve_model(; columns, switches::Switches, solverOptions, outlets=(0,), alg=QNDF(autodiff=false))
	
	# Set up parameter tuple
	jacProto = nothing
	p_jac = nothing
	analytical_jac = nothing
	RHS_q = 0.
	cpp = 0.
	qq = 0.

	# find columns to extract allocation vectors RHS_q, cpp and qq.
	# If there are no columns, the allocations vectors are unnecessary and kept at 0
	for j in eachindex(columns)
		if typeof(columns[j]) <: ModelBase
			RHS_q = columns[j].RHS_q
			cpp = columns[j].cpp
			qq = columns[j].qq
			break
		end 
	end

	x0 = solverOptions.x0
	
	# save measurements if KalmanFilter is used
	original_measurements = nothing
    if isa(alg, KalmanFilter)
        original_measurements = deepcopy(alg.measurement_buffer.measurements)
    end

	#running simulations
	for i = 1: length(switches.section_times) - 1 # corresponds to sections i=1

		# Update i in parameter tuple
		p = (columns, RHS_q, cpp, qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)

		# Update inlet concentrations at for new section for each unit. If they are static, they are not updated in RHS 
		@inbounds for sink = 1:solverOptions.nColumns 
			@inbounds for j=1:columns[sink].nComp
				#switches.inlet_conditions[section, sink, k]
				inlet_concentrations!(columns[sink].cIn, switches, j, i, sink, solverOptions.x0, 0.0, DynamicInlets()) 
			end 
		end

		# If Analytical Jacobian == yes, set analytical Jacobian
		# Is only supported for batch operation! 
		solverOptions.analyticalJacobian = true # for debugging purposes
		if solverOptions.analyticalJacobian == true
		
			# determine static jacobian and allocation matrices that are stored in p_jac
			#p_jac = jac_static(columns[1], switches.ConnectionInstance.u_tot[switches.switchSetup[i], 1], p) 
			#p = (columns, columns[1].RHS_q, columns[1].cpp, columns[1].qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, p_jac)
			#analytical_jac = analytical_jac! #refers to the function

			analytical_jac = (J,x,p,t) -> systemJacobian!(J,x,p,t)
		end

		# If jacobian prototype, compute at every section time as switches might change Jacobian 
		if solverOptions.prototypeJacobian == false

			# determine jacobian prototype using finite differences - saves compilation time but perhaps change
			jacProto = sparse(jac_finite_diff(systemProblem!,p,solverOptions.x0 .+ 1e-6, 1e-8))
			
			# set dq0dq to zero when computing SMA w. formulation 1 as dq/dt is not needed to be computed
			# This makes it slightly faster
			# if typeof(bind)==SMA 
			# 	@. @views jacProto[1 +columns[h].adsStride +columns[h].nComp*units[h].bindStride :columns[h].bindStride +columns[h].adsStride +columns[h].nComp*units[h].bindStride] = 0
			# end
		end

		if isa(alg, KalmanFilter)

			empty!(alg.measurement_buffer.measurements)

			section_start = switches.section_times[i]
			section_end = switches.section_times[i+1]

			for measurement in original_measurements
				if section_start <= measurement.time <= section_end
					local_time = measurement.time - section_start

					add_measurement!(alg, measurement.unit_idx, local_time, measurement.value)
				end
			end
			
			alg.last_update_time = -Inf

		end

		# update the tspan and the inlets through i to the system
		tspan = (switches.section_times[i], switches.section_times[i+1])
		fun = ODEFunction(systemProblem!; jac = analytical_jac)
		prob = ODEProblem(fun, x0, (0, tspan[2]-tspan[1]), p)
		idx_1 = findfirst(==(switches.section_times[i]), solverOptions.solution_times)
		idx_2 = findfirst(==(switches.section_times[i+1]), solverOptions.solution_times)
		sol_times = solverOptions.solution_times[idx_1 : idx_2] .- switches.section_times[i]
		sol = solve(prob, alg, saveat=sol_times, abstol=solverOptions.abstol, reltol=solverOptions.reltol) 
		
		#New initial conditions
		x0 = sol.u[end]
		#Extract solution in solution_outlet in each unit 
		for j = 1: solverOptions.nColumns
			# todo: understand wht solver changes the size of the solution_outlet matrix
			num_new_points = length(sol.t[2:end])
			current_size = size(columns[j].solution_outlet, 1)
			new_size_needed = length(columns[j].solution_times) + num_new_points
			
			if current_size != new_size_needed
				# Create new matrix of exact needed size
				new_matrix = similar(columns[j].solution_outlet, new_size_needed, size(columns[j].solution_outlet, 2))
				# Copy existing data up to the minimum of both sizes
				copy_size = min(current_size, new_size_needed)
				new_matrix[1:copy_size, :] = view(columns[j].solution_outlet, 1:copy_size, :)
				columns[j].solution_outlet = new_matrix
			end
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
				# Get the number of new points to add
				num_new_points = length(sol.t[2:end])

				# Check and resize the outlets solution matrix if needed
				current_size = size(outlets[j].solution_outlet, 1) 
				new_size_needed = length(outlets[j].solution_times) + num_new_points
				
				if current_size != new_size_needed
					# Create new matrix of exact needed size
					new_matrix = similar(outlets[j].solution_outlet, new_size_needed, size(outlets[j].solution_outlet, 2))
					# Copy existing data up to the minimum of both sizes
					copy_size = min(current_size, new_size_needed)
					new_matrix[1:copy_size, :] = view(outlets[j].solution_outlet, 1:copy_size, :)
					outlets[j].solution_outlet = new_matrix
				end
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

	end
			# Write to HDF5 using a function if relevant 
	if isa(alg, KalmanFilter) && !isnothing(original_measurements)
		alg.measurement_buffer.measurements = original_measurements
	end
    return nothing
end


# Define the function representing the differential equations
function systemProblem!(RHS, x, p, t)
    columns, RHS_q, cpp, qq, i, nColumns, idx_units, switches = p
	# i corresponds to section 
	
	@inbounds for h = 1:nColumns 
		compute!(RHS, RHS_q, cpp, qq, x, columns[h], t, i, h, switches, idx_units) 
	end
	nothing
end

# Compute column model RHS
function compute!(RHS, RHS_q, cpp, qq, x, m :: ModelBase, t, section, sink, switches, idx_units) 
	# section = i from call 
	# sink is the unit i.e., h from previous call
	
	# Compute binding term. 
	# The cpp, qq and rhs_q are set up to ease code reading
	cpp = @view x[1 + m.adsStride + idx_units[sink] : m.adsStride + m.bindStride * m.nComp + idx_units[sink]]
	qq = @view x[1 + m.adsStride + m.bindStride*m.nComp + idx_units[sink] : m.adsStride + m.bindStride * m.nComp * 2 + idx_units[sink]]
	RHS_q = @view RHS[1 + m.adsStride + m.bindStride * m.nComp + idx_units[sink] : m.adsStride + m.bindStride * m.nComp * 2 + idx_units[sink]]
	compute_binding!(RHS_q, cpp, qq, m.bind, m.nComp, m.bindStride, t)

	# Compute transport term
	compute_transport!(RHS, RHS_q, cpp, x, m, t, section, sink, switches, idx_units)
    
    nothing
end




