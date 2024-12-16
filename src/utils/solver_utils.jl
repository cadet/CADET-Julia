



mutable struct SolverCache
	# Here, initial conditions and solver options are specified. 

	# Initial conditions: can be a vector of all x0, could be a single x0 or could be different x0 for every component
	x0::Vector{Float64} #if x0 is given as input, defaults to 0
	
	# tolerances for the solver, defaults to 1e-12, 1e-10
	abstol::Float64
	reltol::Float64

	# time step for solution
	dt::Float64
	solution_times::Vector{Float64}

	# solver options
	prototypeJacobian::Bool # Generate a prototype jacobian - makes it run a lot faster
	analyticalJacobian::Bool # Use analytical jacobian - does not necessary improve performance
	
	idx_units::Vector{Int64}
	nColumns::Int64
	# the ODE solver cannot be stored because when changing the solver, the datatype is also changed..

	function SolverCache(; columns::Union{Tuple,ModelBase}, switches::Switches, outlets::Union{Tuple,CreateOutlet} = (0,), x0 = [0], abstol=1e-12, reltol=1e-10, dt=1.0, solution_times=[0], prototypeJacobian=true, analyticalJacobian=false)

		# To have outlets as a tuple
		if typeof(columns)<:ModelBase
			columns = (columns,)
		end

		# To have outlets as a tuple
		if typeof(outlets)==CreateOutlet
			outlets = (outlets,)
		end

		nColumns = length(columns)

		# Determining the indices for when a new unit is starting and setting up the solution_outlet matrices
		idx_units = zeros(Int64, nColumns)
		for i = 2:nColumns 
			idx_units[i] = idx_units[i-1] + columns[i-1].unitStride
		end
		
		# Solution times 
		if solution_times != [0] && length(solution_times) >= 2
			dt = solution_times[2] - solution_times[1]
		elseif solution_times == [0] # if solution times are not specified, determine based on dt 
			solution_times = 0:dt:switches.section_times[end]
		end

		x0len = 0
		for i = 1:nColumns
			x0len += columns[i].unitStride
			
			# Constructing the solution outlet for each unit 
			columns[i].solution_outlet = -ones(Float64,length(solution_times),columns[i].nComp)
			columns[i].solution_times = Float64[]
			#columns[i].solution_outlet[:,end] = solverOptions.solution_times #should only be added once solved 
			
			# Attaching the correct index to the outlet of a unit 
			if outlets == (0,)
				continue
			else 
				for j in eachindex(outlets)
					for k in eachindex(outlets[j].idx_unit)
						if outlets[j].idx_unit[k] == i # if connection between column and outlet 
							outlets[j].idx_outlet[k] = idx_units[i] # set the right index 
						end 
					end
					outlets[j].solution_outlet = -ones(Float64,length(solution_times),columns[i].nComp)
				end
			end
		end

		# Set up the initial conditions, x0 vector, correctly
		if x0 == [0]
			# if x0 is defaulted, take the values for c0, cp0 and q0 for each unit 
			x0 = zeros(Float64,x0len)

			for i = 1: nColumns
			
				if typeof(columns[i]) == cstr
					x0[1 + idx_units[i] : idx_units[i] + columns[i].unitStride] = columns[i].c0
				else
					x0[1 + idx_units[i] : idx_units[i] + columns[i].ConvDispOpInstance.nPoints * columns[i].nComp] = columns[i].c0 
					x0[1 + idx_units[i] + columns[i].adsStride : idx_units[i] + columns[i].adsStride + columns[i].bindStride * columns[i].nComp] = columns[i].cp0 
					x0[1 + idx_units[i] + columns[i].adsStride + columns[i].bindStride*columns[i].nComp : idx_units[i] + columns[i].adsStride + columns[i].bindStride * columns[i].nComp * 2] = columns[i].q0
				end
			end 
		elseif length(x0) == x0len # if the whole x0 is specified 
			nothing 
		else 
			throw("x0 not correctly specified. The length must match the length of all the units initial conditions. If not specifying, the initial conditions from each unit is used.")
	
		end 

		# Setting up the switches such that they follow a repeated pattern. 
		# that means if the number of section times are longer than the number of switches, the switches should be repeated 
		# such that one does not need to specify all the switches in for example in an SMB but only one cycle. 
		# it is done so in a row-wise manner such that if a switch is not determined at a specific section time, it will be overwritten by the cycle 
		# it checks if u_tot = -1 (default value), then that will be the last specified switch and the remaining will be repeated. 
		# Fill in elements in between the non-filled values of the switchSetup 
		# This means replace the -1 
		switches.switchSetup = rearrange_switch_setup(switches)
		
		# If having multiple switches, 
		if switches.nSwitches>1
			# Repeat the concentration specifications 
			switches.ConnectionInstance.cIn_c = repeat_pattern(switches.ConnectionInstance.cIn_c, switches.switchSetup, switches.nSwitches)
			switches.ConnectionInstance.cIn_l = repeat_pattern(switches.ConnectionInstance.cIn_l, switches.switchSetup, switches.nSwitches)
			switches.ConnectionInstance.cIn_q = repeat_pattern(switches.ConnectionInstance.cIn_q, switches.switchSetup, switches.nSwitches)
			switches.ConnectionInstance.cIn_cube = repeat_pattern(switches.ConnectionInstance.cIn_cube, switches.switchSetup, switches.nSwitches)
		end
		
		
		# Fill in initial conditions in solution_outlet matrices
		for i=1:nColumns
			for j=1:columns[1].nComp # Storing initial conditions in output matrix
				if typeof(columns[i]) == cstr
					columns[i].solution_outlet[1,j] = x0[j + idx_units[i]]
				else
					columns[i].solution_outlet[1,j] = x0[j*columns[i].ConvDispOpInstance.nPoints + idx_units[i]]
				end
			end
			append!(columns[i].solution_times, 0)
		end

		# Fill in initial conditions in outlets 
		if outlets != (0,)
			for i in eachindex(outlets)
				if outlets[i].idx_outlet != [-1]
					for j=1:columns[1].nComp # Storing initial conditions in output matrix
						if typeof(columns[outlets[i].idx_unit[switches.switchSetup[i]]]) == cstr
							outlets[i].solution_outlet[1,j] = x0[j + idx_units[i]]
						else
							outlets[i].solution_outlet[1,j] = x0[j*columns[outlets[i].idx_unit[switches.switchSetup[1]]].ConvDispOpInstance.nPoints + outlets[i].idx_outlet[switches.switchSetup[i]]]
						end
					end
					append!(outlets[i].solution_times, 0)
				end
			end
		end

		
		new(x0, abstol, reltol, dt, solution_times, prototypeJacobian, analyticalJacobian, idx_units, nColumns)
	end
end

function rearrange_switch_setup(switches)
	"""
		A function to rearrange the switchSetup such that it follows a repeated pattern.
		Inputs are:
		switches: The switches object

		Outputs are:
		a: The rearranged switchSetup

		If there is only one switch, it copy the switches for the number of section times.
	"""
	a = -ones(Int64,length(switches.switchSetup))

	#if there is only one switch
	if switches.nSwitches==1 
		a[:] .= switches.switchSetup[1]
	#if there are two switches
	elseif switches.nSwitches>1 #
		# The inlet concentrations are following a repetetive pattern 
		a = repeat_pattern(switches.switchSetup)
	end
	return a
end
