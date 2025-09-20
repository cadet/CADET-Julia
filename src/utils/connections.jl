

mutable struct CreateInlet 
	# A struct for creating an inlet 
	# An inlet here is specified as inlet at each section time
	# All the inlet coefficients are zero by default 
	cIn_c::Matrix{Float64} # Concentrations for each switch, unit and component
    cIn_l::Matrix{Float64}
	cIn_q::Matrix{Float64}
	cIn_cube::Matrix{Float64}
	
	"""
    	CreateInlet(; nComp::Int64, nSections::Int64)

	Constructs an inlet unit for the simulation, storing inlet concentration coefficients for each section and component.

	# Arguments
	- `nComp`: Number of components.
	- `nSections`: Number of section times.

	# Fields
	- `cIn_c`, `cIn_l`, `cIn_q`, `cIn_cube`: Matrices of inlet concentration coefficients (constant, linear, quadratic, cubic) for each section and component.
	"""
	function CreateInlet(; nComp::Int64, nSections::Int64)

		# Instantiate the inlets 
		cIn_c = -ones(Float64, nSections, nComp)
		cIn_l = -ones(Float64, nSections, nComp)
		cIn_q = -ones(Float64, nSections, nComp)
		cIn_cube = -ones(Float64, nSections, nComp)

		# Set default values 
		cIn_c[1,1:end] = zeros(Float64, nComp)
		cIn_l[1,1:end] = zeros(Float64, nComp)
		cIn_q[1,1:end] = zeros(Float64, nComp)
		cIn_cube[1,1:end] = zeros(Float64, nComp)

        new(cIn_c, cIn_l, cIn_q, cIn_cube)
	end 

end

"""
    modify_inlet!(; inlet, nComp, section, cIn_c, cIn_l, cIn_q, cIn_cube)

Modifies the inlet concentration coefficients for a given section of a `CreateInlet` object.

# Arguments
- `inlet`: The `CreateInlet` instance to modify.
- `nComp`: Number of components.
- `section`: Section index to modify.
- `cIn_c`, `cIn_l`, `cIn_q`, `cIn_cube`: Vectors of new coefficients for the section (optional, default zeros).

# Details
- Updates the specified section with new coefficients.
- Repeats the pattern for remaining sections if needed.
""" 
function modify_inlet!(; inlet::CreateInlet, nComp::Int64, section, cIn_c = [0.0], cIn_l=[0.0], cIn_q=[0.0], cIn_cube=[0.0])
	if size(cIn_c)[1] == nComp
		inlet.cIn_c[section,:] = cIn_c 
	elseif cIn_c == [0.0]
			inlet.cIn_c[section,:] = zeros(Float64, nComp)
	else 
		throw("Incorrect specification") 
	end 
	
	if size(cIn_l)[1] == nComp
		inlet.cIn_l[section,:] = cIn_l
	elseif cIn_l == [0.0]
		inlet.cIn_l[section,:] = zeros(Float64, nComp)
	else
		throw("Incorrect specification")
	end

	if size(cIn_q)[1] == nComp
		inlet.cIn_q[section,:] = cIn_q
	elseif cIn_q == [0.0]
		inlet.cIn_q[section,:] = zeros(Float64, nComp)
	else
		throw("Incorrect specification")
	end


	if size(cIn_q)[1] == nComp
		inlet.cIn_q[section,:] = cIn_q
	elseif cIn_cube == [0.0]
		inlet.cIn_cube[section,:] = zeros(Float64, nComp)
	else
		throw("Incorrect specification")
	end

	# Repeat pattern that has been specified for the remaining inlets 
	# The reason for this is that one does not want to specicify everytime the inlet changes 
	# also, if the inlet is constant or follows a pattern, this should be specified 
	inlet.cIn_c = repeat_elements(inlet.cIn_c, section)
	inlet.cIn_l = repeat_elements(inlet.cIn_l, section)
	inlet.cIn_q = repeat_elements(inlet.cIn_q, section)
	inlet.cIn_cube = repeat_elements(inlet.cIn_cube, section)
end

"""
    repeat_elements(matrix::Matrix, idx::Int)

Repeats the pattern of the first `idx` rows of a matrix to fill the entire matrix.

# Arguments
- `matrix`: Input matrix.
- `idx`: Number of rows to use as the repeating pattern.

# Returns
A matrix with the pattern repeated to match the original size.
"""
function repeat_elements(matrix::Matrix, idx::Int)
    nrows, ncols = size(matrix)
    pattern = view(matrix, 1:idx, :)
    repetitions = div(nrows, idx)
    remainder = mod(nrows, idx)
    
    repeated_pattern = repeat(pattern, repetitions, 1)
    remainder_pattern = view(pattern, 1:remainder, :)
    
    return vcat(repeated_pattern, remainder_pattern)
end


mutable struct CreateOutlet 
	# A struct for creating an outlet 
	# An outlet is specified as a container that holds specific outputs
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	idx_unit::Vector{Vector{Int64}}
	idx_outlet::Vector{Vector{Int64}}
	u_outlet::Vector{Vector{Float64}} # velocities from connections
	u_tot::Vector{Float64} # total velocities from connections

	"""
    	CreateOutlet(; nComp::Int64)

	Constructs an outlet unit for the simulation, storing outlet solution data and connection indices.

	# Arguments
	- `nComp`: Number of components.

	# Fields
	- `solution_outlet`: Matrix to store outlet solutions.
	- `solution_times`: Vector of time points for outlet solutions.
	- `idx_unit`, `idx_outlet`: Indices of connected units and outlets.
	- `u_outlet`, `u_tot`: Velocities for each switch and total velocities.
	"""
	function CreateOutlet(;nComp::Int64)
		solution_outlet = zeros(Float64,0,nComp)
		solution_times = Float64[]
		idx_unit = [[-1]]
		idx_outlet = [[-1]]
		u_outlet = [[-1.]]
		u_tot = zeros(Float64,1)
		new(solution_outlet, solution_times, idx_unit, idx_outlet, u_outlet, u_tot)
	end 
end 


"""
    Connection(...)

Defines the connections between units (inlets, columns, CSTRs, outlets) in the simulation, including flow rates, concentration coefficients, and connection indices.

# Constructors
- For inlets, columns, CSTRs, and outlets as sources/sinks.
- For default configuration (initializes empty connection matrices).

# Fields
- `u_unit`, `u_inlet`, `u_tot`: Velocities for each switch and unit.
- `c_connect`, `idx_connect`: Connection matrices and indices.
- `cIn_c`, `cIn_l`, `cIn_q`, `cIn_cube`: Inlet concentration coefficients.
- `Q_inlet_*`, `Q_unit_*`: Flow rate coefficients.
- `dynamic_flow`: Dynamic flow specification for each switch and unit.
"""
mutable struct Connection	

	u_unit::Vector{Vector{Vector{Float64}}} # Inlet comming from another unit i.e. switch, unit(sink), value
	u_inlet::Vector{Vector{Vector{Float64}}} # Inlet comming an inlet i.e. switch, unit(sink), value
	u_tot::Matrix{Float64} # sum of u_unit and u_inlet, defaults to -1 to avoid dividing with zero and to check whether a switch is specified 
	c_connect::Vector{Vector{Vector{Vector{Int64}}}} # Connections for each switch, unit(sink), component, value (1 or 0)
	idx_connect::Vector{Vector{Vector{Vector{Int64}}}} # Connection indices for each switch, unit(sink), component, value (index)
	# For the concentration specifications, they should replace the latest used for the switch time. 
	# So if inlet 1 is specified for section 1, it is the same for section 2 unless anything else is specified 
	# for the switches, they should follow a repeated pattern like for an SMB such that you only have to specify a cycle once
	cIn_c::Vector{Vector{Vector{Vector{Float64}}}} # Concentrations for each switch, unit(sink), component, value
    cIn_l::Vector{Vector{Vector{Vector{Float64}}}}
	cIn_q::Vector{Vector{Vector{Vector{Float64}}}}
	cIn_cube::Vector{Vector{Vector{Vector{Float64}}}}

	Q_inlet_c::Vector{Vector{Vector{Float64}}} # Inlet flow rates for each switch, unit(sink), value
    Q_inlet_l::Vector{Vector{Vector{Float64}}}
	Q_inlet_q::Vector{Vector{Vector{Float64}}}
	Q_inlet_cube::Vector{Vector{Vector{Float64}}}
	Q_unit_c::Vector{Vector{Vector{Float64}}} 
    Q_unit_l::Vector{Vector{Vector{Float64}}}
	Q_unit_q::Vector{Vector{Vector{Float64}}}
	Q_unit_cube::Vector{Vector{Vector{Float64}}}
	dynamic_flow 


	# For these cosntructors, we cannot use keyword arguments when having multiple constructors 
	
	# If an inlet is specified as input
	function Connection(switches, switch::Int64, section::Int64, source::CreateInlet, sink::Int64, u::Float64, Q_c, Q_l, Q_q, Q_cube, dynamic_flow_check)
		# Here input is the inlet concentration 

		# The switchSetup is set to a section 
		switches.switchSetup[section] = switch
		
		# The inlet flux is and u_tot are modified 
		push!(switches.ConnectionInstance.u_inlet[switch][sink], u)
		switches.ConnectionInstance.u_tot[switch, sink] += u

		# The section specifications from the inlet should be copied 
		# if there are multiple sections but only one switch, fill out all the sections with the same values
		if switches.nSwitches<2
			# Set inlet concentration
			for i in 1:length(switches.ConnectionInstance.cIn_c) # For each section
				for j in 1:length(switches.ConnectionInstance.cIn_c[i][sink]) # for each component
						push!(switches.ConnectionInstance.cIn_c[i][sink][j], source.cIn_c[i, j])
						push!(switches.ConnectionInstance.cIn_l[i][sink][j], source.cIn_l[i, j])
						push!(switches.ConnectionInstance.cIn_q[i][sink][j], source.cIn_q[i, j])
						push!(switches.ConnectionInstance.cIn_cube[i][sink][j], source.cIn_cube[i, j])
				end
			end

		else
			# In between two section times and 1 switch, it will assume previous switch and repeat that. 
			# This is how it is in CADET-Core - see switchTest.py that tests this. 
			# For example, if switch is constant throughout two section times, this should only be specified once. 
			# In this case, what happens if you specify 2 switches in a 5 section. Switch 2 begins at section 2. 
			# Hence [switch1, ?, switch2, ?, ?]
			# The result is that it assumes swicht1 for whole period i.e., [switch1, switch1, switch1, switch1, switch1] 
			# If [switch1, switch2, ?, ?, ?] is specified, 
			# the result is [switch1, switch2, switch1, switch2, switch1]

			# Set inlet concentration 
			for j in 1:length(switches.ConnectionInstance.cIn_c[section][sink]) # for each component
				push!(switches.ConnectionInstance.cIn_c[section][sink][j], source.cIn_c[section, j])
				push!(switches.ConnectionInstance.cIn_l[section][sink][j], source.cIn_l[section, j])
				push!(switches.ConnectionInstance.cIn_q[section][sink][j], source.cIn_q[section, j])
				push!(switches.ConnectionInstance.cIn_cube[section][sink][j], source.cIn_cube[section, j])
			end
			
		end

		# If dynamic flow has been specified 
		if dynamic_flow_check == true
			switches.ConnectionInstance.dynamic_flow[switch, sink] = YesDynamicFlow()
			push!(switches.ConnectionInstance.Q_inlet_c[switch][sink], Q_c)
			push!(switches.ConnectionInstance.Q_inlet_l[switch][sink], Q_l)
			push!(switches.ConnectionInstance.Q_inlet_q[switch][sink], Q_q)
			push!(switches.ConnectionInstance.Q_inlet_cube[switch][sink], Q_cube)
		end
	end
	
	# If a column is specified as input i.e., in series 
	function Connection(switches, switch::Int64, section::Int64, source::Tuple{Int64,ModelBase}, sink::Int64, u::Float64, Q_c, Q_l, Q_q, Q_cube, dynamic_flow_check)
		# Unpack source 
		columnNumber = source[1]
		model = source[2]

		# The switchSetup is set to a section 
		switches.switchSetup[section] = switch
		
		# Set inlet velocities 
		push!(switches.ConnectionInstance.u_unit[switch][sink], u)
		switches.ConnectionInstance.u_tot[switch, sink] += u

		# Connection matrix and indices should be edited
		for j =1:model.nComp
			push!(switches.ConnectionInstance.c_connect[switch][sink][j], 1.0) # connection matrix, should be edited if some components are not connected

			# Indexing of the connection matrix
			push!(switches.ConnectionInstance.idx_connect[switch][sink][j], model.ConvDispOpInstance.nPoints + (j-1) * model.ConvDispOpInstance.nPoints + switches.idx_units[columnNumber])
		end

		# If dynamic flow has been specified 
		if dynamic_flow_check == true
			switches.ConnectionInstance.dynamic_flow[switch, sink] = YesDynamicFlow()
			push!(switches.ConnectionInstance.Q_unit_c[switch][sink], Q_c)
			push!(switches.ConnectionInstance.Q_unit_l[switch][sink], Q_l)
			push!(switches.ConnectionInstance.Q_unit_q[switch][sink], Q_q)
			push!(switches.ConnectionInstance.Q_unit_cube[switch][sink], Q_cube)
		end
		
	end

	# If a cstr is specified as input i.e., in series 
	function Connection(switches, switch::Int64, section::Int64, source::Tuple{Int64,cstr}, sink::Int64, u::Float64, Q_c, Q_l, Q_q, Q_cube, dynamic_flow_check)
		# Unpack source 
		columnNumber = source[1]
		model = source[2]

		# The switchSetup is set to a section 
		switches.switchSetup[section] = switch
		
		# Here the matrix should be edited
		# edit u_total, Q and C_connect 
		push!(switches.ConnectionInstance.u_unit[switch][sink], u)
		switches.ConnectionInstance.u_tot[switch, sink] += u

		# Connection matrix and indices should be edited
		for j =1:model.nComp
			push!(switches.ConnectionInstance.c_connect[switch][sink][j], 1.0) # connection matrix, should be edited if some components are not connected

			# Indexing of the connection matrix 
			push!(switches.ConnectionInstance.idx_connect[switch][sink][j], j + switches.idx_units[columnNumber])
		end
		
		# If dynamic flow has been specified 
		if dynamic_flow_check == true
			switches.ConnectionInstance.dynamic_flow[switch, sink] = YesDynamicFlow()
			switches.ConnectionInstance.Q_unit_c[switch, sink] = Q_c
			switches.ConnectionInstance.Q_unit_l[switch, sink] = Q_l
			switches.ConnectionInstance.Q_unit_q[switch, sink] = Q_q
			switches.ConnectionInstance.Q_unit_cube[switch, sink] = Q_cube
		end
		
		# If dynamic flow has been specified 
		if dynamic_flow_check == true
			switches.ConnectionInstance.dynamic_flow[switch, sink] = YesDynamicFlow()
			switches.ConnectionInstance.Q_unit_c[switch, sink] = Q_c
			switches.ConnectionInstance.Q_unit_l[switch, sink] = Q_l
			switches.ConnectionInstance.Q_unit_q[switch, sink] = Q_q
			switches.ConnectionInstance.Q_unit_cube[switch, sink] = Q_cube
		end
		
	end

	# Default configuration 
	function Connection(nSections::Int64, nSwitches::Int64, nComp::Int64, nColumns::Int64)
		
		# Here the matrix should be edited
		u_unit = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
        u_inlet = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		u_tot = zeros(Float64, nSwitches, nColumns) 
		c_connect =  [ [ [ Float64[] for _ in 1:nComp ] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		idx_connect = [ [ [ Int64[] for _ in 1:nComp ] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		cIn_c = [ [ [ Float64[] for _ in 1:nComp ] for _ in 1:nColumns ] for _ in 1:nSections ]
		cIn_l = [ [ [ Float64[] for _ in 1:nComp ] for _ in 1:nColumns ] for _ in 1:nSections ]
		cIn_q = [ [ [ Float64[] for _ in 1:nComp ] for _ in 1:nColumns ] for _ in 1:nSections ]
		cIn_cube = [ [ [ Float64[] for _ in 1:nComp ] for _ in 1:nColumns ] for _ in 1:nSections ]
		
		# For dynamic flow specifications 
		Q_inlet_c = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ] # Inlet flow rates for each switch, unit, defaults to zeros 
		Q_inlet_l = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		Q_inlet_q = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		Q_inlet_cube = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		Q_unit_c = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		Q_unit_l = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		Q_unit_q = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		Q_unit_cube = [ [ Float64[] for _ in 1:nColumns ] for _ in 1:nSwitches ]
		dynamic_flow = Array{Any}(undef, nSwitches, nColumns)
		for i in 1:nSwitches, j in 1:nColumns
			dynamic_flow[i,j] = NoDynamicFlow()
		end
		
		new(u_unit, u_inlet, u_tot, c_connect, idx_connect, cIn_c, cIn_l, cIn_q, cIn_cube, Q_inlet_c, Q_inlet_l, Q_inlet_q, Q_inlet_cube, Q_unit_c, Q_unit_l, Q_unit_q, Q_unit_cube, dynamic_flow)
	end
	
	
	# If an outlet is specified as sink
	function Connection(switches, switch::Int64, section::Int64, source::Int64, sink::CreateOutlet, u::Float64, Q_in_c, Q_in_l, Q_in_q, Q_in_cube, dynamic_flow_check)
		# The switchSetup is set to a section 
		switches.switchSetup[section] = switch

		# if first time using this function, create outlet indicies
		if sink.idx_unit == [[-1]]

			# idx_unit and idx_outlet will have connections for each switch
			# Multiple units can connect to each outlet
			sink.idx_unit = [ Int[] for _ in 1:switches.nSwitches ]
			sink.idx_outlet = [ Int[] for _ in 1:switches.nSwitches ]

			sink.u_outlet = [ Int[] for _ in 1:switches.nSwitches ]
			sink.u_tot = zeros(Float64, switches.nSwitches)
		end

		# Modify the inlet velocities
		push!(sink.u_outlet[switch], u)
		sink.u_tot[switch] += u

		# It should determine the right indicies from which the output from the source should be written to. 
		push!(sink.idx_unit[switch], source)
	end
	
end

mutable struct Switches
	# The general purpose of the switches is to determine the inlet concentrations
	# The inlet formulation is: cin[switch, col, comp] = (u_inlet * (c_c + c_l*t +...+) + u_unit c_unit) / (u_inlet + u_unit)
	# c_unit = c_connect * x[idx_connect] 
	# Hence c_connect and Finlets are zero by default

	nSections::Int64
	section_times::Vector{Float64}
	nSwitches::Int64
	switchSetup::Vector{Int64}
	

	ConnectionInstance::Connection 
	idx_units::Vector{Int64}
	inlet_conditions::Array{InletConditions, 3}

	"""
		Switches(; nSections, section_times, nSwitches, nColumns, nComp, idx_units)

	Stores and manages the switching logic for inlet concentrations and connections in the simulation.

	# Arguments
	- `nSections`: Number of section times.
	- `section_times`: Vector of section time points.
	- `nSwitches`: Number of switches.
	- `nColumns`: Number of columns/units.
	- `nComp`: Number of components.
	- `idx_units`: Indices for each unit in the global state vector.

	# Fields
	- `ConnectionInstance`: The `Connection` object managing all connections.
	- `switchSetup`: Vector specifying the switch used at each section.
	- `inlet_conditions`: Array of inlet condition types for each section, unit, and component.
	"""	
	function Switches(; nSections::Int64, section_times::Vector{Float64}, nSwitches::Int64, nColumns::Int64, nComp::Int64, idx_units::Vector{Int64})
	
		# Establish default zero configuration
		ConnectionInstance = Connection(nSections, nSwitches, nComp, nColumns) # connection(nSections = nSections, nComp = nComp, nColumns = nColumns)
		switchSetup = -ones(Int64,nSections)
		inlet_conditions = fill(StaticInlets(), nSections, nColumns, nComp)
	
		new(nSections, section_times, nSwitches, switchSetup, ConnectionInstance, idx_units, inlet_conditions)
	end
end
