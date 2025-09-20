





mutable struct cstr
	
	# Inlets are stored in switches.connectionInstance

	nComp::Int64
	V::Float64
	cIn::Vector{Float64}
	
	c0::Union{Float64, Vector{Float64}, Int64} # defaults to 0
	unitStride::Int64
	
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	
	
	"""
		cstr(; nComp, V, c0=0)

	Constructs a continuously stirred tank reactor (CSTR) unit for use in CADET-Julia simulations.

	# Arguments
	- `nComp`: Number of components in the reactor.
	- `V`: Volume of the reactor.
	- `c0`: (Optional) Initial concentration vector for each component. Defaults to zeros if not provided.

	# Returns
	A `cstr` struct containing all parameters and initial state for a constant-volume CSTR.

	# Details
	- Supports inlets from both external sources and other units (columns or CSTRs).
	- Assumes constant reactor volume.
	- The governing equation is: 
		dc/dt = Q_tot/V (Q_col * c_col + Q_inlet * c_inlet - c) 
	"""
	function cstr(; nComp, V, c0 = 0)
		
		# Set initial condition vector 
		if c0 == 0
			c0 = zeros(Float64, nComp)
		elseif length(c0) == nComp #if initial conditions for each component is given
			nothing
		else 
			throw(error("Initial concentrations incorrectly written"))
		end
		cIn = [0.0]
		unitStride = nComp
		
		# Solution_outlet as well 
		solution_outlet = zeros(Float64,1,nComp)
		solution_times = Float64[]


		new(nComp, V, cIn, c0, unitStride, solution_outlet, solution_times)
	end
end



"""
    compute!(RHS, RHS_q, cpp, qq, x, m::cstr, t, section, sink, switches, idx_units)

Computes the right-hand side (RHS) of the ODE system for a continuously stirred tank reactor (CSTR) unit.

# Arguments
- `RHS`: Vector to store the computed derivatives (right-hand side) for concentrations.
- `RHS_q`, `cpp`, `qq`: Additional arguments for compatibility (not used in CSTR).
- `x`: State vector containing current concentrations.
- `m::cstr`: The CSTR unit instance.
- `t`: Current simulation time.
- `section`: Current section index (for switching conditions).
- `sink`: Index of the current unit (CSTR).
- `switches`: Switches object containing flow and inlet condition information.
- `idx_units`: Vector of starting indices for each unit in the global state vector.

# Details
- Updates the inlet flow and concentration for the CSTR based on current time and switching conditions.
- Computes the time derivative of the concentration for each component according to the CSTR mass balance
"""

function compute!(RHS, RHS_q, cpp, qq, x, m::cstr, t, section, sink, switches, idx_units) 
	
	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)
    
	# Loop over components
	@inbounds for j = 1:m.nComp

		# Determining inlet concentration 
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j]) 

		RHS[j + idx_units[sink]] = switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink] / m.V * (m.cIn[1] - x[j + idx_units[sink]])
	end
	
    nothing
end





