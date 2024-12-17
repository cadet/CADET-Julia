





mutable struct cstr
	
	# Inlets are stored in switches.connectionInstance

	nComp::Int64
	V::Float64
	cIn::Vector{Float64}
	
	c0::Union{Float64, Vector{Float64}, Int64} # defaults to 0
	unitStride::Int64
	
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	
	

	function cstr(; nComp, V, c0 = 0)
		"""
		A struct containing the parameters for the continuously stirred tank reactor 
		Supports inlets from inlet and other units such as columns or cstr 
		Support constant volumne cstr only 
		The equation for that is 
		dc/dt = Q_tot/V (Q_col * c_col + Q_inlet * c_inlet - c) 
		"""
		
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



# Define a function to compute the transport term for the cstr
function compute!(RHS, RHS_q, cpp, qq, x, m::cstr, t, section, sink, switches, idx_units) 
	# section = i from call 
	# sink is the unit i.e., h from previous call
	
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





