mutable struct RadialConvDispOp
	# A struct containing the DG variables for the convection dispersion operator which is the same for LRM, LRMP and GRM. 
	# DG properties for the convection dispersion operator
	polyDeg::Int64
	nCells::Int64
	nNodes::Int64
    nPoints::Int64
    strideNode::Int64
    strideCell::Int64

	nodes::Vector{Float64}
    invWeights::Vector{Float64}
    invMM::Matrix{Float64}
	MM00::Matrix{Float64}
	MM01::Matrix{Float64}
	polyDerM::Matrix{Float64}
    deltarho::Float64
    rho_i::Vector{Float64}
    rho_ip1::Vector{Float64}

	# Allocation vectors and matrices
	mul1::Vector{Float64}
    c_star::Vector{Float64}
    g_star::Vector{Float64}
	Dc::Vector{Float64}
    Dg::Vector{Float64}
    h::Vector{Float64}

	function RadialConvDispOp(polyDeg, nCells, col_inner_radius, col_outer_radius)

		nNodes = polyDeg + 1
		nPoints = nNodes * nCells	#Number of cells times number of nodes per component
		strideNode = 1				#Number of points to next concentration in the state vector, here 1
		strideCell = nNodes * strideNode

		# delta rho [m]
		deltarho = (col_outer_radius - col_inner_radius) / nCells
		rho_i  = [col_inner_radius + (cell - 1) * deltarho for cell in 1:nCells] # face radii
        rho_ip1 = [col_inner_radius + cell * deltarho for cell in 1:nCells]

		# Obtain LGL nodes and weights
		nodes, invWeights = DGElements.lglnodes(polyDeg) #LGL nodes & weights
		invMM = DGElements.invMMatrix(nodes, polyDeg) #Inverse mass matrix
		MM00 = DGElements.MMatrix(nodes, polyDeg)
		MM01 = DGElements.MMatrix(nodes, polyDeg, 0, 1)
		polyDerM = DGElements.derivativeMatrix(polyDeg, nodes) #derivative matrix

		# allocation vectors and matrices
		mul1 = zeros(Float64, nNodes)
		c_star = zeros(Float64, nCells + 1)
		g_star = zeros(Float64, nCells + 1)
		Dc = zeros(Float64, nPoints)
		Dg = zeros(Float64, nPoints)
		h = zeros(Float64, nPoints)

		new(polyDeg, nCells, nNodes, nPoints, strideNode, strideCell, nodes, invWeights, invMM, MM00, MM01, polyDerM, deltarho, rho_i, rho_ip1, mul1, c_star, g_star, Dc, Dg, h)
	end
end

################################# TRANSPORT MODELS #################################
abstract type RadialModelBase 
	# From here, the transport models are found
end

################################# LUMPED RATE MODEL (rLRM) #################################
	mutable struct rLRM <: RadialModelBase
	# Check parameters
	# These parameters are the minimum to be specified for the LRM
	nComp::Int64 
    col_inner_radius::Float64
    col_outer_radius::Float64
	#col_height::Float64
	cross_section_area::Float64
    d_rad::Union{Float64, Vector{Float64}}
    eps_c::Float64
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
    
	cIn::Vector{Float64}

	polyDeg::Int64
	nCells::Int64	
	
    # Convection Dispersion properties
	RadialConvDispOpInstance::RadialConvDispOp

	# based on the input, remaining properties are calculated in the function LRM
	#Determined properties 
	bindStride::Int64
	adsStride::Int64
	unitStride::Int64

	# Allocation vectors and matrices
    idx::UnitRange{Int64}
	Fc::Float64
	Fjac::Float64
	cpp::Vector{Float64}
    RHS_q::Vector{Float64}
    qq::Vector{Float64}
	RHS::Vector{Float64}
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	bind::bindingBase
	
	

	# Default variables go in the arguments in the LRM
	function rLRM(; nComp, col_inner_radius, col_outer_radius, d_rad, eps_c, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, cross_section_area=1.0)
		
		# Get necessary variables for convection dispersion DG 
		RadialConvDispOpInstance = RadialConvDispOp(polyDeg, nCells, col_inner_radius, col_outer_radius)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = RadialConvDispOpInstance.nPoints
		adsStride = 0  #stride to guide to liquid adsorption concentrations i.e., cpp
		unitStride = adsStride + bindStride*nComp*2

		# allocation vectors and matrices
		idx = 1:RadialConvDispOpInstance.nPoints
		Fc = 0.0
		Fjac = 0.0
		cpp = zeros(Float64, RadialConvDispOpInstance.nPoints * nComp)
		RHS_q = zeros(Float64, RadialConvDispOpInstance.nPoints * nComp)
		qq = zeros(Float64, RadialConvDispOpInstance.nPoints * nComp)
		RHS = zeros(Float64, adsStride + 2*nComp*bindStride)
		cIn = zeros(Float64, nComp)
		
		# Set initial condition vectors 
		c0, cp0, q0 = initial_condition_specification(nComp, RadialConvDispOpInstance, bindStride, c0, cp0, q0)
		
		# Solution_outlet as well 
		solution_outlet = zeros(Float64,1,nComp)
		solution_times = Float64[]

		# Default binding - assumes linear with zero binding 
		bind = Linear(
					ka = zeros(Float64,nComp),
					kd = zeros(Float64,nComp),
					is_kinetic = true, #if false, a high kkin is set to approximate rapid eq. if true, kkin=1
					nBound = zeros(Bool,nComp), # Number of bound components, specify non-bound states by a zero, defaults to assume all bound states e.g., [1,0,1]
					bindStride = bindStride, # Not necessary for Linear model, only for Langmuir and SMA
					# nBound =  [1,0,1,1] # Specify non-bound states by a zero, defaults to assume all bound states
					)
		
		# The new commando must match the order of the elements in the struct!
		new(nComp, col_inner_radius, col_outer_radius, cross_section_area, d_rad, eps_c, c0, cp0, q0, cIn, polyDeg, nCells, RadialConvDispOpInstance, bindStride, adsStride, unitStride, idx, Fc, Fjac, cpp, RHS_q, qq, RHS, solution_outlet, solution_times, bind)
	end
end


# Define a function to compute the transport term for the rLRM
function compute_transport!(RHS, RHS_q, cpp, x, m::rLRM, t, section, sink, switches, idx_units) 
	# section = i from call 
	# sink is the unit i.e., h from previous call
	
	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)
	
	nCells = m.RadialConvDispOpInstance.nCells
	rho_i   = m.RadialConvDispOpInstance.rho_i
	rho_ip1 = m.RadialConvDispOpInstance.rho_ip1

	@inbounds for j = 1:m.nComp
		# inlet conc for comp j
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j])

		# inlet velocity for this section/sink
		v_in = switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink]

		# face velocities (constant volumetric flow ⇒ v ∝ 1/r)
		radial_v = Vector{Float64}(undef, nCells + 1)
		radial_v[1] = v_in
		for f in 2:nCells
			radial_v[f] = v_in * (rho_i[1] / rho_i[f])
		end
		radial_v[nCells + 1] = v_in * (rho_i[1] / rho_ip1[nCells])

		# diffusive face scales for comp j
		d = (m.d_rad isa Vector) ? m.d_rad[j] : m.d_rad
		left_scale_vec  = rho_i   .* d
		right_scale_vec = rho_ip1 .* d

		# slice for comp j
		m.idx = 1 + (j - 1) * m.RadialConvDispOpInstance.nPoints : j * m.RadialConvDispOpInstance.nPoints

		cpp_block = @view x[1 + idx_units[sink] : idx_units[sink] + m.RadialConvDispOpInstance.nPoints * m.nComp]

		RadialConvDispOperatorDG.radialresidualImpl!(m.RadialConvDispOpInstance.Dc, cpp_block, m.idx, m.RadialConvDispOpInstance.strideNode, m.RadialConvDispOpInstance.strideCell, m.RadialConvDispOpInstance.nNodes, m.RadialConvDispOpInstance.nCells, m.RadialConvDispOpInstance.deltarho, m.polyDeg, m.RadialConvDispOpInstance.invWeights, m.RadialConvDispOpInstance.polyDerM, m.RadialConvDispOpInstance.invMM, m.RadialConvDispOpInstance.MM00, m.RadialConvDispOpInstance.MM01, m.RadialConvDispOpInstance.nodes, v_in, d, m.cIn[j], m.RadialConvDispOpInstance.c_star, m.RadialConvDispOpInstance.g_star, m.RadialConvDispOpInstance.Dg, m.RadialConvDispOpInstance.h, m.RadialConvDispOpInstance.mul1, rho_i, rho_ip1, radial_v, left_scale_vec, right_scale_vec)

		@views RHS[m.idx .+ idx_units[sink]] .= m.RadialConvDispOpInstance.Dc
	end
	
    nothing
end