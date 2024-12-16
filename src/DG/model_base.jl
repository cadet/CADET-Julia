





mutable struct ConvDispOp
	# A struct containing the DG variables for the convection dispersion operator which is the same for LRM, LRMP and GRM. 
	# The variables are determined based on polyDeg (in mobile phase), nCells and colLength.  

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
    polyDerM::Matrix{Float64}
    deltaZ::Float64

	# Allocation vectors and matrices
	mul1::Vector{Float64}
    c_star::Vector{Float64}
    h_star::Vector{Float64}
	Dc::Vector{Float64}
    Dh::Vector{Float64}
    h::Vector{Float64}

	function ConvDispOp(polyDeg,nCells,colLength)

		nNodes = polyDeg + 1
		nPoints = nNodes*nCells #Number of axial cells times number of nodes per component
		strideNode = 1  # number of points to next concentration in the state vector, here 1
		strideCell = nNodes * strideNode

		# Obtain LGL nodes and weights
		nodes,invWeights = DGElements.lglnodes(polyDeg) #LGL nodes & weights
		invMM = DGElements.invMMatrix(nodes,polyDeg) #Inverse mass matrix
		polyDerM = DGElements.derivativeMatrix(polyDeg,nodes) #derivative matrix
		
		# delta z [m]
		deltaZ = colLength / nCells
		
		# allocation vectors and matrices
		mul1 = zeros(Float64, length(polyDerM*nodes))
		c_star = zeros(Float64,nCells+1)
		h_star = zeros(Float64,nCells+1)
		Dc = zeros(Float64, nPoints)
		Dh = zeros(Float64, nPoints)
		h = zeros(Float64, nPoints)

		new(polyDeg, nCells, nNodes, nPoints, strideNode, strideCell, nodes, invWeights, invMM, polyDerM, deltaZ, mul1, c_star, h_star, Dc, Dh, h)
	end
end


# Pore phase 
mutable struct PoreOp
	# A struct containing the DG variables for the convection dispersion operator which is the same for LRM, LRMP and GRM. 
	# The variables are determined based on polyDeg (in mobile phase), nCells and colLength.  

	# DG properties for the convection dispersion operator
	stridePore::Int64
	LiftMatrixRed::Vector{Float64}
	invMM_Ar::Matrix{Float64}
	nNodesPore::Int64
	boundaryPore::Vector{Float64}
	term1::Vector{Float64}

	function PoreOp(polyDegPore,nPoints,nComp,Rc, Rp, deltaR, eps_p)

		nNodesPore = polyDegPore + 1
		stridePore = nPoints*nNodesPore*nComp #stride pore phase

		# Obtain LGL nodes, weights for pore phase
		nodesPore,invWeightsPore = DGElements.lglnodes(polyDegPore)              #LGL nodes & weights

		## Matrices for term 1 - volume integral
		MMatrix = DGElements.MMatrix(nodesPore, polyDegPore, 0, 2) * (deltaR/2)^2
		MMatrix += DGElements.MMatrix(nodesPore, polyDegPore, 0, 1) * (deltaR) * Rc
		MMatrix += DGElements.MMatrix(nodesPore, polyDegPore, 0, 0) * Rc^2
		
		# M^-1 
		invMM = inv(MMatrix)
		secondOrdersteifMatrix = Rc^2 * DGElements.second_order_stiff_matrix(nodesPore,polyDegPore, 0, 0)
		secondOrdersteifMatrix += Rc * deltaR * DGElements.second_order_stiff_matrix(nodesPore,polyDegPore, 0, 1)
		secondOrdersteifMatrix += (deltaR/2)^2 * DGElements.second_order_stiff_matrix(nodesPore,polyDegPore, 0, 2)

		# M^-1 * A^r * (2/deltaR)^2
		invMM_Ar = invMM * secondOrdersteifMatrix .* (2/deltaR)^2

		# # Ir vector 
		# Ir = zeros(Float64, nNodesPore) #_disc.Ir
		# for i=1:nNodesPore
		# 	Ir[i] = (deltaR/2) * (nodesPore[i] + 1) + Rc #
		# end
		# Ir = Ir.^2 # For spheres

		# Lifting matrix, L = M^-1 eps
		LiftMatrix = zeros(nNodesPore,2)
		LiftMatrix[1,1] = 1
		LiftMatrix[nNodesPore,2] = 1
		LiftMatrix = invMM*LiftMatrix
		LiftMatrixRed = LiftMatrix[:,2]  #Since only the second column is used according to the boundary conditions
		LiftMatrixRed *= (2/deltaR) * Rp^2 / eps_p # L (M^-1) (2/deltaR) Rp^2/eps_p


		boundaryPore = zeros(Float64,nPoints)

		term1 = zeros(nNodesPore)

		new(stridePore, LiftMatrixRed, invMM_Ar, nNodesPore, boundaryPore, term1)
	end
end

# A function that determines the inlet concentrations needed for the transport equations 
# Depends on the inlets specified in the switches
function inlet_concentrations!(cIn, switches, j, section, sink, x, t)
	"""
		A function that determines the inlet concentrations needed for the transport equations. 
		Depends on the inlets specified in the switches. 
		input is: 
		cIn: the inlet concentration
		switches: the switches struct
		j: the component index
		section: the section index
		sink: the sink (unit)
		x: the state vector
		t: the time

		output is a modification of:
		cIn: the inlet concentration
		
	"""
	# Restarting the inlet concentration
	cIn[1] = 0.0
    
	# Inlets from inlets
	for l in eachindex(switches.ConnectionInstance.cIn_c[section][sink][j])
		cIn[1] += (switches.ConnectionInstance.cIn_c[section][sink][j][l] + 
			switches.ConnectionInstance.cIn_l[section][sink][j][l]*t +
			switches.ConnectionInstance.cIn_q[section][sink][j][l]*t^2 +
			switches.ConnectionInstance.cIn_cube[section][sink][j][l]*t^3) * switches.ConnectionInstance.u_inlet[switches.switchSetup[section]][sink][l]
	end

	# Inlets from other units 
	for l in eachindex(switches.ConnectionInstance.c_connect[switches.switchSetup[section]][sink][j])
		# println("section = $section, sink = $sink, j = $j, l = $l")
		cIn[1] += switches.ConnectionInstance.u_unit[switches.switchSetup[section]][sink][l] * switches.ConnectionInstance.c_connect[switches.switchSetup[section]][sink][j][l] * x[switches.ConnectionInstance.idx_connect[switches.switchSetup[section]][sink][j][l]]
	end 

	# Divide by total velocity according to mass balance 
	cIn[1] /= switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink]
	nothing
end



################################# TRANSPORT MODELS #################################
abstract type ModelBase 
	# From here, the transport models are found
end

################################# LUMPED RATE MODEL (LRM) #################################
mutable struct LRM <: ModelBase
	# Check parameters
	# These parameters are the minimum to be specified for the LRM
	nComp::Int64 
    colLength::Float64
	cross_section_area::Float64
    d_ax::Union{Float64, Vector{Float64}}
    eps_c::Float64
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
	
    
	cIn::Vector{Float64}
	exactInt::Int64

	polyDeg::Int64
	nCells::Int64	
	
    # Convection Dispersion properties
	ConvDispOpInstance::ConvDispOp

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
    RHS_q::Vector{Real}
    qq::Vector{Float64}
	RHS::Vector{Real}
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	bind::bindingBase
	
	

	# Default variables go in the arguments in the LRM(..)
	function LRM(; nComp, colLength, d_ax, eps_c, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, exactInt=1, cross_section_area=1.0)
		
		# Get necessary variables for convection dispersion DG 
		ConvDispOpInstance = ConvDispOp(polyDeg,nCells,colLength)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = ConvDispOpInstance.nPoints
		adsStride = 0  #stride to guide to liquid adsorption concentrations i.e., cpp
		unitStride = adsStride + bindStride*nComp*2


		# allocation vectors and matrices
		idx = 1:ConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fjac = Fc # The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cpp = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS_q = zeros(Real, ConvDispOpInstance.nPoints * nComp)
		qq = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS = zeros(Real, adsStride + 2*nComp*bindStride)
		cIn = [0.0]
	
		# if the axial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_ax) == Float64 
			d_ax = ones(Float64,nComp)*d_ax
		end
		
		# Set initial condition vectors 
		c0, cp0, q0 = initial_condition_specification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)
		
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
		new(nComp, colLength, cross_section_area, d_ax, eps_c, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, ConvDispOpInstance, bindStride, adsStride, unitStride, idx, Fc, Fjac, cpp, RHS_q, qq, RHS, solution_outlet, solution_times, bind)
	end
end


# Define a function to compute the transport term for the LRM
function compute_transport!(RHS, RHS_q, cpp, x, m::LRM, t, section, sink, switches, idx_units) 
	# section = i from call 
	# sink is the unit i.e., h from previous call
	
	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)
    
	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	@inbounds for j = 1:m.nComp
		
		# Indices
		# For the indicies regarding mobile phase, + idx_units[sink] must be added to get the right column 
		# For the stationary phase, RHS_q is already a slice of the stationary phase of the right column
		m.idx =  1 + (j-1) * m.ConvDispOpInstance.nPoints : m.ConvDispOpInstance.nPoints + (j-1) * m.ConvDispOpInstance.nPoints

		# Determining inlet concentration 
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t) 

		# Convection Dispersion term
		cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase #
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn[1], m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)

		# Mobile phase RHS 
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dh - m.Fc * RHS_q[m.idx]
	end
	
    nothing
end


################################# LUMPED RATE MODEL WITH PORES (LRMP) #################################
mutable struct LRMP <: ModelBase
	# These parameters are the minimum to be specified for the LRMP
	nComp::Int64 
    colLength::Float64
	cross_section_area::Float64
    d_ax::Union{Float64, Vector{Float64}}
    eps_c::Float64
	eps_p::Float64
	kf::Union{Float64, Vector{Float64}}
	Rp::Float64
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
	cIn::Vector{Float64}
	exactInt::Int64


	polyDeg::Int64
	nCells::Int64	
	
    # Convection Dispersion properties
	ConvDispOpInstance::ConvDispOp

	# based on the input, remaining properties are calculated in the function LRM
	#Determined properties 
	bindStride::Int64
	adsStride::Int64
	unitStride::Int64

	# Allocation vectors and matrices
    idx::UnitRange{Int64}
	idx_p::UnitRange{Int64}
	Fc::Float64
	Fp::Float64
	Fjac::Float64
	cpp::Vector{Float64}
    RHS_q::Vector{Float64}
    qq::Vector{Float64}
	RHS::Vector{Float64}
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	bind::bindingBase

	# Default variables go in the arguments in the LRM(..)
	function LRMP(; nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, exactInt=1, cross_section_area=1.0)
		
		# Get necessary variables for convection dispersion DG 
		ConvDispOpInstance = ConvDispOp(polyDeg,nCells,colLength)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = ConvDispOpInstance.nPoints
		adsStride = nComp*ConvDispOpInstance.nPoints  #stride to guide to liquid adsorption concentrations i.e., cpp
		unitStride = adsStride + bindStride*nComp*2


		# allocation vectors and matrices
		idx = 1:ConvDispOpInstance.nPoints
		idx_p = 1:ConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fp = (1-eps_p)/eps_p
		Fjac = Fp 				# The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cpp = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS_q = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS = zeros(Float64, adsStride + 2*nComp*bindStride)
		qq = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		cIn = [0.0]
	
		# if the axial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_ax) == Float64 
			d_ax = ones(Float64,nComp)*d_ax
		end

		# if the kf is specified for a single component, assume they are the same for all components
		if typeof(kf) == Float64 
			kf = ones(Float64,nComp)*kf
		end
		
		# Set initial condition vectors 
		c0, cp0, q0 = initial_condition_specification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)
		
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
		new(nComp, colLength, cross_section_area, d_ax, eps_c, eps_p, kf, Rp, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, ConvDispOpInstance, bindStride, adsStride, unitStride, idx, idx_p, Fc, Fp, Fjac, cpp, RHS_q, qq, RHS, solution_outlet, solution_times, bind)
	end
end

# Define a function to compute the transport term for the LRMP
function compute_transport!(RHS, RHS_q, cpp, x, m::LRMP, t, section, sink, switches, idx_units)

	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)
	
	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	@inbounds for j = 1:m.nComp

		# Indices
		# For the indicies regarding mobile phase and pore phase, + idx_units[sink] must be added to get the right column 
		# For the stationary phase, RHS_q is already a slice of the stationary phase of the right column
		m.idx =  1 + (j-1) * m.ConvDispOpInstance.nPoints : m.ConvDispOpInstance.nPoints + (j-1) * m.ConvDispOpInstance.nPoints
		m.idx_p = m.idx .+ m.ConvDispOpInstance.nPoints*m.nComp

		
		# Determining inlet concentration 
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t) 
		
		# Convection Dispersion term	
		# m.cpp = @view x[m.idx .+ idx_units[sink]] # mobile phase
		cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase		
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn[1], m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)

		# Mobile phase
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dh - m.Fc * 3 / m.Rp * m.kf[j] * (x[m.idx .+ idx_units[sink]] - x[m.idx_p .+ idx_units[sink]])

		# Pore phase dcp/dt = MT/eps_p - Fp dq/dt
		@. @views RHS[m.idx_p .+ idx_units[sink]] = 3 / m.Rp / m.eps_p * m.kf[j] * (x[m.idx .+ idx_units[sink]]- x[m.idx_p .+ idx_units[sink]]) - m.Fp * RHS_q[m.idx] 

	end
	
	nothing
end


################################# GENERAL RATE MODEL (GRM) #################################

mutable struct GRM <: ModelBase
	# These parameters are the minimum to be specified for the LRMP
	nComp::Int64 
    colLength::Float64
	cross_section_area::Float64
    d_ax::Union{Float64, Vector{Float64}}

    eps_c::Float64
	eps_p::Float64
	kf::Union{Float64, Vector{Float64}}
	Rp::Float64
	Rc::Float64
	Dp::Union{Float64, Vector{Float64}}
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
	cIn::Vector{Float64}
	exactInt::Int64


	polyDeg::Int64
	nCells::Int64	
	polyDegPore::Int64

    # Convection Dispersion properties
	ConvDispOpInstance::ConvDispOp

	# Pore phase GRM properties 
	PoreOpInstance::PoreOp

	# based on the input, remaining properties are calculated in the function LRM
	#Determined properties 
	bindStride::Int64
	adsStride::Int64
	unitStride::Int64

	# Allocation vectors and matrices
    idx::UnitRange{Int64}
	idx_p::StepRange{Int64, Int64}
	idx_q::UnitRange{Int64}
	Fc::Float64
	Fp::Float64
	Fjac::Float64
	Jr::Float64
	invRi::Vector{Float64}
	cpp::Vector{Float64}
    RHS_q::Vector{Float64}
    qq::Vector{Float64}
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	bind::bindingBase

	# Default variables go in the arguments in the LRM(..)
	function GRM(; nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, Dp, Rc=0.0, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, polyDegPore=4, nCells=8, exactInt=1, cross_section_area=1.0)
		
		# Get necessary variables for convection dispersion DG 
		ConvDispOpInstance = ConvDispOp(polyDeg,nCells,colLength)

		# allocation vectors and matrices
		idx = 1:ConvDispOpInstance.nPoints
		idx_p = 1:ConvDispOpInstance.nPoints:2*ConvDispOpInstance.nPoints
		idx_q = 1:ConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fp = (1-eps_p)/eps_p
		Fjac = Fp 				# The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cIn = [0.0]
	
		# if the axial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_ax) == Float64 
			d_ax = ones(Float64,nComp)*d_ax
		end
		# if the kf is specified for a single component, assume they are the same for all components
		if typeof(kf) == Float64 
			kf = ones(Float64,nComp)*kf
		end
		# if the Dp is specified for a single component, assume they are the same for all components
		if typeof(Dp) == Float64 
			Dp = ones(Float64,nComp)*Dp
		end

		#The pore grid goes from -1 to 1, to the physical domain that is
		# r = Rc + 1/2 (xi+1)*(Rp-Rc) - such that xi=-1 -> r=Rc, assuming particle is from [Rc,Rp]
		#Particle radius [m]
		Rc = 0.0 # for this model, particle is set at the Rc=0
		ri = Rc .+ 1/2 .* (DGElements.lglnodes(polyDegPore)[1] .+ 1) .* (Rp-Rc) #DGElements.lglnodes(polyDegPore)[1] = nodesPore

		#Inverse Jacobian of the mapping
		Jr = 2/(Rp-Rc) 
		deltaR = Rp-Rc

		#Inverse ri
		invRi = 1 ./ ri
		invRi[1] = 1.0 #WIll not be used anyway - to avoid inf

		# Get necessary variables for pore phase DG 
		PoreOpInstance = PoreOp(polyDegPore, ConvDispOpInstance.nPoints, nComp, Rc, Rp, deltaR, eps_p)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = ConvDispOpInstance.nPoints*PoreOpInstance.nNodesPore
		adsStride = nComp*ConvDispOpInstance.nPoints  #stride to guide to liquid adsorption concentrations i.e., cpp
		unitStride = adsStride + bindStride*nComp*2
		
		# vector allocations
		cpp = zeros(Float64, PoreOpInstance.stridePore)
		RHS_q = zeros(Float64, PoreOpInstance.stridePore)
		qq = zeros(Float64, PoreOpInstance.stridePore)

		# Solution_outlet as well 
		solution_outlet = zeros(Float64,1,nComp)
		solution_times = Float64[]
		
		# Set initial condition vectors 
		c0, cp0, q0 = initial_condition_specification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)
		
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
		new(nComp, colLength, cross_section_area, d_ax, eps_c, eps_p, kf, Rp, Rc, Dp, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, polyDegPore, ConvDispOpInstance, PoreOpInstance, bindStride, adsStride, unitStride, idx, idx_p, idx_q, Fc, Fp, Fjac, Jr, invRi, cpp, RHS_q, qq, solution_outlet, solution_times, bind)
	end
end



# Define a function to compute the transport term for the GRM
function compute_transport!(RHS, RHS_q, cpp, x, m::GRM, t, section, sink, switches, idx_units)

	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)
	
	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	@inbounds for j = 1:m.nComp

		#Indices
		m.idx =  1 + (j-1) * m.ConvDispOpInstance.nPoints : m.ConvDispOpInstance.nPoints + (j-1) * m.ConvDispOpInstance.nPoints
		
		# Determining inlet concentration 
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t) 

		
		# Convection Dispersion term	
		cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn[1], m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)


		#Surface flux to the particles 
		# idx_p is a step range 
		m.idx_p = m.nComp*m.ConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + m.PoreOpInstance.nNodesPore + + idx_units[sink] : m.PoreOpInstance.nNodesPore : m.nComp*m.ConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + + idx_units[sink]
		m.idx = m.idx .+ idx_units[sink]

		# Mobile phase
		@. @views RHS[m.idx] = m.ConvDispOpInstance.Dh - m.Fc * 3 / m.Rp * m.kf[j] * (x[m.idx] - x[m.idx_p])

		# Pore phase - Each idx has _nPolyPore number of states Adding the boundary flux
		# The boundary flux term is computed as kf(x-x*) 
		# The whole term is L (M^-1) (2/deltaR) Rp^2/eps_p * kf (x-x*)
		@. @views m.PoreOpInstance.boundaryPore =  m.kf[j] * (x[m.idx] - x[m.idx_p])


		#Now the rest is computed
		#dcp/dt = 4/Rp Dp Dr/ri cp - (2/Rp)^2 Dp M^-1 A cp + 2/Rp Dp L[0,J/eps/Dp]  - Fp dq/dt
		@inbounds for i in 1:m.ConvDispOpInstance.nPoints

			#Pore phase for each component starts after all mobile phases
			#through all mobile phase, through pore phase j, at porephase i
			m.idx = 1 + m.nComp*m.ConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + (i-1) * m.PoreOpInstance.nNodesPore + idx_units[sink] : m.PoreOpInstance.nNodesPore + m.nComp*m.ConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints  + (i-1) * m.PoreOpInstance.nNodesPore + idx_units[sink]
			m.idx_q = 1 + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + (i-1) * m.PoreOpInstance.nNodesPore : m.PoreOpInstance.nNodesPore + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints  + (i-1) * m.PoreOpInstance.nNodesPore

			# Term 1 is 
			# (2/deltaR)^2 .* M^-1 A^r Dp cp
			mul!(m.PoreOpInstance.term1, m.PoreOpInstance.invMM_Ar,@view(x[m.idx])) 
			broadcast!(*, m.PoreOpInstance.term1, m.Dp[j], m.PoreOpInstance.term1) # 

			#Pore phase boundary conditions 
			#Corresponds to the boundary term - L (M^-1) (2/deltaR) Rp^2/eps_p * kf (x-x*)
			@. @views RHS[m.idx] = m.PoreOpInstance.boundaryPore[i] * m.PoreOpInstance.LiftMatrixRed 

			#Assembling RHS
			@. @views RHS[m.idx] += - m.PoreOpInstance.term1 - m.Fp * RHS_q[m.idx_q]

		end

	end
	
	nothing
end