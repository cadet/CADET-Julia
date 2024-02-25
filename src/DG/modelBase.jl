





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
	polyDerMPore2Jr::Matrix{Float64}
	polyDerMPoreRed::Matrix{Float64}
	invMMPoreAMatrix::Matrix{Float64}
	nNodesPore::Int64
	boundaryPore::Vector{Float64}
	term1::Vector{Float64}
	term11::Vector{Float64}
	term2::Vector{Float64}

	function PoreOp(polyDegPore,nPoints,nComp,Jr)

		nNodesPore = polyDegPore + 1
		stridePore = nPoints*nNodesPore*nComp #stride pore phase

		# Obtain LGL nodes, weights and matrices for pore phase
		nodesPore,invWeightsPore = DGElements.lglnodes(polyDegPore)              #LGL nodes & weights
		invMMPore = DGElements.invMMatrix(nodesPore,polyDegPore)                 #Inverse mass matrix, M^-1
		polyDerMPore = DGElements.derivativeMatrix(polyDegPore,nodesPore)        #derivative matrix, D

		# Obtain matrices for pore phase
		polyDerMPore2Jr = polyDerMPore * 2 * Jr    								#derivative matrix, D*2*Jr
		polyDerMPoreRed = Transpose(polyDerMPore[1,1:nNodesPore])*Jr             # The reduced D matrix for when r->0 * Jr
		AMatrix = Transpose(polyDerMPore) * inv(invMMPore) * polyDerMPore        #A matrix, D^T M D
		invMMPoreAMatrix = invMMPore * AMatrix*Jr^2                               #M^-1*A which is used in RHS * Jr^2

		# Lifting matrix, L = M^-1 eps
		LiftMatrix = zeros(nNodesPore,2)
		LiftMatrix[1,1] = 1
		LiftMatrix[nNodesPore,2] = 1
		LiftMatrix = invMMPore*LiftMatrix
		LiftMatrixRed = LiftMatrix[:,2] #Since only the second column is used according to the boundary conditions

		boundaryPore = zeros(Float64,nPoints)

		term1 = zeros(nNodesPore)
		term11 = zeros( 1)
		term2 = zeros(nNodesPore)

		new(stridePore, LiftMatrixRed, polyDerMPore2Jr, polyDerMPoreRed, invMMPoreAMatrix, nNodesPore, boundaryPore, term1, term11, term2)
	end
end



################################# TRANSPORT MODELS #################################
abstract type modelBase 
	# From here, the transport models are found
end

################################# LUMPED RATE MODEL (LRM) #################################
mutable struct LRM <: modelBase
	# Check parameters
	# These parameters are the minimum to be specified for the LRM
	nComp::Int64 
    colLength::Float64
    d_ax::Union{Float64, Vector{Float64}}
    eps_c::Float64
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
	
    
	cIn::Float64
	exactInt::Int64

	polyDeg::Int64
	nCells::Int64	
	
    # Convection Dispersion properties
	ConvDispOpInstance::ConvDispOp

	# based on the input, remaining properties are calculated in the function LRM
	#Determined properties 
	bindStride::Int64
	adsStride::Int64

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
	

	# Default variables go in the arguments in the LRM(..)
	function LRM(; nComp, colLength, d_ax, eps_c, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, exactInt=1)
		
		# Get necessary variables for convection dispersion DG 
		ConvDispOpInstance = ConvDispOp(polyDeg,nCells,colLength)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = ConvDispOpInstance.nPoints
		adsStride = 0  #stride to guide to liquid adsorption concentrations i.e., cpp


		# allocation vectors and matrices
		idx = 1:ConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fjac = Fc # The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cpp = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS_q = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		qq = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS = zeros(Float64, adsStride + 2*nComp*bindStride)
		cIn = 0.0
	
		# if the axial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_ax) == Float64 
			d_ax = ones(Float64,nComp)*d_ax
		end
		
		# Set initial condition vectors 
		c0, cp0, q0 = initialConditionSpecification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)
		
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
		new(nComp, colLength, d_ax, eps_c, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, ConvDispOpInstance, bindStride, adsStride, idx, Fc, Fjac, cpp, RHS_q, qq, RHS, solution_outlet, solution_times,bind)
	end
end


# Define a function to compute the transport term for the LRM
function computeTransport!(RHS, RHS_q, x, m::LRM, t, section, sink, switches, idx_units) 
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
		m.cIn = ((switches.connectionInstance.cIn_c[section,sink, j] + 
					switches.connectionInstance.cIn_l[section,sink, j]*t +
					switches.connectionInstance.cIn_q[section,sink, j]*t^2 +
					switches.connectionInstance.cIn_cube[section,sink, j]*t^3) * switches.connectionInstance.u_inlet[switches.switchSetup[section], sink] +
					switches.connectionInstance.u_unit[switches.switchSetup[section], sink] * switches.connectionInstance.c_connect[switches.switchSetup[section], sink, j] * x[switches.connectionInstance.idx_connect[switches.switchSetup[section], sink, j]]) / switches.connectionInstance.u_tot[switches.switchSetup[section], sink]

		# Convection Dispersion term
		m.cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, m.cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.connectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn, m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)

		# Mobile phase RHS 
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dh - m.Fc * RHS_q[m.idx]
	end
	
    nothing
end

# A function that determines the inlet concentrations needed for the transport equations 
# Depends on the inlets specified in the switches
# function inletConcentrations!(cIn, switches, j, switch, sink, x, t, idx_units)
    
	# cIn = ((switches.connectionInstance.cIn_c[switch,sink, j] + 
			# switches.connectionInstance.cIn_l[switch,sink, j]*t +
			# switches.connectionInstance.cIn_q[switch,sink, j]*t^2 +
			# switches.connectionInstance.cIn_cube[switch,sink, j]*t^3) * switches.connectionInstance.u_inlet[switch, sink] +
			# switches.connectionInstance.u_unit[switch, sink] * switches.connectionInstance.c_connect[switch, sink, j] * x[idx_units[sink] + switches.connectionInstance.idx_connect[switch, sink, j]]) / switches.connectionInstance.u_tot[switch, sink]
# end

################################# LUMPED RATE MODEL WITH PORES (LRMP) #################################
mutable struct LRMP <: modelBase
	# These parameters are the minimum to be specified for the LRMP
	nComp::Int64 
    colLength::Float64
    d_ax::Union{Float64, Vector{Float64}}
    eps_c::Float64
	eps_p::Float64
	kf::Union{Float64, Vector{Float64}}
	Rp::Float64
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
	cIn::Float64
	exactInt::Int64


	polyDeg::Int64
	nCells::Int64	
	
    # Convection Dispersion properties
	ConvDispOpInstance::ConvDispOp

	# based on the input, remaining properties are calculated in the function LRM
	#Determined properties 
	bindStride::Int64
	adsStride::Int64

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
	function LRMP(; nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, exactInt=1)
		
		# Get necessary variables for convection dispersion DG 
		ConvDispOpInstance = ConvDispOp(polyDeg,nCells,colLength)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = ConvDispOpInstance.nPoints
		adsStride = nComp*ConvDispOpInstance.nPoints  #stride to guide to liquid adsorption concentrations i.e., cpp


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
		cIn = 0.0
	
		# if the axial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_ax) == Float64 
			d_ax = ones(Float64,nComp)*d_ax
		end

		# if the kf is specified for a single component, assume they are the same for all components
		if typeof(kf) == Float64 
			kf = ones(Float64,nComp)*kf
		end
		
		# Set initial condition vectors 
		c0, cp0, q0 = initialConditionSpecification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)
		
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
		new(nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, ConvDispOpInstance, bindStride, adsStride, idx, idx_p, Fc, Fp, Fjac, cpp, RHS_q, qq, RHS, solution_outlet, solution_times, bind)
	end
end

# Define a function to compute the transport term for the LRMP
function computeTransport!(RHS, RHS_q, x, m::LRMP, t, section, sink, switches, idx_units)
	
	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	@inbounds for j = 1:m.nComp

		# Indices
		# For the indicies regarding mobile phase and pore phase, + idx_units[sink] must be added to get the right column 
		# For the stationary phase, RHS_q is already a slice of the stationary phase of the right column
		m.idx =  1 + (j-1) * m.ConvDispOpInstance.nPoints : m.ConvDispOpInstance.nPoints + (j-1) * m.ConvDispOpInstance.nPoints
		m.idx_p = m.idx .+ m.ConvDispOpInstance.nPoints*m.nComp

		
		# Determining inlet concentration 
		# inletConcentrations!(m.cIn, switches, j, switch, sink, x, t, idx_units) 
		m.cIn = ((switches.connectionInstance.cIn_c[section,sink, j] + 
					switches.connectionInstance.cIn_l[section,sink, j]*t +
					switches.connectionInstance.cIn_q[section,sink, j]*t^2 +
					switches.connectionInstance.cIn_cube[section,sink, j]*t^3) * switches.connectionInstance.u_inlet[switches.switchSetup[section], sink] +
					switches.connectionInstance.u_unit[switches.switchSetup[section], sink] * switches.connectionInstance.c_connect[switches.switchSetup[section], sink, j] * x[switches.connectionInstance.idx_connect[switches.switchSetup[section], sink, j]]) / switches.connectionInstance.u_tot[switches.switchSetup[section], sink]
		
		# Convection Dispersion term	
		# m.cpp = @view x[m.idx .+ idx_units[sink]] # mobile phase
		m.cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase		
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, m.cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.connectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn, m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)

		# Mobile phase
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dh - m.Fc * 3 / m.Rp * m.kf[j] * (x[m.idx .+ idx_units[sink]] - x[m.idx_p .+ idx_units[sink]])

		# Pore phase dcp/dt = MT/eps_p - Fp dq/dt
		@. @views RHS[m.idx_p .+ idx_units[sink]] = 3 / m.Rp / m.eps_p * m.kf[j] * (x[m.idx .+ idx_units[sink]]- x[m.idx_p .+ idx_units[sink]]) - m.Fp * RHS_q[m.idx] 

	end
	
	nothing
end


################################# GENERAL RATE MODEL (GRM) #################################

mutable struct GRM <: modelBase
	# These parameters are the minimum to be specified for the LRMP
	nComp::Int64 
    colLength::Float64
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
	cIn::Float64
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
	function GRM(; nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, Dp, Rc=0.0, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, polyDegPore=4, nCells=8, exactInt=1)
		
		# Get necessary variables for convection dispersion DG 
		ConvDispOpInstance = ConvDispOp(polyDeg,nCells,colLength)

		# allocation vectors and matrices
		idx = 1:ConvDispOpInstance.nPoints
		idx_p = 1:ConvDispOpInstance.nPoints:2*ConvDispOpInstance.nPoints
		idx_q = 1:ConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fp = (1-eps_p)/eps_p
		Fjac = Fp 				# The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cIn = 0.0
	
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

		#Inverse ri
		invRi = 1 ./ ri
		invRi[1] = 1.0 #WIll not be used anyway - to avoid inf

		# Get necessary variables for pore phase DG 
		PoreOpInstance = PoreOp(polyDegPore,ConvDispOpInstance.nPoints,nComp,Jr)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = ConvDispOpInstance.nPoints*PoreOpInstance.nNodesPore
		adsStride = nComp*ConvDispOpInstance.nPoints  #stride to guide to liquid adsorption concentrations i.e., cpp
		
		# vector allocations
		cpp = zeros(Float64, PoreOpInstance.stridePore)
		RHS_q = zeros(Float64, PoreOpInstance.stridePore)
		qq = zeros(Float64, PoreOpInstance.stridePore)

		# Solution_outlet as well 
		solution_outlet = zeros(Float64,1,nComp)
		solution_times = Float64[]
		
		# Set initial condition vectors 
		c0, cp0, q0 = initialConditionSpecification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)
		
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
		new(nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, Rc, Dp, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, polyDegPore, ConvDispOpInstance, PoreOpInstance, bindStride, adsStride, idx, idx_p, idx_q, Fc, Fp, Fjac, Jr, invRi, cpp, RHS_q, qq, solution_outlet, solution_times, bind)
	end
end



# Define a function to compute the transport term for the LRMP
function computeTransport!(RHS, RHS_q, x, m::GRM, t, section, sink, switches, idx_units)
	
	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	@inbounds for j = 1:m.nComp

		#Indices
		m.idx =  1 + (j-1) * m.ConvDispOpInstance.nPoints : m.ConvDispOpInstance.nPoints + (j-1) * m.ConvDispOpInstance.nPoints
		
		# Determining inlet concentration 
		# inletConcentrations!(m.cIn, switches, j, switch, sink, x, t, idx_units) 
		m.cIn = ((switches.connectionInstance.cIn_c[section,sink, j] + 
					switches.connectionInstance.cIn_l[section,sink, j]*t +
					switches.connectionInstance.cIn_q[section,sink, j]*t^2 +
					switches.connectionInstance.cIn_cube[section,sink, j]*t^3) * switches.connectionInstance.u_inlet[switches.switchSetup[section], sink] +
					switches.connectionInstance.u_unit[switches.switchSetup[section], sink] * switches.connectionInstance.c_connect[switches.switchSetup[section], sink, j] * x[switches.connectionInstance.idx_connect[switches.switchSetup[section], sink, j]]) / switches.connectionInstance.u_tot[switches.switchSetup[section], sink]
		
		# Convection Dispersion term	
		m.cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, m.cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.connectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn, m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)


		#Surface flux to the particles 
		# idx_p is a step range 
		m.idx_p = m.nComp*m.ConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + m.PoreOpInstance.nNodesPore : m.PoreOpInstance.nNodesPore : m.nComp*m.ConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints

		# Mobile phase
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dh - m.Fc * 3 / m.Rp * m.kf[j] * (x[m.idx .+ idx_units[sink]] - x[m.idx_p .+ idx_units[sink]])

		#Pore phase - Each idx has _nPolyPore number of states Adding the boundary flux
		#dcp/dt = 2/Rp Dp Dr/ri cp - (2/Rp)^2 Dp M^-1 A cp + 2/Rp Dp L[0,J/eps/Dp]  - Fp dq/dt
		#First the 2/Rp Dp L[0,J/eps/Dp] is added
		@. @views m.PoreOpInstance.boundaryPore =  m.Jr * m.kf[j] / m.eps_p * (x[m.idx .+ idx_units[sink]] - x[m.idx_p .+ idx_units[sink]])


		#Now the rest is computed
		#dcp/dt = 4/Rp Dp Dr/ri cp - (2/Rp)^2 Dp M^-1 A cp + 2/Rp Dp L[0,J/eps/Dp]  - Fp dq/dt
		@inbounds for i in 1:m.ConvDispOpInstance.nPoints

			#Pore phase for each component starts after all mobile phases
			#through all mobile phase, through pore phase j, at porephase i
			m.idx = m.nComp*m.ConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + 1 + (i-1) * m.PoreOpInstance.nNodesPore : m.PoreOpInstance.nNodesPore + m.nComp*m.ConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints  + (i-1) * m.PoreOpInstance.nNodesPore
			
			#Stationary phase starts after all pore phases have been determined
			m.idx_q = m.nComp*m.ConvDispOpInstance.nPoints + m.PoreOpInstance.stridePore + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + 1 + (i-1) * m.PoreOpInstance.nNodesPore : m.PoreOpInstance.nNodesPore + m.nComp*m.ConvDispOpInstance.nPoints + m.PoreOpInstance.stridePore + (j-1) *m.PoreOpInstance.nNodesPore*m.ConvDispOpInstance.nPoints + (i-1) * m.PoreOpInstance.nNodesPore

			#Since the term1 needs special treatment at r->0, it is computed separately
			#we cannot compute 1/r when r=0 but rewriting using L'Hôpital's rule gives
			#1/r dc/dr = d^2/dr^2 for r->0
			# 1/r(xi) dc/dxi 2/Rp = 2/Rp d^2/dr^2 for r-> 0
			#term1 = 1/ri*2 * Jr * Dp * D*c

			mul!(m.PoreOpInstance.term1, m.PoreOpInstance.polyDerMPore2Jr,@view(x[m.idx .+ idx_units[sink]])) #term = 2 * Jr * D*c
			broadcast!(*,m.PoreOpInstance.term1, m.Dp[j] ,m.PoreOpInstance.term1) #term1 = 2 * Jr * Dp * D*c
			mul!(m.PoreOpInstance.term11, m.PoreOpInstance.polyDerMPoreRed,m.PoreOpInstance.term1) #term 11 = 2 * Jr^2 * Dp * D*D*c
			broadcast!(*,m.PoreOpInstance.term1,m.invRi,m.PoreOpInstance.term1) #term1 = 1/ri*2 * Jr * Dp * D*c
			m.PoreOpInstance.term1[1] = m.PoreOpInstance.term11[1] #correction at r->0

			#Pore phase boundary conditions 
			#Corresponds to the boundary term - 2/Rp Dp L[0,J/eps/Dp]
			@. @views RHS[m.idx .+ idx_units[sink]] = m.PoreOpInstance.boundaryPore[i] * m.PoreOpInstance.LiftMatrixRed 

			#Term 2 - Jr^2 .* Dp .* M^-1 A c
			mul!(m.PoreOpInstance.term2, m.PoreOpInstance.invMMPoreAMatrix,@view(x[m.idx .+ idx_units[sink]]))
			broadcast!(*,m.PoreOpInstance.term2,m.Dp[j],m.PoreOpInstance.term2)

			#Assembling RHS
			@. @views RHS[m.idx .+ idx_units[sink]] += m.PoreOpInstance.term1 - m.PoreOpInstance.term2 - m.Fp * RHS[m.idx_q]

		end

	end
	
	nothing
end