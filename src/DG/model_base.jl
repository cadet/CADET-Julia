"""
    ConvDispOp

A struct containing the Discontinuous Galerkin (DG) variables for the convection-dispersion operator, used in LRM, LRMP, and GRM models.

# Fields
- `polyDeg`, `nCells`, `nNodes`, `nPoints`, `strideNode`, `strideCell`: Discretization parameters.
- `nodes`, `invWeights`, `invMM`, `polyDerM`, `deltaZ`: DG node and matrix data.
- `mul1`, `c_star`, `h_star`, `Dc`, `Dh`, `h`: Allocation vectors for intermediate calculations.

# Constructor
- `ConvDispOp(polyDeg, nCells, colLength)`
"""
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
	weights::Vector{Float64}
    invWeights::Vector{Float64}
    invMM::Matrix{Float64}
	MM01::Matrix{Float64}
	MM00::Matrix{Float64}
	rMM::Vector{Matrix{Float64}}
	invrMM::Vector{Matrix{Float64}}
	polyDerM::Matrix{Float64}
	S_g::Vector{Matrix{Float64}}
	d_rad_i::Vector{Float64}  # Dispersion coefficient at each cell interface
    deltarho::Float64
    rho_i::Vector{Float64}
	rho_nodes::Vector{Float64}  # Physical radial coordinates at each DG node

	# Allocation vectors and matrices
	mul1::Vector{Float64}
	mul2::Vector{Float64}
    c_star::Vector{Float64}
    g_star::Vector{Float64}
	Dc::Vector{Float64}
    Dg::Vector{Float64}
    g::Vector{Float64}

	function RadialConvDispOp(polyDeg, nCells, col_inner_radius, col_outer_radius, d_rad)

		nNodes = polyDeg + 1
		nPoints = nNodes * nCells	#Number of cells times number of nodes per component
		strideNode = 1				#Number of points to next concentration in the state vector, here 1
		strideCell = nNodes * strideNode

		nodes, invWeights = DGElements.lglnodes(polyDeg)
		weights = 1 ./ invWeights  # Compute weights from invWeights
		deltarho = (col_outer_radius - col_inner_radius) / nCells
		rho_i  = [col_inner_radius + (cell - 1 ) * deltarho for cell in 1:(nCells + 1)] # face radii

		# Physical radial coordinates at each DG node: ρ̂(ξ) = ρ_i + (Δρ/2)(1 + ξ)
		rho_nodes = zeros(Float64, nPoints)
		for cell in 1:nCells
			for node in 1:nNodes
				rho_nodes[(cell-1) * nNodes + node] = rho_i[cell] + (deltarho/2) * (1 + nodes[node])
			end
		end

		invMM = DGElements.invLagrangeMMatrix(nodes, polyDeg) #Inverse mass matrix
		MM01 = DGElements.LagrangeMMatrix(nodes, polyDeg, 0, 1) # Mass matrix 0,1
		MM00 = DGElements.LagrangeMMatrix(nodes, polyDeg, 0, 0) # Mass matrix 0,0
		rMM = DGElements.weightedMMatrix(nodes, polyDeg, nCells, rho_i, deltarho) # weighted mass matrix
		invrMM = DGElements.invweightedMMatrix(nodes, polyDeg, nCells, rho_i, deltarho) # inverse weighted mass matrix
		polyDerM = DGElements.derivativeMatrix(polyDeg, nodes) #derivative matrix
    	S_g = DGElements.dispMMatrix(nodes, polyDeg, rho_i, deltarho, d_rad, polyDerM, rMM)

		# Compute D values at each cell interface
		if isa(d_rad, Function)
			d_rad_i = [d_rad(rho_i[i]) for i in 1:nCells+1]
		else
			d_rad_i = fill(Float64(d_rad), nCells + 1)
		end

		# allocation vectors and matrices
		mul1 = zeros(Float64, nNodes)
		mul2 = zeros(Float64, nNodes)
		c_star = zeros(Float64, nCells + 1)
		g_star = zeros(Float64, nCells + 1)
		Dc = zeros(Float64, nPoints)
		Dg = zeros(Float64, nPoints)
		g = zeros(Float64, nPoints)

		new(polyDeg, nCells, nNodes, nPoints, strideNode, strideCell, nodes, weights, invWeights, invMM, MM01, MM00, rMM, invrMM, polyDerM, S_g, d_rad_i, deltarho, rho_i, rho_nodes, mul1, mul2, c_star, g_star, Dc, Dg, g)
	end
end

"""
    PoreOp

A struct containing the DG variables for the pore phase operator, used in LRMP and GRM models.

# Fields
- `stridePore`, `LiftMatrixRed`, `invMM_Ar`, `nNodesPore` are DG stuff for the pore phase.
- `term1`, `boundaryPore`: are for allocations.

# Constructor
- `PoreOp(polyDegPore, nPoints, nComp, Rc, Rp, deltaR, eps_p)`
"""
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

"""
    FilmDiffOp

# Fields
- `Q`: Base coefficient (3 / R_p)
- `M_K`: Mass matrices`, `invrMM`, `nComp`: DG matrix data.
- `temp`: buffer

# Constructor
    FilmDiffusion(R_p, k_f, nPoints, nComp, polyDeg, nodes, rho_i, deltarho; quadrature_type=:gauss)
"""
mutable struct FilmDiffOp
	Q::Float64       					   # (3/R_p)
	M_K::Vector{Matrix{Float64}}           # Mass matrices for mass transfer
	invrMM::Vector{Matrix{Float64}}        # Inverse weighted mass matrices
	nComp::Int64                           # Number of components
	temp::Vector{Float64}         		   # Temporary buffer

	function FilmDiffOp(R_p::Float64, k_f::Union{Float64, Function}, nPoints::Int64, nComp::Int64, polyDeg::Int64, nodes::Vector{Float64}, rho_i::Vector{Float64}, deltarho::Float64; quadrature_type::Symbol=:gauss)
		# Compute base mass transfer coefficient
		Q = (3.0 / R_p)
		nCells = length(rho_i) - 1
		rMM = DGElements.weightedMMatrix(nodes, polyDeg, nCells, rho_i, deltarho)
		invrMM = DGElements.invweightedMMatrix(nodes, polyDeg, nCells, rho_i, deltarho)
		M_K = DGElements.filmDiffMMatrix(nodes, polyDeg, rho_i, deltarho, k_f, rMM)

		# Allocate buffer
		temp = zeros(Float64, nPoints * nComp)

		new(Q, M_K, invrMM, nComp, temp)
	end
end

"""
    InletConditions

Abstract type for specifying inlet condition types (static or dynamic) for a unit.
"""
abstract type InletConditions end
# If only the inlet velocity has been specified, it does nothing 
# If dynamic flow rates has been activated, it determines the flow rates based on the time and the values of the coefficients. 

struct DynamicInlets <: InletConditions end
struct StaticInlets <: InletConditions end

"""
    inlet_concentrations!(cIn, switches, j, section, sink, x, t, inlet_type::InletConditions)

Computes the static/dynamic inlet concentration for a given component and section, based on time-dependent coefficients and connections.

# Arguments
- `cIn`: Inlet concentration vector (modified in place).
- `switches`: Switches object.
- `j`: Component index.
- `section`: Section index.
- `sink`: Unit index.
- `x`: State vector.
- `t`: Current time.
- `inlet_type`: Should be `DynamicInlets`.

# Details for dynamic inlets
- Sums contributions from inlets and connected units.
- Divides by total velocity to ensure mass balance.

# Details for static inlets
- Leaves `cIn` unchanged.
"""
function inlet_concentrations!(cIn, switches, j, section, sink, x, t, inlet_type::StaticInlets)
    nothing 
end

function inlet_concentrations!(cIn, switches, j, section, sink, x, t, inlet_type::DynamicInlets)
	# Restarting the inlet concentration
	cIn[j] = 0.0
    
	# Inlets from inlets
	@inbounds for l in eachindex(switches.ConnectionInstance.cIn_c[section][sink][j])
		cIn[j] += (switches.ConnectionInstance.cIn_c[section][sink][j][l] + 
			switches.ConnectionInstance.cIn_l[section][sink][j][l]*t +
			switches.ConnectionInstance.cIn_q[section][sink][j][l]*t^2 +
			switches.ConnectionInstance.cIn_cube[section][sink][j][l]*t^3) * switches.ConnectionInstance.u_inlet[switches.switchSetup[section]][sink][l]
	end

	# Inlets from other units 
	@inbounds for l in eachindex(switches.ConnectionInstance.c_connect[switches.switchSetup[section]][sink][j])
		# println("section = $section, sink = $sink, j = $j, l = $l")
		cIn[j] += switches.ConnectionInstance.u_unit[switches.switchSetup[section]][sink][l] * switches.ConnectionInstance.c_connect[switches.switchSetup[section]][sink][j][l] * x[switches.ConnectionInstance.idx_connect[switches.switchSetup[section]][sink][j][l]]
	end 

	# Divide by total velocity according to mass balance 
	cIn[j] /= switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink]
	nothing
end



################################# TRANSPORT MODELS #################################
"""
    ModelBase

Abstract type for all transport model types (e.g., LRM, LRMP, GRM).
"""
abstract type ModelBase 
end

################################# LUMPED RATE MODEL (LRM) #################################
"""
    LRM <: ModelBase

Lumped Rate Model (LRM) struct for a chromatographic column.

# Fields
- `nComp`, `colLength`, `cross_section_area`, `d_ax`, `eps_c`, `Fc`: Physical and geometric parameters.
- `c0`, `cp0`, `q0`: Initial condition vectors.
- `cIn`, `bind`: Inlet concentration vector and binding model.
- `exactInt`, `polyDeg`, `nCells`, `ConvDispOpInstance`, `Fjac`, `idx`: DG stuff.
- `bindStride`, `adsStride`, `unitStride`: Strides and indexing parameters.
- `cpp`, `RHS_q`, `qq`, `RHS`, `solution_outlet`, `solution_times`: State and allocation variables.

# Constructor
- `LRM(; nComp, colLength, d_ax, eps_c, c0, cp0, q0, polyDeg, nCells, exact_integration, cross_section_area)`
"""
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
	exactInt

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
	function LRM(; nComp, colLength, d_ax, eps_c, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, exact_integration=1, cross_section_area=1.0)
		
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
		cIn = zeros(Float64, nComp)
	
		# if the axial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_ax) == Float64 
			d_ax = ones(Float64,nComp)*d_ax
		end
		
		# Set initial condition vectors 
		c0, cp0, q0 = initial_condition_specification(nComp, ConvDispOpInstance, bindStride, c0, cp0, q0)
		
		# Solution_outlet as well 
		solution_outlet = zeros(Float64,1,nComp)
		solution_times = Float64[]

		# specifying the integration method for the DGSEM
		if exact_integration == 1 || exact_integration == true
			exactInt = ConvDispOperatorDG.exact_integration()
		else
			exactInt = ConvDispOperatorDG.collocation()
		end

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


"""
    compute_transport!(RHS, RHS_q, cpp, x, m::ModelBase, t, section, sink, switches, idx_units)

Computes the transport term for the LRM/LRMP/GRM.

# Arguments
- `RHS`: Vector to store the computed derivatives.
- `RHS_q`: View into `RHS` for stationary phase variables.
- `cpp`: View into `x` for pore phase variables.
- `x`: State vector.
- `m::ModelBase`: The model instance.
- `t`: Current simulation time.
- `section`: Current section index.
- `sink`: Index of the current unit.
- `switches`: Switches object.
- `idx_units`: Vector of starting indices for each unit in the global state vector.

# Details
- Updates inlet concentrations and computes the convection-dispersion term for each component.
- Stores the change in concentration in `RHS`.
"""
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
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j]) 

		# Convection Dispersion term
		cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase #
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn[j], m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)

		# Mobile phase RHS 
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dh - m.Fc * RHS_q[m.idx]
	end
	
    nothing
end


################################# LUMPED RATE MODEL WITH PORES (LRMP) #################################
"""
    LRMP <: ModelBase

Lumped Rate Model with Pores (LRMP) struct for a chromatographic column.

# Fields
- `nComp`, `colLength`, `cross_section_area`, `d_ax`, `Fc`, `Fp`, `Fjac`, `eps_c`, `eps_p`, `kf`, `Rp`: Physical and geometric parameters.
- `c0`, `cp0`, `q0`: Initial condition vectors.
- `cIn`, `bind`: Inlet concentration vector and binding model.
- `exactInt`, `polyDeg`, `nCells`, `ConvDispOpInstance`: DG stuff.
- `bindStride`, `adsStride`, `unitStride`, `idx`, `idx_p`: Stride and indexing parameters.
- `cpp`, `RHS_q`, `qq`, `RHS`, `solution_outlet`, `solution_times`: State and allocation variables.

# Constructor
- `LRMP(; nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, c0, cp0, q0, polyDeg, nCells, exact_integration, cross_section_area)`
"""
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
	exactInt


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
	function LRMP(; nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, exact_integration=1, cross_section_area=1.0)
		
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
		cIn = zeros(Float64, nComp)
	
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

		# specifying the integration method for the DGSEM
		if exact_integration==1 || exact_integration == true
			exactInt = ConvDispOperatorDG.exact_integration()
		else
			exactInt = ConvDispOperatorDG.collocation()
		end
		
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
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j]) 
		
		# Convection Dispersion term	
		# m.cpp = @view x[m.idx .+ idx_units[sink]] # mobile phase
		cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase		
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn[j], m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)

		# Mobile phase
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dh - m.Fc * 3 / m.Rp * m.kf[j] * (x[m.idx .+ idx_units[sink]] - x[m.idx_p .+ idx_units[sink]])

		# Pore phase dcp/dt = MT/eps_p - Fp dq/dt
		@. @views RHS[m.idx_p .+ idx_units[sink]] = 3 / m.Rp / m.eps_p * m.kf[j] * (x[m.idx .+ idx_units[sink]]- x[m.idx_p .+ idx_units[sink]]) - m.Fp * RHS_q[m.idx] 

	end
	
	nothing
end


################################# GENERAL RATE MODEL (GRM) #################################
"""
    GRM <: ModelBase

General Rate Model (GRM) struct for a chromatographic column.

# Fields
- `nComp`, `colLength`, `cross_section_area`, `d_ax`, `eps_c`, `eps_p`, `Fc`, `Fp`, `Fjac`, `kf`, `Rp`, `Rc`, `Dp`: Physical and geometric parameters.
- `c0`, `cp0`, `q0`: Initial condition vectors.
- `cIn`, `bind`: Inlet concentration vector and binding model.
- `exactInt`, `polyDeg`, `nCells`, `Jr`, `invRi`, `polyDegPore`, `ConvDispOpInstance`, `PoreOpInstance`: DG stuff.
- `bindStride`, `adsStride`, `unitStride`, `idx`, `idx_p`, `idx_q`: Stride and indexing parameters.
- `cpp`, `RHS_q`, `qq`, `solution_outlet`, `solution_times`: State and allocation variables.

# Constructor
- `GRM(; nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, Dp, Rc, c0, cp0, q0, polyDeg, polyDegPore, nCells, exact_integration, cross_section_area)`
"""
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
	exactInt


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
	function GRM(; nComp, colLength, d_ax, eps_c, eps_p, kf, Rp, Dp, Rc=0.0, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, polyDegPore=4, nCells=8, exact_integration=1, cross_section_area=1.0)
		
		# Get necessary variables for convection dispersion DG 
		ConvDispOpInstance = ConvDispOp(polyDeg,nCells,colLength)

		# allocation vectors and matrices
		idx = 1:ConvDispOpInstance.nPoints
		idx_p = 1:ConvDispOpInstance.nPoints:2*ConvDispOpInstance.nPoints
		idx_q = 1:ConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fp = (1-eps_p)/eps_p
		Fjac = Fp 				# The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cIn = zeros(Float64, nComp)
	
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

		# specifying the integration method for the DGSEM
		if exact_integration==1 || exact_integration == true
			exactInt = ConvDispOperatorDG.exact_integration()
		else
			exactInt = ConvDispOperatorDG.collocation()
		end
		
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
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j]) 

		
		# Convection Dispersion term	
		cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp] # mobile phase
		ConvDispOperatorDG.residualImpl!(m.ConvDispOpInstance.Dh, cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nPoints, m.ConvDispOpInstance.nNodes, m.nCells, m.ConvDispOpInstance.deltaZ, m.polyDeg, m.ConvDispOpInstance.invWeights, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.d_ax[j], m.cIn[j], m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.h_star, m.ConvDispOpInstance.Dc, m.ConvDispOpInstance.h, m.ConvDispOpInstance.mul1, m.exactInt)


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

################################# RADIAL LUMPED RATE MODEL (rLRM) #################################
"""
    rLRM <: ModelBase

Lumped Rate Model (LRM) struct for a chromatographic radial column.

# Fields
- `nComp`, `col_Rho_c`, `col_Rho`, `colHeight`, `cross_section_area`, `d_rad`, `eps_c`, `Fc`: Physical and geometric parameters.
- `c0`, `cp0`, `q0`: Initial condition vectors.
- `cIn`, `bind`: Inlet concentration vector and binding model.
- `exactInt`, `polyDeg`, `nCells`, `ConvDispOpInstance`, `Fjac`, `idx`: DG stuff.
- `bindStride`, `adsStride`, `unitStride`: Strides and indexing parameters.
- `cpp`, `RHS_q`, `qq`, `RHS`, `solution_outlet`, `solution_times`: State and allocation variables.

# Constructor
- `LRM(; nComp, col_height, d_rad, eps_c, c0, cp0, q0, polyDeg, nCells, exact_integration, cross_section_area)`
"""	
mutable struct rLRM <: ModelBase
	# Check parameters
	# These parameters are the minimum to be specified for the LRM
	nComp::Int64
    col_Rho_c::Float64
    col_Rho::Float64
	col_height::Float64
	cross_section_area::Float64
    d_rad::Union{Float64, Vector{Float64}, Vector{Function}}
    eps_c::Float64
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
    
	cIn::Vector{Float64}

	polyDeg::Int64
	nCells::Int64	
	
    # Convection Dispersion properties
	ConvDispOpInstance::RadialConvDispOp

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
	function rLRM(; nComp, col_Rho_c, col_Rho, d_rad, eps_c, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, col_height=1.0)

		# Normalize d_rad first (before creating ConvDispOpInstance)
		if typeof(d_rad) == Float64
			d_rad_init = d_rad  # Use scalar for S_g initialization
			d_rad = ones(Float64,nComp) * d_rad
		elseif isa(d_rad, Function)
			d_rad_init = d_rad  # Use function for S_g initialization
			d_rad = Vector{Function}([d_rad for _ in 1:nComp])
		else
			d_rad_init = d_rad[1]  # Use first component for S_g initialization
		end

		# Get necessary variables for convection dispersion DG
		ConvDispOpInstance = RadialConvDispOp(polyDeg, nCells, col_Rho_c, col_Rho, d_rad_init)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = ConvDispOpInstance.nPoints
		adsStride = 0  #stride to guide to liquid adsorption concentrations i.e., cpp
		unitStride = adsStride + bindStride*nComp*2

		# allocation vectors and matrices
		idx = 1:ConvDispOpInstance.nPoints
		cross_section_area = 2 * π * col_height
		Fc = (1-eps_c)/eps_c
		Fjac = Fc
		cpp = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS_q = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		qq = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS = zeros(Float64, adsStride + 2*nComp*bindStride)
		cIn = zeros(Float64, nComp)

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
		new(nComp, col_Rho_c, col_Rho, col_height, cross_section_area, d_rad, eps_c, c0, cp0, q0, cIn, polyDeg, nCells, ConvDispOpInstance, bindStride, adsStride, unitStride, idx, Fc, Fjac, cpp, RHS_q, qq, RHS, solution_outlet, solution_times, bind)
	end
end

"""
    compute_transport!(RHS, RHS_q, cpp, x, m::ModelBase, t, section, sink, switches, idx_units)

Computes the transport term for the LRM/LRMP/GRM.

# Arguments
- `RHS`: Vector to store the computed derivatives.
- `RHS_q`: View into `RHS` for stationary phase variables.
- `cpp`: View into `x` for pore phase variables.
- `x`: State vector.
- `m::ModelBase`: The model instance.
- `t`: Current simulation time.
- `section`: Current section index.
- `sink`: Index of the current unit.
- `switches`: Switches object.
- `idx_units`: Vector of starting indices for each unit in the global state vector.

# Details
- Updates inlet concentrations and computes the convection-dispersion term for each component.
- Stores the change in concentration in `RHS`.
"""
# Define a function to compute the transport term for the rLRM
function compute_transport!(RHS, RHS_q, cpp, x, m::rLRM, t, section, sink, switches, idx_units)
	# section = i from call
	# sink is the unit i.e., h from previous call

	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)

	@inbounds for j = 1:m.nComp

		# Indices
		# For the indicies regarding mobile phase, + idx_units[sink] must be added to get the right column
		# For the stationary phase, RHS_q is already a slice of the stationary phase of the right column
		m.idx =  1 + (j-1) * m.ConvDispOpInstance.nPoints : m.ConvDispOpInstance.nPoints + (j-1) * m.ConvDispOpInstance.nPoints

		# inlet conc for comp j
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j])

		# S_g and d_rad_i are precomputed in constructor (assumes same d_rad for all components)

		# Convection Dispersion term
		cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp]
		RadialConvDispOperatorDG.radialresidualImpl!(m.ConvDispOpInstance.Dc, cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nNodes, m.ConvDispOpInstance.nCells, m.ConvDispOpInstance.deltarho, m.polyDeg, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, m.ConvDispOpInstance.MM01, m.ConvDispOpInstance.MM00, m.ConvDispOpInstance.rMM, m.ConvDispOpInstance.invrMM, m.ConvDispOpInstance.S_g, m.ConvDispOpInstance.nodes, m.ConvDispOpInstance.weights, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.ConvDispOpInstance.d_rad_i, m.ConvDispOpInstance.rho_i, m.cIn[j], m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.g_star, m.ConvDispOpInstance.Dg, m.ConvDispOpInstance.g, m.ConvDispOpInstance.mul1, m.ConvDispOpInstance.mul2)

		# Mobile phase RHS
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dc - m.Fc * RHS_q[m.idx]
	end

    nothing
end

################################# RADIAL LUMPED RATE MODEL WITH PORES (rLRMP) #################################
"""
    rLRMP <: ModelBase

Lumped Rate Model with Pores (rLRM) struct for a chromatographic radial column.

# Fields
- `nComp`, `col_Rho_c`, `col_Rho`, `colHeight`, `cross_section_area`, `d_rad`, `eps_c`, `Fc`: Physical and geometric parameters.
- `c0`, `cp0`, `q0`: Initial condition vectors.
- `cIn`, `bind`: Inlet concentration vector and binding model.
- `exactInt`, `polyDeg`, `nCells`, `ConvDispOpInstance`, `Fjac`, `idx`: DG stuff.
- `bindStride`, `adsStride`, `unitStride`: Strides and indexing parameters.
- `cpp`, `RHS_q`, `qq`, `RHS`, `solution_outlet`, `solution_times`: State and allocation variables.

# Constructor
- `LRMP(; nComp, col_height, d_rad, eps_c, c0, cp0, q0, polyDeg, nCells, exact_integration, cross_section_area)`
"""	
mutable struct rLRMP <: ModelBase
	# Check parameters
	# These parameters are the minimum to be specified for the rLRMP
	nComp::Int64
    col_Rho_c::Float64
    col_Rho::Float64
	col_height::Float64
    d_rad::Union{Float64, Vector{Float64}, Vector{Function}}
    eps_c::Float64
	eps_p::Float64
	kf::Union{Float64, Vector{Float64}, Vector{Function}}
	Rp::Float64
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0

	cIn::Vector{Float64}

	polyDeg::Int64
	nCells::Int64

    # Convection Dispersion properties
	ConvDispOpInstance::RadialConvDispOp
	FilmDiffOpInstance::FilmDiffOp

	# based on the input, remaining properties are calculated in the function LRMP
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
	
	

	# Default variables go in the arguments in the rLRMP
	function rLRMP(; nComp, col_Rho_c, col_Rho, d_rad, eps_c, eps_p, kf, Rp, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, col_height=1.0)

		# Normalize d_rad first (before creating ConvDispOpInstance)
		if typeof(d_rad) == Float64
			d_rad_init = d_rad  # Use scalar for S_g initialization
			d_rad = ones(Float64,nComp) * d_rad
		elseif isa(d_rad, Function)
			d_rad_init = d_rad  # Use function for S_g initialization
			d_rad = Vector{Function}([d_rad for _ in 1:nComp])
		else
			d_rad_init = d_rad[1]  # Use first component for S_g initialization
		end

		# Get necessary variables for convection dispersion DG
		ConvDispOpInstance = RadialConvDispOp(polyDeg, nCells, col_Rho_c, col_Rho, d_rad_init)

		# The bind stride is the stride between each component for binding. For LRMP, it is nPoints=(polyDeg + 1) * nCells
		bindStride = ConvDispOpInstance.nPoints
		adsStride = nComp * ConvDispOpInstance.nPoints  #stride to guide to liquid adsorption concentrations i.e., cpp
		unitStride = adsStride + bindStride*nComp*2

		# allocation vectors and matrices
		idx = 1:ConvDispOpInstance.nPoints
		idx_p = 1:ConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fp = (1-eps_p)/eps_p
		Fjac = Fp  # The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cpp = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS_q = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		qq = zeros(Float64, ConvDispOpInstance.nPoints * nComp)
		RHS = zeros(Float64, adsStride + 2*nComp*bindStride)
		cIn = zeros(Float64, nComp)

		# if the kf is specified for a single component, assume they are the same for all components
		if typeof(kf) == Float64
			kf = ones(Float64,nComp) * kf
		elseif isa(kf, Function)
			kf = Vector{Function}([kf for _ in 1:nComp])
		end

		# Initialize FilmDiffOp for radial geometry
		# Note: For now, we'll create this with a scalar kf (or first component)
		# Individual components will recompute their M_K matrices in compute_transport!
		kf_init = isa(kf, Vector) ? kf[1] : kf
		FilmDiffOpInstance = FilmDiffOp(Rp, kf_init, ConvDispOpInstance.nPoints, nComp, polyDeg, ConvDispOpInstance.nodes, ConvDispOpInstance.rho_i, ConvDispOpInstance.deltarho)

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
		new(nComp, col_Rho_c, col_Rho, col_height, d_rad, eps_c, eps_p, kf, Rp, c0, cp0, q0, cIn, polyDeg, nCells, ConvDispOpInstance, FilmDiffOpInstance, bindStride, adsStride, unitStride, idx, idx_p, Fc, Fp, Fjac, cpp, RHS_q, qq, RHS, solution_outlet, solution_times, bind)
	end
end

"""
    compute_transport!(RHS, RHS_q, cpp, x, m::ModelBase, t, section, sink, switches, idx_units)

Computes the transport term for the LRM/LRMP/GRM.

# Arguments
- `RHS`: Vector to store the computed derivatives.
- `RHS_q`: View into `RHS` for stationary phase variables.
- `cpp`: View into `x` for pore phase variables.
- `x`: State vector.
- `m::ModelBase`: The model instance.
- `t`: Current simulation time.
- `section`: Current section index.
- `sink`: Index of the current unit.
- `switches`: Switches object.
- `idx_units`: Vector of starting indices for each unit in the global state vector.

# Details
- Updates inlet concentrations and computes the convection-dispersion term for each component.
- Stores the change in concentration in `RHS`.
"""
# Define a function to compute the transport term for the rLRMP
function compute_transport!(RHS, RHS_q, cpp, x, m::rLRMP, t, section, sink, switches, idx_units)
	# section = i from call
	# sink is the unit i.e., h from previous call

	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)

	@inbounds for j = 1:m.nComp

		# Indices
		# For the indicies regarding mobile phase and pore phase, + idx_units[sink] must be added to get the right column
		# For the stationary phase, RHS_q is already a slice of the stationary phase of the right column
		m.idx =  1 + (j-1) * m.ConvDispOpInstance.nPoints : m.ConvDispOpInstance.nPoints + (j-1) * m.ConvDispOpInstance.nPoints
		m.idx_p = m.idx .+ m.ConvDispOpInstance.nPoints*m.nComp

		# inlet conc for comp j
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j])

		# S_g and d_rad_i are precomputed in constructor (assumes same d_rad for all components)
		# M_K is precomputed in FilmDiffOp constructor (assumes same kf for all components)

		# Convection Dispersion term
		cpp = @view x[1 + idx_units[sink] : idx_units[sink] + m.ConvDispOpInstance.nPoints * m.nComp]
		RadialConvDispOperatorDG.radialresidualImpl!(m.ConvDispOpInstance.Dc, cpp, m.idx, m.ConvDispOpInstance.strideNode, m.ConvDispOpInstance.strideCell, m.ConvDispOpInstance.nNodes, m.ConvDispOpInstance.nCells, m.ConvDispOpInstance.deltarho, m.polyDeg, m.ConvDispOpInstance.polyDerM, m.ConvDispOpInstance.invMM, m.ConvDispOpInstance.MM01, m.ConvDispOpInstance.MM00, m.ConvDispOpInstance.rMM, m.ConvDispOpInstance.invrMM, m.ConvDispOpInstance.S_g, m.ConvDispOpInstance.nodes, m.ConvDispOpInstance.weights, switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.ConvDispOpInstance.d_rad_i, m.ConvDispOpInstance.rho_i, m.cIn[j], m.ConvDispOpInstance.c_star, m.ConvDispOpInstance.g_star, m.ConvDispOpInstance.Dg, m.ConvDispOpInstance.g, m.ConvDispOpInstance.mul1, m.ConvDispOpInstance.mul2)

		# Film diffusion term with radial geometry weighting
		# Compute M_ρ^{-1} * M_K * (c - cp) for each cell
		fill!(m.FilmDiffOpInstance.temp, 0.0)
		for cell in 1:m.ConvDispOpInstance.nCells
			cell_idx = (cell-1)*m.ConvDispOpInstance.nNodes + 1 : cell*m.ConvDispOpInstance.nNodes
			c_cell = @view x[m.idx[cell_idx] .+ idx_units[sink]]
			cp_cell = @view x[m.idx_p[cell_idx] .+ idx_units[sink]]
			# temp = M_K * (c - cp)
			temp_diff = c_cell .- cp_cell
			temp_MK = m.FilmDiffOpInstance.M_K[cell] * temp_diff
			# Apply M_ρ^{-1} and geometric factor
			m.FilmDiffOpInstance.temp[cell_idx] .= m.FilmDiffOpInstance.Q .* (m.ConvDispOpInstance.invrMM[cell] * temp_MK)
		end

		# Mobile phase: dc/dt = Dc - Fc * filmDiff
		@. @views RHS[m.idx .+ idx_units[sink]] = m.ConvDispOpInstance.Dc - m.Fc * m.FilmDiffOpInstance.temp[m.idx]

		# Pore phase: dcp/dt = filmDiff/eps_p - Fp * dq/dt
		@. @views RHS[m.idx_p .+ idx_units[sink]] = m.FilmDiffOpInstance.temp[m.idx] / m.eps_p - m.Fp * RHS_q[m.idx]
	end

    nothing
end