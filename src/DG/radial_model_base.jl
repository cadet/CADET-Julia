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
    MM::Matrix{Float64}
    polyDerM::Matrix{Float64}
	invMrhoM::Vector{Matrix{Float64}}
    SgMatrix::Union{Nothing, Vector{Matrix{Float64}}}
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
    Mr_mass_r::Vector{Matrix{Float64}}

	function RadialConvDispOp(polyDeg, nCells, col_inner_radius, col_outer_radius; d_rad_const::Union{Nothing,Float64}=nothing)

		nNodes = polyDeg + 1
		nPoints = nNodes * nCells	#Number of cells times number of nodes per component
		strideNode = 1				#Number of points to next concentration in the state vector, here 1
		strideCell = nNodes * strideNode

		# delta rho [m]
		deltarho = (col_outer_radius - col_inner_radius) / nCells

		# Precompute face radii (left faces)
		rho_i  = [col_inner_radius + (cell-1) * deltarho for cell in 1:nCells]
        # Precompute right face radii
        rho_ip1 = [col_inner_radius + cell * deltarho for cell in 1:nCells]

		# Obtain LGL nodes and weights
		nodes, invWeights = DGElements.lglnodes(polyDeg) #LGL nodes & weights
		invMM = DGElements.invMMatrix(nodes, polyDeg) #Inverse mass matrix
		MM = DGElements.MMatrix(nodes, polyDeg)
		polyDerM = DGElements.derivativeMatrix(polyDeg, nodes) #derivative matrix
		invMrhoM = [
			DGElements.invMrhoMatrix(nodes, polyDeg, deltarho, col_inner_radius + (cell-1)*deltarho)
		    for cell in 1:nCells
		]

		if d_rad_const === nothing
			SgMatrix = nothing
		else
			SgMatrix = [
				DGElements.weighted_stiff_Matrix(nodes, polyDeg, col_inner_radius + (cell-1)*deltarho, deltarho, _ -> d_rad_const)
				for cell in 1:nCells
			]
		end

		# allocation vectors and matrices
		mul1 = zeros(Float64, nNodes)
		c_star = zeros(Float64, nCells + 1)
		g_star = zeros(Float64, nCells + 1)
		Dc = zeros(Float64, nPoints)
		Dg = zeros(Float64, nPoints)
		h = zeros(Float64, nPoints)

        # Cache per-cell weighted mass for w(r)=ρ. Later we scale by k_f if constant.
        Mr_mass_r = [
            DGElements.weightedQuadrature(nodes, col_inner_radius + (cell-1)*deltarho, deltarho, r -> r)
            for cell in 1:nCells
        ]

		new(polyDeg, nCells, nNodes, nPoints, strideNode, strideCell, nodes, invWeights, invMM, MM, polyDerM, invMrhoM, SgMatrix, deltarho, rho_i, rho_ip1, mul1, c_star, g_star, Dc, Dg, h, Mr_mass_r)
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
	col_height::Float64
	cross_section_area::Float64
    d_rad::Union{Float64, Vector{Float64}}
    eps_c::Float64
	c0::Union{Float64, Vector{Float64}} # defaults to 0
	cp0::Union{Float64, Vector{Float64}} # if not specified, defaults to c0
	q0::Union{Float64, Vector{Float64}} # defaults to 0
	
    
	cIn::Vector{Float64}
	exactInt

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
	function rLRM(; nComp, col_inner_radius, col_outer_radius, col_height, d_rad, eps_c, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, exact_integration=1, cross_section_area=1.0)
		
		# Get necessary variables for convection dispersion DG 
		RadialConvDispOpInstance = RadialConvDispOp(polyDeg, nCells, col_inner_radius, col_outer_radius, d_rad_const = (d_rad isa Float64 ? d_rad : nothing))

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
	
		# if the radial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_rad) == Float64 
			d_rad = ones(Float64,nComp)*d_rad
		end
		
		# Set initial condition vectors 
		c0, cp0, q0 = initial_condition_specification(nComp, RadialConvDispOpInstance, bindStride, c0, cp0, q0)
		
		# Solution_outlet as well 
		solution_outlet = zeros(Float64,1,nComp)
		solution_times = Float64[]

		exactInt = RadialConvDispOperatorDG.exact_integration()

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
		new(nComp, col_inner_radius, col_outer_radius, col_height, cross_section_area, d_rad, eps_c, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, RadialConvDispOpInstance, bindStride, adsStride, unitStride, idx, Fc, Fjac, cpp, RHS_q, qq, RHS, solution_outlet, solution_times, bind)
	end
end


# Define a function to compute the transport term for the rLRM
function compute_transport!(RHS, RHS_q, cpp, x, m::rLRM, t, section, sink, switches, idx_units) 
	# section = i from call 
	# sink is the unit i.e., h from previous call
	
	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)
    
	# Precompute face velocities once per call (same for all components)
	faces_v = RadialConvDispOperatorDG.compute_faces_v(switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.RadialConvDispOpInstance.rho_i, m.RadialConvDispOpInstance.rho_ip1, m.col_inner_radius, m.RadialConvDispOpInstance.nCells)
	
	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	@inbounds for j = 1:m.nComp
		# Indices
		# For the indicies regarding mobile phase, + idx_units[sink] must be added to get the right column 
		# For the stationary phase, RHS_q is already a slice of the stationary phase of the right column
		m.idx =  1 + (j-1) * m.RadialConvDispOpInstance.nPoints : m.RadialConvDispOpInstance.nPoints + (j-1) * m.RadialConvDispOpInstance.nPoints

		# Determining inlet concentration 
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j]) 

		# Convection Dispersion term
		cpp_block = @view x[1 + idx_units[sink] : idx_units[sink] + m.RadialConvDispOpInstance.nPoints * m.nComp]

		# per-component diffusion scales (m.d_rad[j] can be scalar here)
		D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(m.d_rad[j], m.RadialConvDispOpInstance.rho_i, m.RadialConvDispOpInstance.rho_ip1)
		left_scale_vec  = @. m.RadialConvDispOpInstance.rho_i  * D_left
		right_scale_vec = @. m.RadialConvDispOpInstance.rho_ip1 * D_right

		RadialConvDispOperatorDG.radialresidualImpl!(
			m.RadialConvDispOpInstance.Dc,
			cpp_block, m.idx,
			m.RadialConvDispOpInstance.strideNode,
			m.RadialConvDispOpInstance.strideCell,
			m.RadialConvDispOpInstance.nNodes,
			m.RadialConvDispOpInstance.nCells,
			m.RadialConvDispOpInstance.deltarho,
			m.polyDeg,
			m.RadialConvDispOpInstance.invWeights,
			m.RadialConvDispOpInstance.nodes,
			m.RadialConvDispOpInstance.polyDerM,
			m.RadialConvDispOpInstance.invMM,
			m.RadialConvDispOpInstance.MM,
			m.RadialConvDispOpInstance.invMrhoM,
			m.RadialConvDispOpInstance.SgMatrix,
			switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink],
			m.d_rad[j],
			m.cIn[j],
			m.RadialConvDispOpInstance.c_star,
			m.RadialConvDispOpInstance.g_star,
			m.RadialConvDispOpInstance.Dg,
			m.RadialConvDispOpInstance.h,
			m.RadialConvDispOpInstance.mul1,
			m.exactInt,
			m.RadialConvDispOpInstance.rho_i,
			m.RadialConvDispOpInstance.rho_ip1,
			faces_v,
			left_scale_vec,
			right_scale_vec
		)

		# Mobile phase RHS 
		@. @views RHS[m.idx .+ idx_units[sink]] = m.RadialConvDispOpInstance.Dc
	end
	
    nothing
end


################################# LUMPED RATE MODEL WITH PORES (rLRMP) #################################
mutable struct rLRMP <: RadialModelBase
	# These parameters are the minimum to be specified for the LRMP
	nComp::Int64 
    col_inner_radius::Float64
    col_outer_radius::Float64
	col_height::Float64
	cross_section_area::Float64
    d_rad::Union{Float64, Vector{Float64}}
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
	RadialConvDispOpInstance::RadialConvDispOp

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
	film_vec::Vector{Float64}
	RHS::Vector{Float64}
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	bind::bindingBase

	# Default variables go in the arguments in the LRM
	function rLRMP(; nComp, col_inner_radius, col_outer_radius, col_height, d_rad, eps_c, eps_p, kf, Rp, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, nCells=8, exact_integration=1, cross_section_area=1.0)
		
		# Get necessary variables for convection dispersion DG 
		RadialConvDispOpInstance = RadialConvDispOp(polyDeg, nCells, col_inner_radius, col_outer_radius, d_rad_const = (d_rad isa Float64 ? d_rad : nothing))

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = RadialConvDispOpInstance.nPoints
		adsStride = nComp*RadialConvDispOpInstance.nPoints  #stride to guide to liquid adsorption concentrations i.e., cpp
		unitStride = adsStride + bindStride*nComp*2

		# allocation vectors and matrices
		idx = 1:RadialConvDispOpInstance.nPoints
		idx_p = 1:RadialConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fp = (1-eps_p)/eps_p
		Fjac = Fp 				# The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cpp = zeros(Float64, RadialConvDispOpInstance.nPoints * nComp)
		RHS_q = zeros(Float64, RadialConvDispOpInstance.nPoints * nComp)
		RHS = zeros(Float64, adsStride + 2*nComp*bindStride)
		qq = zeros(Float64, RadialConvDispOpInstance.nPoints * nComp)
		film_vec = zeros(Float64, RadialConvDispOpInstance.nPoints)
		cIn = zeros(Float64, nComp)
	
		# if the axial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_rad) == Float64 
			d_rad = ones(Float64,nComp)*d_rad
		end

		# if the kf is specified for a single component, assume they are the same for all components
		if typeof(kf) == Float64 
			kf = ones(Float64,nComp)*kf
		end
		
		# Set initial condition vectors 
		c0, cp0, q0 = initial_condition_specification(nComp, RadialConvDispOpInstance, bindStride, c0, cp0, q0)

		# specifying the integration method for the DGSEM
		exactInt = RadialConvDispOperatorDG.exact_integration()
		
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
		new(nComp, col_inner_radius, col_outer_radius, col_height, cross_section_area, d_rad, eps_c, eps_p, kf, Rp, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, RadialConvDispOpInstance, bindStride, adsStride, unitStride, idx, idx_p, Fc, Fp, Fjac, cpp, RHS_q, qq, film_vec, RHS, solution_outlet, solution_times, bind)
	end
end

# Define a function to compute the transport term for the LRMP
function compute_transport!(RHS, RHS_q, cpp, x, m::rLRMP, t, section, sink, switches, idx_units)

	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)
	
	# Precompute face velocities once per call (same for all components)
	faces_v = RadialConvDispOperatorDG.compute_faces_v(switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink], m.RadialConvDispOpInstance.rho_i, m.RadialConvDispOpInstance.rho_ip1, m.col_inner_radius, m.RadialConvDispOpInstance.nCells)

	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	@inbounds for j = 1:m.nComp
		# Indices
		# For the indicies regarding mobile phase and pore phase, + idx_units[sink] must be added to get the right column 
		# For the stationary phase, RHS_q is already a slice of the stationary phase of the right column
		m.idx =  1 + (j-1) * m.RadialConvDispOpInstance.nPoints : m.RadialConvDispOpInstance.nPoints + (j-1) * m.RadialConvDispOpInstance.nPoints
		m.idx_p = m.idx .+ m.RadialConvDispOpInstance.nPoints*m.nComp

		# Determining inlet concentration 
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j]) 
		# Convection Dispersion term	
		# m.cpp = @view x[m.idx .+ idx_units[sink]] # mobile phase
		cpp_block = @view x[1 + idx_units[sink] : idx_units[sink] + m.RadialConvDispOpInstance.nPoints * m.nComp]

		D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(m.d_rad[j], m.RadialConvDispOpInstance.rho_i, m.RadialConvDispOpInstance.rho_ip1)
		left_scale_vec  = @. m.RadialConvDispOpInstance.rho_i  * D_left
		right_scale_vec = @. m.RadialConvDispOpInstance.rho_ip1 * D_right

		RadialConvDispOperatorDG.radialresidualImpl!(
			m.RadialConvDispOpInstance.Dc,
			cpp_block, m.idx,
			m.RadialConvDispOpInstance.strideNode,
			m.RadialConvDispOpInstance.strideCell,
			m.RadialConvDispOpInstance.nNodes,
			m.RadialConvDispOpInstance.nCells,
			m.RadialConvDispOpInstance.deltarho,
			m.polyDeg,
			m.RadialConvDispOpInstance.invWeights,
			m.RadialConvDispOpInstance.nodes,
			m.RadialConvDispOpInstance.polyDerM,
			m.RadialConvDispOpInstance.invMM,
			m.RadialConvDispOpInstance.MM,
			m.RadialConvDispOpInstance.invMrhoM,
			m.RadialConvDispOpInstance.SgMatrix,
			switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink],
			m.d_rad[j],
			m.cIn[j],
			m.RadialConvDispOpInstance.c_star,
			m.RadialConvDispOpInstance.g_star,
			m.RadialConvDispOpInstance.Dg,
			m.RadialConvDispOpInstance.h,
			m.RadialConvDispOpInstance.mul1,
			m.exactInt,
			m.RadialConvDispOpInstance.rho_i,
			m.RadialConvDispOpInstance.rho_ip1,
			faces_v,
			left_scale_vec,
			right_scale_vec
		)

        # Matrix-weighted film transfer using M_rho^{-1} M_K (see derivation eqs. (14a), (17a))
        let nN = m.RadialConvDispOpInstance.nNodes,
            nC = m.RadialConvDispOpInstance.nCells,
            nodes = m.RadialConvDispOpInstance.nodes,
            invMrhoM = m.RadialConvDispOpInstance.invMrhoM,
            deltarho = m.RadialConvDispOpInstance.deltarho,
            rho_i_vec = m.RadialConvDispOpInstance.rho_i

            # Temporary work vectors per cell
            tmp_cell = similar(m.RadialConvDispOpInstance.mul1)
            diff_cell = similar(m.RadialConvDispOpInstance.mul1)

            for cell in 1:nC
                # Slices for this component & cell in the big state vector
                cell_start = m.idx.start + (cell-1)*nN
                cell_end   = cell_start + nN - 1

                c_cell  = @view x[cell_start + idx_units[sink] : cell_end + idx_units[sink]]
                cp_cell = @view x[(cell_start + m.RadialConvDispOpInstance.nPoints*m.nComp) + idx_units[sink] : (cell_end   + m.RadialConvDispOpInstance.nPoints*m.nComp) + idx_units[sink]]

                # diff = (c - cp)
                @. diff_cell = c_cell - cp_cell

                kf_val = m.kf[j]
                if !(kf_val isa Function)
                    # Reuse cached M_r (weight = ρ) and scale by constant k_f
                    MK_r = m.RadialConvDispOpInstance.Mr_mass_r[cell]
                    mul!(tmp_cell, MK_r, diff_cell)
                    @. tmp_cell = kf_val * tmp_cell
                else
                    # Fallback: build full M_K with ρ * k_f(r)
                    MK_full = RadialConvDispOperatorDG.film_mass_matrix(nodes, m.RadialConvDispOpInstance.nNodes,
                                                                        rho_i_vec[cell], deltarho, kf_val)
                    mul!(tmp_cell, MK_full, diff_cell)
                end
                # tmp_cell = M_rho^{-1} * tmp_cell
                mul!(tmp_cell, invMrhoM[cell], tmp_cell)

                # Coefficients
                bulk_coeff = m.Fc * 3.0 / m.Rp                 # (1-ε_c)/ε_c * 3/Rp
                pore_coeff = 3.0 / (m.Rp * m.eps_p)            # 3/(Rp ε_p)

                # Accumulate into RHS: bulk -= bulk_coeff * tmp; pore += pore_coeff * tmp
                @views @. RHS[cell_start + idx_units[sink] : cell_end + idx_units[sink]] = m.RadialConvDispOpInstance.Dc[cell_start - m.idx.start + 1 : cell_end - m.idx.start + 1] - bulk_coeff * tmp_cell

                @views @. RHS[(cell_start + m.RadialConvDispOpInstance.nPoints*m.nComp) + idx_units[sink] : (cell_end   + m.RadialConvDispOpInstance.nPoints*m.nComp) + idx_units[sink]] = pore_coeff * tmp_cell - m.Fp * RHS_q[m.idx][cell_start - m.idx.start + 1 : cell_end - m.idx.start + 1]
            end
        end
	end
	
	nothing
end


################################# GENERAL RATE MODEL (GRM) #################################

mutable struct rGRM <: RadialModelBase
	# These parameters are the minimum to be specified for the LRMP
	nComp::Int64
	col_inner_radius::Float64
    col_outer_radius::Float64
	col_height::Float64
	cross_section_area::Float64
    d_rad::Union{Float64, Vector{Float64}}

    eps_c::Float64
	eps_p::Float64
	kf::Union{Float64, Vector{Float64}}
	Rp::Float64
	Rc::Float64
	Dp::Union{Float64, Vector{Float64}}
	c0::Union{Float64, Vector{Float64}}
	cp0::Union{Float64, Vector{Float64}}
	q0::Union{Float64, Vector{Float64}}
	cIn::Vector{Float64}
	exactInt


	polyDeg::Int64
	nCells::Int64	
	polyDegPore::Int64

    # Convection Dispersion properties
	RadialConvDispOpInstance::RadialConvDispOp

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
	film_vec::Vector{Float64}
	solution_outlet::Matrix{Float64}
	solution_times::Vector{Float64}
	bind::bindingBase

	# Default variables go in the arguments in the LRM(..)
	function rGRM(; nComp, col_inner_radius, col_outer_radius, col_height, d_rad, eps_c, eps_p, kf, Rp, Dp, Rc=0.0, c0 = 0.0, cp0 = -1, q0 = 0, polyDeg=4, polyDegPore=4, nCells=8, exact_integration=1, cross_section_area=1.0)
		
		# Get necessary variables for convection dispersion DG 
		RadialConvDispOpInstance = RadialConvDispOp(polyDeg, nCells, col_inner_radius, col_outer_radius, d_rad_const = (d_rad isa Float64 ? d_rad : nothing))

		# allocation vectors and matrices
		idx = 1:RadialConvDispOpInstance.nPoints
		idx_p = 1:RadialConvDispOpInstance.nPoints:2*RadialConvDispOpInstance.nPoints
		idx_q = 1:RadialConvDispOpInstance.nPoints
		Fc = (1-eps_c)/eps_c
		Fp = (1-eps_p)/eps_p
		Fjac = Fp 				# The phase ratio used for Jacobian i.e., dcdc = Fjac * dqdc
		cIn = zeros(Float64, nComp)
	
		# if the radial dispersion is specified for a single component, assume they are the same for all components
		if typeof(d_rad) == Float64 
			d_rad = ones(Float64,nComp)*d_rad
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
		ri = Rc .+ 1/2 .* (DGElements.lglnodes(polyDegPore)[1] .+ 1) .* (Rp-Rc) #DGElements.lglnodes(polyDegPore)[1] = nodesPore

		#Inverse Jacobian of the mapping
		Jr = 2/(Rp-Rc) 
		deltarho = Rp-Rc

		#Inverse ri
		invRi = 1 ./ ri
		invRi[1] = 1.0

		# Get necessary variables for pore phase DG 
		PoreOpInstance = PoreOp(polyDegPore, RadialConvDispOpInstance.nPoints, nComp, Rc, Rp, deltarho, eps_p)

		# The bind stride is the stride between each component for binding. For LRM, it is nPoints=(polyDeg + 1) * nCells
		bindStride = RadialConvDispOpInstance.nPoints*PoreOpInstance.nNodesPore
		adsStride = nComp*RadialConvDispOpInstance.nPoints  #stride to guide to liquid adsorption concentrations i.e., cpp
		unitStride = adsStride + bindStride*nComp*2
		
		# vector allocations
		cpp = zeros(Float64, PoreOpInstance.stridePore)
		RHS_q = zeros(Float64, PoreOpInstance.stridePore)
		qq = zeros(Float64, PoreOpInstance.stridePore)
		film_vec = zeros(Float64, RadialConvDispOpInstance.nPoints)

		# Solution_outlet as well 
		solution_outlet = zeros(Float64,1,nComp)
		solution_times = Float64[]
		
		# Set initial condition vectors 
		c0, cp0, q0 = initial_condition_specification(nComp, RadialConvDispOpInstance, bindStride, c0, cp0, q0)

		# specifying the integration method for the DGSEM
		exactInt = RadialConvDispOperatorDG.exact_integration()
		
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
		new(nComp, col_inner_radius, col_outer_radius, col_height, cross_section_area, d_rad, eps_c, eps_p, kf, Rp, Rc, Dp, c0, cp0, q0, cIn, exactInt, polyDeg, nCells, polyDegPore, RadialConvDispOpInstance, PoreOpInstance, bindStride, adsStride, unitStride, idx, idx_p, idx_q, Fc, Fp, Fjac, Jr, invRi, cpp, RHS_q, qq, film_vec, solution_outlet, solution_times, bind)
	end
end



# Define a function to compute the transport term for the GRM
function compute_transport!(RHS, RHS_q, cpp, x, m::rGRM, t, section, sink, switches, idx_units)

	# Determining inlet velocity if specified dynamically
	get_inlet_flows!(switches, switches.ConnectionInstance.dynamic_flow[switches.switchSetup[section], sink], section, sink, t, m)
	
	# Precompute face velocities once per call (same for all components)
	faces_v = RadialConvDispOperatorDG.compute_faces_v(switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink],m.RadialConvDispOpInstance.rho_i,m.RadialConvDispOpInstance.rho_ip1,m.col_inner_radius,m.RadialConvDispOpInstance.nCells)

	# Loop over components where convection dispersion term is determined and the isotherm term is subtracted
	# Convection Dispersion term	
	cpp_block = @view x[1 + idx_units[sink] : idx_units[sink] + m.RadialConvDispOpInstance.nPoints * m.nComp]
	@inbounds for j = 1:m.nComp
		# local “in-block” indices for component j
		idx_bulk = (1 + (j-1)*m.RadialConvDispOpInstance.nPoints) : (j*m.RadialConvDispOpInstance.nPoints)
		idx_bulk_abs = idx_bulk .+ idx_units[sink]

		# inlet concentration (unchanged)
		inlet_concentrations!(m.cIn, switches, j, section, sink, x, t, switches.inlet_conditions[section, sink, j])
		# per-component diffusion scales
		D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(m.d_rad[j], m.RadialConvDispOpInstance.rho_i, m.RadialConvDispOpInstance.rho_ip1)
		left_scale_vec  = @. m.RadialConvDispOpInstance.rho_i  * D_left
		right_scale_vec = @. m.RadialConvDispOpInstance.rho_ip1 * D_right

		# transport residual for mobile phase of comp j
		RadialConvDispOperatorDG.radialresidualImpl!(
			m.RadialConvDispOpInstance.Dc,
			cpp_block, idx_bulk,
			m.RadialConvDispOpInstance.strideNode,
			m.RadialConvDispOpInstance.strideCell,
			m.RadialConvDispOpInstance.nNodes,
			m.RadialConvDispOpInstance.nCells,
			m.RadialConvDispOpInstance.deltarho,
			m.polyDeg,
			m.RadialConvDispOpInstance.invWeights,
			m.RadialConvDispOpInstance.nodes,
			m.RadialConvDispOpInstance.polyDerM,
			m.RadialConvDispOpInstance.invMM,
			m.RadialConvDispOpInstance.MM,
			m.RadialConvDispOpInstance.invMrhoM,
			m.RadialConvDispOpInstance.SgMatrix,
			switches.ConnectionInstance.u_tot[switches.switchSetup[section], sink],
			m.d_rad[j],
			m.cIn[j],
			m.RadialConvDispOpInstance.c_star,
			m.RadialConvDispOpInstance.g_star,
			m.RadialConvDispOpInstance.Dg,
			m.RadialConvDispOpInstance.h,
			m.RadialConvDispOpInstance.mul1,
			m.exactInt,
			m.RadialConvDispOpInstance.rho_i,
			m.RadialConvDispOpInstance.rho_ip1,
			faces_v, left_scale_vec,
			right_scale_vec
			)

			# Film exchange uses the pore concentration at the particle surface (outer pore node)
			# Build cp_surf indices for the j-th component at each bulk radial point
			let Np = m.PoreOpInstance.nNodesPore
				# base offset to the beginning of j-th component pore block (absolute indices)
				base_pore_j = idx_units[sink] + m.nComp * m.RadialConvDispOpInstance.nPoints + (j-1) * (Np * m.RadialConvDispOpInstance.nPoints)
				@inbounds for i in 1:m.RadialConvDispOpInstance.nPoints
					# bulk index for i-th radial point (absolute)
					idx_c_i = idx_bulk_abs[i]
					# pore surface node is the last node within the pore slice of size Np
					idx_cp_surf_i = base_pore_j + (i-1) * Np + Np
					m.film_vec[i] = m.kf[j] * (x[idx_c_i] - x[idx_cp_surf_i])
				end
			end
		@. @views RHS[idx_bulk_abs] = m.RadialConvDispOpInstance.Dc - m.Fc * 3 / m.Rp * m.film_vec

		# Pore phase - Each idx has _nPolyPore number of states Adding the boundary flux
		# The boundary flux term is computed as kf(x-x*) 
		# The whole term is L (M^-1) (2/deltarho) Rp^2/eps_p * film_vec
		@. @views m.PoreOpInstance.boundaryPore = m.film_vec

		#Now the rest is computed
		#dcp/dt = 4/Rp Dp Dr/ri cp - (2/Rp)^2 Dp M^-1 A cp + 2/Rp Dp L[0,J/eps/Dp]  - Fp dq/dt
		@inbounds for i in 1:m.RadialConvDispOpInstance.nPoints
			#Pore phase for each component starts after all mobile phases
			#through all mobile phase, through pore phase j, at porephase i
			m.idx = 1 + m.nComp*m.RadialConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.RadialConvDispOpInstance.nPoints + (i-1) * m.PoreOpInstance.nNodesPore + idx_units[sink] : m.PoreOpInstance.nNodesPore + m.nComp*m.RadialConvDispOpInstance.nPoints + (j-1) *m.PoreOpInstance.nNodesPore*m.RadialConvDispOpInstance.nPoints  + (i-1) * m.PoreOpInstance.nNodesPore + idx_units[sink]
			m.idx_q = 1 + (j-1) *m.PoreOpInstance.nNodesPore*m.RadialConvDispOpInstance.nPoints + (i-1) * m.PoreOpInstance.nNodesPore : m.PoreOpInstance.nNodesPore + (j-1) *m.PoreOpInstance.nNodesPore*m.RadialConvDispOpInstance.nPoints  + (i-1) * m.PoreOpInstance.nNodesPore

			# Term 1 is 
			# (2/deltarho)^2 .* M^-1 A^r Dp cp
			mul!(m.PoreOpInstance.term1, m.PoreOpInstance.invMM_Ar,@view(x[m.idx])) 
			broadcast!(*, m.PoreOpInstance.term1, m.Dp[j], m.PoreOpInstance.term1) # 

			#Pore phase boundary conditions 
			#Corresponds to the boundary term - L (M^-1) (2/deltarho) Rp^2/eps_p * kf (x-x*)
			@. @views RHS[m.idx] = m.PoreOpInstance.boundaryPore[i] * m.PoreOpInstance.LiftMatrixRed 

			#Assembling RHS
			@. @views RHS[m.idx] += - m.PoreOpInstance.term1 - m.Fp * RHS_q[m.idx_q]
		end
	end
	
	nothing
end