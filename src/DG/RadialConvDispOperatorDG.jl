module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    export enable_debug, disable_debug

    # Global debug flag
    const DEBUG_RADIAL_RESIDUAL = Ref(false)
    const DEBUG_LOG_FILE = Ref{Union{Nothing, IOStream}}(nothing)

    function enable_debug(logfile::String="radial_residual_debug.log")
        DEBUG_RADIAL_RESIDUAL[] = true
        DEBUG_LOG_FILE[] = open(logfile, "w")
        println(DEBUG_LOG_FILE[], "="^80)
        println(DEBUG_LOG_FILE[], "RADIAL RESIDUAL DEBUG LOG")
        println(DEBUG_LOG_FILE[], "="^80)
        flush(DEBUG_LOG_FILE[])
        println("Debug logging enabled. Writing to: $logfile")
    end

    function disable_debug()
        DEBUG_RADIAL_RESIDUAL[] = false
        if DEBUG_LOG_FILE[] !== nothing
            close(DEBUG_LOG_FILE[])
            DEBUG_LOG_FILE[] = nothing
        end
        println("Debug logging disabled.")
    end

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _polyDerM, _invMM, _MM01, _S, _rMM, _invrMM, _nodes, _weights, v, d_rad, rho_i, cIn, c_star, g_star, Dg, _h, mul1)
        fill!(Dg, 0.0)   # reset auxiliary buffer used to build g
        fill!(Dc, 0.0)   # reset residual accumulator for mobile phase
        scale  = rho_i .* d_rad
        map = 2.0 / _deltarho

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n" * "="^80)
            println(io, "RESIDUAL FUNCTION CALL - DETAILED DEBUG OUTPUT")
            println(io, "="^80)

            # 1. Input Parameters
            println(io, "\n--- INPUT PARAMETERS ---")
            println(io, "  nNodes  = $_nNodes")
            println(io, "  nCells  = $_nCells")
            println(io, "  polyDeg = $_polyDeg")
            println(io, "  deltarho = $_deltarho")
            println(io, "  v (inlet velocity) = $v")
            println(io, "  d_rad (dispersion) = $d_rad")
            println(io, "  cIn (inlet concentration) = $cIn")
            println(io, "  map = 2.0/deltarho = $map")

            # 2. Initial State
            println(io, "\n--- INITIAL STATE VECTOR y ---")
            println(io, "  idx range: $idx")
            println(io, "  State values by cell:")
            for cell in 1:_nCells
                cell_vals = y[idx[1] + (cell-1)*_nNodes : idx[1] + cell*_nNodes - 1]
                println(io, "    Cell $cell: $cell_vals")
            end

            # 3. Radial positions
            println(io, "\n--- RADIAL POSITIONS (rho_i) ---")
            println(io, "  Cell interfaces:")
            for i in eachindex(rho_i)
                println(io, "    rho_i[$i] = ", rho_i[i])
            end

            # 4. Nodes and Weights
            println(io, "\n--- QUADRATURE NODES AND WEIGHTS ---")
            println(io, "  Reference nodes (ξ ∈ [-1, 1]):")
            for i in eachindex(_nodes)
                println(io, "    node[$i] = ", _nodes[i], ",  weight[$i] = ", _weights[i])
            end

            # 5. Scale vector
            println(io, "\n--- SCALE VECTOR (rho_i .* d_rad) ---")
            println(io, "  scale = $scale")

            # 6. Matrices
            println(io, "\n--- DIFFERENTIATION MATRIX (_polyDerM) ---")
            println(io, "  Size: $(size(_polyDerM))")
            for i in axes(_polyDerM, 1)
                println(io, "  Row $i: ", _polyDerM[i,:])
            end

            _polyDerMTranspose = transpose(_polyDerM)
            println(io, "\n--- transpose DIFFERENTIATION MATRIX (_polyDerMTranspose) ---")
            println(io, "  Size: $(size(_polyDerM))")
            for i in axes(_polyDerMTranspose, 1)
                println(io, "  Row $i: ", _polyDerMTranspose[i,:])
            end

            _MM = inv(_invMM)
            println(io, "\n--- MASS MATRIX (_MM) ---")
            println(io, "  Size: $(size(_MM))")
            for i in axes(_MM, 1)
                println(io, "  Row $i: ", _MM[i,:])
            end

            println(io, "\n--- MASS MATRIX 0,1 (_MM01) ---")
            println(io, "  Size: $(size(_MM01))")
            for i in axes(_MM01, 1)
                println(io, "  Row $i: ", _MM01[i,:])
            end

            println(io, "\n--- INVERSE MASS MATRIX (_invMM) ---")
            println(io, "  Size: $(size(_invMM))")
            for i in axes(_invMM, 1)
                println(io, "  Row $i: ", _invMM[i,:])
            end

            X = _polyDerMTranspose * _MM
            println(io, "\n--- OPERATOR Dᵀ * MM ---")
            println(io, "  Size: $(size(X))")
            for i in axes(X, 1)
                println(io, "  Row $i: ", X[i,:])
            end

            println(io, "\n--- SKEW-SYMMETRIC MATRIX (_S) ---")
            println(io, "  Size: $(size(_S))")
            for i in axes(_S, 1)
                println(io, "  Row $i: ", _S[i,:])
            end

            println(io, "\n--- RADIAL MASS MATRICES (_rMM) ---")
            for cell in 1:_nCells
                println(io, "  Cell $cell:")
                for i in axes(_rMM[cell], 1)
                    println(io, "    Row $i: ", _rMM[cell][i,:])
                end
            end

            println(io, "\n--- INVERSE RADIAL MASS MATRICES (_invrMM) ---")
            for cell in 1:_nCells
                println(io, "  Cell $cell:")
                for i in axes(_invrMM[cell], 1)
                    println(io, "    Row $i: ", _invrMM[cell][i,:])
                end
            end

            flush(io)
        end

        # Build strong derivative in ξ: Dg ← -D * c
        volumeIntegraly!(y, idx, Dg, _nCells, _nNodes, _polyDerM, mul1)

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n--- STEP 1: VOLUME DERIVATIVE (Dg = -D*c) ---")
            println(io, "  After volumeIntegraly!:")
            for cell in 1:_nCells
                cell_vals = Dg[(cell-1)*_nNodes + 1 : cell*_nNodes]
                println(io, "    Cell $cell Dg: $cell_vals")
            end
            flush(io)
        end

        # Numerical fluxes c*
        interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells)

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n--- STEP 2: AUXILIARY NUMERICAL FLUXES (c_star) ---")
            println(io, "  Interface fluxes c_star (averaged concentrations):")
            for i in eachindex(c_star)
                println(io, "    c_star[$i] = ", c_star[i])
            end
            flush(io)
        end

        # Lift face terms into g and scale to physical derivative
        surfaceIntegraly!(Dg, y, idx, _strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg)

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n--- STEP 3: SURFACE LIFTING (Dg after surfaceIntegraly!) ---")
            println(io, "  After adding surface contributions:")
            for cell in 1:_nCells
                cell_vals = Dg[(cell-1)*_nNodes + 1 : cell*_nNodes]
                println(io, "    Cell $cell Dg: $cell_vals")
            end
            flush(io)
        end

        # g = (2/Δρ) * [ M^{-1}B(c* - c) - D c ]
        @. _h = map * Dg

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n--- STEP 4: SCALE TO PHYSICAL DERIVATIVE (_h = map * Dg) ---")
            println(io, "  map = $map")
            println(io, "  _h (auxiliary variable g):")
            for cell in 1:_nCells
                cell_vals = _h[(cell-1)*_nNodes + 1 : cell*_nNodes]
                println(io, "    Cell $cell _h: $cell_vals")
            end
            flush(io)
        end

        # Auxiliary fluxes g*
        interfaceFluxCentral!(g_star, _h, _strideNode, _strideCell, _nCells)

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n--- STEP 5: AUXILIARY FLUXES g_star ---")
            println(io, "  Interface fluxes g_star:")
            for i in eachindex(g_star)
                println(io, "    g_star[$i] = ", g_star[i])
            end
            flush(io)
        end

        # Volume Integral: M_ρ^{-1} ( Sᵀ (v c) - S_g g ),
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _polyDerM, _S, _rMM, _invrMM, _h, mul1, rho_i, _deltarho, d_rad, _nodes, _weights, v)

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n--- STEP 6: VOLUME INTEGRAL ---")
            println(io, "  Dc after volumeIntegral! [Sᵀ(v c) - Dᵀ M_ρ (d_rad * g), then M_ρ^{-1}]:")
            for cell in 1:_nCells
                cell_vals = Dc[(cell-1)*_nNodes + 1 : cell*_nNodes]
                println(io, "    Cell $cell Dc: $cell_vals")
            end
            flush(io)
        end

        # Surface Integral: M_ρ^{-1} ( B(v c*) - B_g(scale * g*) )
        surfaceIntegral!(Dc, _nNodes, _nCells, c_star, g_star, scale, v, rho_i, _invrMM)

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n--- STEP 7: SURFACE INTEGRAL ---")
            println(io, "  Dc after surfaceIntegral! [B(v c*) - B_g(scale * g*)]:")
            for cell in 1:_nCells
                cell_vals = Dc[(cell-1)*_nNodes + 1 : cell*_nNodes]
                println(io, "    Cell $cell Dc: $cell_vals")
            end
            flush(io)
        end

        @. Dc = map * Dc

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n--- STEP 8: FINAL SCALING (Dc = map * Dc) ---")
            println(io, "  map = $map")
            println(io, "  FINAL RESIDUAL Dc:")
            for cell in 1:_nCells
                cell_vals = Dc[(cell-1)*_nNodes + 1 : cell*_nNodes]
                println(io, "    Cell $cell Dc: $cell_vals")
            end

            println(io, "\n  Residual statistics:")
            println(io, "    Max |Dc| = ", maximum(abs.(Dc)))
            println(io, "    Min Dc   = ", minimum(Dc))
            println(io, "    Max Dc   = ", maximum(Dc))
            println(io, "    Mean Dc  = ", sum(Dc)/length(Dc))

            println(io, "\n" * "="^80)
            flush(io)
        end

        return nothing
    end

    # compute strong derivative: - D c
    @inline function volumeIntegraly!(state, idx, stateDer::Vector{Float64}, _nCells::Int64, _nNodes::Int64, _polyDerM::Matrix{Float64}, mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[idx[1] + (Cell-1) * _nNodes : idx[1] - 1 + _nNodes + (Cell-1) * _nNodes]))
            broadcast!(-, @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]),@view(stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)]), mul1)
        end
        nothing
    end

    # Compute strong derivative: M_ρ^{-1} (Sᵀ(v*c) - Dᵀ M_ρ (d_rad * g)),
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _S::Matrix{Float64}, _rMM::Vector{Matrix{Float64}}, _invrMM::Vector{Matrix{Float64}}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Union{Float64, Function}, _nodes::Vector{Float64}, _weights::Vector{Float64}, v_in::Float64)
        base = first(idx)

        if DEBUG_RADIAL_RESIDUAL[]
            io = DEBUG_LOG_FILE[]
            println(io, "\n  --- VOLUME INTEGRAL DETAILED COMPUTATION ---")
        end

        @inbounds for Cell in 1:_nCells
            c_cell  = @view(y[base + (Cell - 1) * _nNodes : base - 1 + Cell * _nNodes])
            g_cell  = @view(g[1 + (Cell - 1) * _nNodes : Cell * _nNodes])
            Dc_cell = @view(Dc[1 + (Cell - 1) * _nNodes : Cell * _nNodes])
            rho_nodes = @. rho_i[Cell] + (_deltarho/2) * (1.0 + _nodes)
            v_nodes = v_in * rho_i[1] ./ rho_nodes
            @. mul1 = v_nodes * c_cell

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "\n    Cell $Cell:")
                println(io, "      rho_left  = ", rho_i[Cell])
                println(io, "      rho_right = ", rho_i[Cell+1])
                println(io, "      deltarho  = ", _deltarho)

                println(io, "      rho at nodes:")
                for i in eachindex(rho_nodes)
                    println(io, "        rho_nodes[$i] = ", rho_nodes[i])
                end

                println(io, "      velocity at nodes (v = v_in * rho_in / rho):")
                for i in eachindex(v_nodes)
                    println(io, "        v_nodes[$i] = ", v_nodes[i])
                end

                println(io, "      c_cell = ", collect(c_cell))
                println(io, "      g_cell = ", collect(g_cell))
                println(io, "      mul1 (v*c) = $mul1")
            end

            # convective volume: Sᵀ (v*c)
            tmp_conv = transpose(_S) * mul1

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "      tmp_conv (Sᵀ * v*c) = $tmp_conv")
            end

            # dispersion volume: S_g g
            temp_vec_diff = _rMM[Cell] * (d_rad .* g_cell)
            tmp_diff = transpose(_polyDerM) * temp_vec_diff

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "      d_rad * g_cell = ", d_rad .* g_cell)
                println(io, "      M_ρ * (d_rad * g) = $temp_vec_diff")
                println(io, "      tmp_diff (Dᵀ * M_ρ * d_rad * g) = $tmp_diff")
                println(io, "      Dc_cell (before adding volume terms) = ", collect(Dc_cell))
            end

            volume_contrib = tmp_conv - tmp_diff

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "      volume_contrib (before M_ρ^{-1}) = $volume_contrib")
            end

            # Apply M_ρ^{-1}
            volume_contrib_invM = _invrMM[Cell] * volume_contrib

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "      volume_contrib (after M_ρ^{-1}) = $volume_contrib_invM")
            end

            @. Dc_cell += volume_contrib_invM

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "      Dc_cell (after adding volume contribution) = ", collect(Dc_cell))
                flush(io)
            end
        end

        return nothing
    end

    # Exact-integration lifting of Surface contributions
    @inline function surfaceIntegraly!(stateDer,state,idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg)
        for Cell in 1:_nCells
            @inbounds @simd for Node in 1:_nNodes
                stateDer[1 + ( (Cell - 1) * strideCell_stateDer) + ( (Node - 1) * strideNode_stateDer)] -=
                    (_invMM[Node, 1])  * ( (state[idx[1] + ((Cell - 1) * strideCell_state)]) - (_surfaceFlux[Cell]) ) -
                    (_invMM[Node, end]) * ( (state[idx[1] + ((Cell - 1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (_surfaceFlux[Cell + 1]) )
            end
        end
        nothing
    end

    # Exact Integration lifting of Surface - M_ρ{-1} [B(v c*) - B_g(scale g*)]
    @inline function surfaceIntegral!(Dc::Vector{Float64}, nNodes::Int, _nCells::Int, c_star::Vector{Float64}, g_star::Vector{Float64}, scale::Vector{Float64}, v_in::Float64, rho_i::Vector{Float64}, _invrMM::Vector{Matrix{Float64}})
        @inbounds for Cell in 1:_nCells
            Dc_cell = @view(Dc[(Cell - 1) * nNodes + 1 : Cell * nNodes])

            v_left = v_in * rho_i[1] / rho_i[Cell]
            v_right = v_in * rho_i[1] / rho_i[Cell + 1]

            # Compute surface contributions: B(v c*) - B_g(scale * g*)
            # Initialize surface
            surface_contrib = zeros(nNodes)

            # Left boundary
            surface_contrib[1] = - (v_left * c_star[Cell] - (scale[Cell] * g_star[Cell]))

            # Right boundary
            surface_contrib[nNodes] = - (v_right * c_star[Cell + 1] - (scale[Cell + 1] * g_star[Cell + 1]))

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "\n    Surface Cell $Cell:")
                println(io, "      v_left = $v_left, v_right = $v_right")
                println(io, "      c_star[left=$(Cell)] = $(c_star[Cell]), c_star[right=$(Cell+1)] = $(c_star[Cell+1])")
                println(io, "      g_star[left=$(Cell)] = $(g_star[Cell]), g_star[right=$(Cell+1)] = $(g_star[Cell+1])")
                println(io, "      scale[left=$(Cell)] = $(scale[Cell]), scale[right=$(Cell+1)] = $(scale[Cell+1])")
                println(io, "      surface_contrib (before M_ρ^{-1}) = $surface_contrib")
            end

            # Apply M_ρ^{-1}
            surface_contrib_invM = _invrMM[Cell] * surface_contrib

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "      surface_contrib (after M_ρ^{-1}) = $surface_contrib_invM")
            end

            @. Dc_cell += surface_contrib_invM

            if DEBUG_RADIAL_RESIDUAL[]
                io = DEBUG_LOG_FILE[]
                println(io, "      Dc_cell (after adding surface) = ", collect(Dc_cell))
                flush(io)
            end
        end

        return nothing
    end

    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C,idx, _strideNode::Int64, _strideCell::Int64,_nCells::Int64)
        # Auxiliary flux: c* = 0.5 (c_l + c_r) for g - Determines the interfaces (because of lifting matrix, B) between the cells - hence the length is nCells +1
        
        # _surfaceFlux = zeros(_nCells+1)
        # calculate inner interface fluxes
        @inbounds for Cell in 2:_nCells
            _surfaceFlux[Cell] = 0.5 * ((C[idx[1] + (Cell-1) * _strideCell - _strideNode]) + (C[idx[1] + (Cell-1) * _strideCell]))
        end    
        # calculate boundary interface fluxes
        _surfaceFlux[1] = 0.5 * ((C[idx[1] ]) + (C[idx[1] ]))  # left boundary interface
        _surfaceFlux[_nCells+1] = 0.5 * ((C[idx[1] + _nCells * _strideCell - _strideNode]) + (C[idx[1] + _nCells * _strideCell - _strideNode]))  # right boundary interface
        nothing
    end

    # Numerical fluxes (auxiliary central g*)
    @inline function interfaceFluxCentral!(_surfaceFlux::Vector{Float64}, G, _strideNode::Int, _strideCell::Int, _nCells::Int)
        @inbounds @simd for Cell in 2:_nCells
            _surfaceFlux[Cell] = 0.5 * (G[(Cell - 1) * _strideCell - _strideNode + 1] + G[(Cell - 1) * _strideCell + 1])
        end
        _surfaceFlux[1] = G[1]
        _surfaceFlux[_nCells + 1] = G[_nCells * _strideCell - _strideNode + 1]
        return nothing
    end

end