module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    export set_debug_log, close_debug_log

    # Global log file handle
    const DEBUG_LOG = Ref{Union{IO,Nothing}}(nothing)

    function set_debug_log(filename::String)
        if DEBUG_LOG[] !== nothing
            close(DEBUG_LOG[])
        end
        DEBUG_LOG[] = open(filename, "w")
        println("Debug log file created: ", filename)
    end

    function close_debug_log()
        if DEBUG_LOG[] !== nothing
            close(DEBUG_LOG[])
            DEBUG_LOG[] = nothing
        end
    end

    # Helper function to print to both console and log file
    function debug_println(args...)
        println(args...)
        if DEBUG_LOG[] !== nothing
            println(DEBUG_LOG[], args...)
            flush(DEBUG_LOG[])
        end
    end

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _polyDerM, _invMM, _MM, _rMM, _invrMM, _nodes, v, d_rad, rho_i, cIn, c_star, g_star, Dg, _h, mul1)
        fill!(Dg, 0.0)   # reset auxiliary buffer used to build g
        fill!(Dc, 0.0)   # reset residual accumulator for mobile phase

		# diffusive face scales for comp j
		scale  = rho_i .* d_rad

        map = 2.0 / _deltarho

        # Build strong derivative in ξ: Dg ← -D * c
        volumeIntegraly!(y, idx, Dg, _nCells, _nNodes, _polyDerM, mul1)

        # Numerical fluxes c*
        interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells)

        # Lift face terms into g and scale to physical derivative
        surfaceIntegraly!(Dg, y, idx, _strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg)

        # g = (2/Δρ) * [ M^{-1}B(c* - c) - D c ]
        @. _h = map * Dg

        # Auxiliary fluxes g*
        interfaceFluxCentral!(g_star, _h, _strideNode, _strideCell, _nCells)

        # ========== DEBUGGING OUTPUT ==========
        debug_println("\n" * "="^80)
        debug_println("RADIAL RESIDUAL DEBUG - AUXILIARY PROBLEM")
        debug_println("="^80)
        debug_println("Input parameters:")
        debug_println("  v_in = ", v, " m/s")
        debug_println("  d_rad = ", d_rad, " m²/s")
        debug_println("  cIn = ", cIn)
        debug_println("  deltarho = ", _deltarho, " m")
        debug_println("  nCells = ", _nCells, ", nNodes = ", _nNodes)
        debug_println("\nRadial face positions (rho_i):")
        for i = 1:length(rho_i)
            debug_println("  Face ", i-1, ": rho_i[", i, "] = ", rho_i[i], " m")
        end
        debug_println("\nAuxiliary variable g at each node (after all operations):")
        for Cell = 1:_nCells
            debug_println("  Cell $Cell:")
            for node = 1:_nNodes
                idx_g = (Cell-1)*_nNodes + node
                debug_println("    Node ", node, ": _h[", idx_g, "] = ", _h[idx_g])
            end
        end
        debug_println("\nAuxiliary interface fluxes g* (central):")
        for face = 1:(_nCells+1)
            debug_println("  Face ", face-1, ": g_star[", face, "] = ", g_star[face])
        end

        # Numerical fluxes v*c* (convective conservative flux, upwind)
        #interfaceFluxUpwind!(c_star, y, idx, _strideNode, _strideCell, _nCells, v, rho_i, cIn)

        debug_println("\n" * "="^80)
        debug_println("RADIAL RESIDUAL DEBUG - MAIN PROBLEM")
        debug_println("="^80)
        debug_println("Concentration at each node:")
        base = first(idx)
        for Cell = 1:_nCells
            debug_println("  Cell $Cell:")
            for node = 1:_nNodes
                idx_c = base + (Cell-1)*_nNodes + (node-1)
                debug_println("    Node ", node, ": y[", idx_c, "] = ", y[idx_c])
            end
        end
        debug_println("\nConvective interface fluxes c*:")
        for face = 1:(_nCells+1)
            debug_println("  Face ", face-1, ": c_star[", face, "] = ", c_star[face])
        end
        debug_println("\nDiffusive face scales (rho_i * d_rad):")
        for face = 1:(_nCells+1)
            debug_println("  Face ", face-1, ": scale[", face, "] = ", scale[face])
        end

        # Surface Integral: B(v c*) - B_g(scale * g*)
        surfaceIntegral!(Dc, _nNodes, _nCells, idx, _strideNode, _strideCell, c_star, g_star, scale, v, rho_i)

        debug_println("\nAfter SURFACE integrals Dc:")
        for Cell = 1:_nCells
            debug_println("  Cell $Cell:")
            for node = 1:_nNodes
                idx_dc = (Cell-1)*_nNodes + node
                debug_println("    Node ", node, ": Dc[", idx_dc, "] = ", Dc[idx_dc])
            end
        end

        # Volume Integral: Sᵀ (v c) - S_g g, then apply M_ρ^{-1} and (2/Δρ) scaling
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _polyDerM, _MM, _rMM, _invrMM, _h, mul1, rho_i, _deltarho, d_rad, _nodes, map, v)

        debug_println("\nAfter VOLUME integrals (final Dc):")
        for Cell = 1:_nCells
            debug_println("  Cell $Cell:")
            for node = 1:_nNodes
                idx_dc = (Cell-1)*_nNodes + node
                debug_println("    Node ", node, ": Dc[", idx_dc, "] = ", Dc[idx_dc])
            end
        end
        debug_println("="^80)
        debug_println()

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

    # Compute strong derivative: Dᵀ M (v*c) - Dᵀ M_ρ (d_rad * g), then apply M_ρ^{-1}
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _MM::Matrix{Float64}, _rMM::Vector{Matrix{Float64}}, _invrMM::Vector{Matrix{Float64}}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Float64, _nodes::Vector{Float64}, map::Float64, v_in::Float64)
        base = first(idx)
        debug_println("\n" * "-"^80)
        debug_println("VOLUME INTEGRAL DEBUG")
        debug_println("-"^80)

        # Print matrices once (they're the same for all cells, except rMM and invrMM)
        debug_println("\nLegendre-Gauss-Lobatto (LGL) Nodes in reference coordinates [-1, 1]:")
        for i = 1:_nNodes
            debug_println("  Node ", i, ": ξ = ", _nodes[i])
        end

        debug_println("\nDG Polynomial Derivative Matrix (polyDerM) [", _nNodes, "x", _nNodes, "]:")
        debug_println("  (This computes dφ/dξ where φ are the Lagrange basis functions)")
        for i = 1:_nNodes
            debug_println("  Row ", i, ": ", _polyDerM[i,:])
        end
        debug_println("  (This computes dφ/dξ where φ are the Lagrange basis functions)")
        for i = 1:_nNodes
            debug_println("  Row ", i, ": ", transpose(_polyDerM[i,:]))
        end

        debug_println("\nStandard Mass Matrix M (MM) [", _nNodes, "x", _nNodes, "]:")
        debug_println("  (This is M(0,0) for standard integration without weighting)")
        for i = 1:_nNodes
            debug_println("  Row ", i, ": ", _MM[i,:])
        end
        @inbounds for Cell in 1:_nCells
            c_cell  = @view(y[base + (Cell - 1) * _nNodes : base - 1 + Cell * _nNodes])
            g_cell  = @view(g[1 + (Cell - 1) * _nNodes : Cell * _nNodes])
            Dc_cell = @view(Dc[1 + (Cell - 1) * _nNodes : Cell * _nNodes])
            rho_nodes = @. rho_i[Cell] + (_deltarho/2) * (1.0 + _nodes)
            v_nodes = v_in * rho_i[1] ./ rho_nodes
            # rho * v = v_in * rho_in = constant
            # v = v_in * rho_in / rho
            # v * c directly (the physical flux)
            @. mul1 = v_nodes * c_cell 

            debug_println("\n" * "="^80)
            debug_println("Cell $Cell:")
            debug_println("="^80)

            debug_println("\nRadially-Weighted Mass Matrix M_ρ for Cell ", Cell, " [", _nNodes, "x", _nNodes, "]:")
            for i = 1:_nNodes
                debug_println("  Row ", i, ": ", _rMM[Cell][i,:])
            end

            debug_println("\nInverse M_ρ^{-1} for Cell ", Cell, " [", _nNodes, "x", _nNodes, "]:")
            for i = 1:_nNodes
                debug_println("  Row ", i, ": ", _invrMM[Cell][i,:])
            end

            debug_println("\nRadial node positions (rho_nodes):")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": rho = ", rho_nodes[i], " m")
            end
            debug_println("  Velocity at nodes (v_nodes = v_in * rho_i[1] / rho):")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": v = ", v_nodes[i], " m/s")
            end
            debug_println("  Concentration at nodes (c_cell):")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": c = ", c_cell[i])
            end
            debug_println("  Product v*c at nodes (mul1 = v_nodes * c_cell):")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": v*c = ", mul1[i])
            end

            # convective volume: Dᵀ M (v*c)
            debug_println("\n--- Convective Term Computation ---")
            temp_vec_conv = _MM * mul1
            debug_println("  Intermediate: M * (v*c) =")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": ", temp_vec_conv[i])
            end
            tmp_conv = transpose(_polyDerM) * temp_vec_conv
            debug_println("  Final: tmp_conv = Dᵀ * M * (v*c) =")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": tmp_conv = ", tmp_conv[i])
            end

            # dispersion volume: S_g g with S_g = Dᵀ * M_ρ * d_rad
            debug_println("\n--- Dispersion Term Computation ---")
            debug_println("  d_rad * g_cell =")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": ", d_rad * g_cell[i])
            end
            temp_vec_diff = _rMM[Cell] * (d_rad .* g_cell)
            debug_println("  Intermediate: M_ρ * (d_rad * g) =")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": ", temp_vec_diff[i])
            end
            tmp_diff = transpose(_polyDerM) * temp_vec_diff
            debug_println("  Final: tmp_diff = Dᵀ * [M_ρ * (d_rad * g)] =")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": tmp_diff = ", tmp_diff[i])
            end

            # Show comparison between convective and dispersion terms
            debug_println("\n--- Comparison: tmp_conv vs tmp_diff ---")
            debug_println("  For free-stream preservation (c=1 everywhere), these should be EQUAL!")
            for i = 1:_nNodes
                diff_val = tmp_conv[i] - tmp_diff[i]
                debug_println("    Node ", i, ": tmp_conv = ", tmp_conv[i],
                             ", tmp_diff = ", tmp_diff[i], ", difference = ", diff_val)
            end

            # Accumulate: Dc += (tmp_conv - tmp_diff)
            @. Dc_cell += (tmp_conv - tmp_diff)
            debug_println("\n--- After adding volume terms (Dc += tmp_conv - tmp_diff) ---")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": Dc = ", Dc_cell[i])
            end

            # Apply M_ρ^{-1} and scaling (2/Δρ) in one operation
            @views Dc_cell .= map .* (_invrMM[Cell] * Dc_cell)
            debug_println("  After applying M_ρ^{-1} and (2/Δρ) scaling:")
            for i = 1:_nNodes
                debug_println("    Node ", i, ": Dc = ", Dc_cell[i])
            end
        end
        debug_println("-"^80)

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

    # Exact Integration lifting of Surface - [B(v*c*) - B_g(scale*g*)]
    @inline function surfaceIntegral!(_surfaceFlux::Vector{Float64}, nNodes::Int, _nCells::Int, idx::UnitRange{Int}, _strideNode::Int, _strideCell::Int, c_star::Vector{Float64}, g_star::Vector{Float64}, scale::Vector{Float64}, v_in::Float64, rho_i::Vector{Float64})
        # accumulate: B(v * c*) - B_g(scale * g*)
        base = first(idx)
        @inbounds for Cell in 1:_nCells
                if (v_in * rho_i[1] / rho_i[Cell]) >= 0  # Check velocity direction
                    _surfaceFlux[(Cell - 1) * nNodes + 1] -= (v_in * rho_i[1] ./ rho_i[Cell]) * c_star[base + (Cell - 1) * _strideCell - _strideNode] - (scale[Cell] * g_star[Cell])
                    _surfaceFlux[Cell * nNodes] -= (v_in * rho_i[1] ./ rho_i[Cell + 1]) * c_star[base + (Cell - 1) * _strideCell] - (scale[Cell + 1] * g_star[Cell + 1])
                end
        end
        return nothing
    end

    # Auxiliary fluxes c* = 0.5 (c_L + c_R)
    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C, idx::UnitRange{Int}, _strideNode::Int, _strideCell::Int, _nCells::Int)
        base = first(idx)
        @inbounds for Cell in 2:_nCells
                _surfaceFlux[Cell] = 0.5 * (C[base + (Cell - 1) * _strideCell - _strideNode] + C[base + (Cell - 1) * _strideCell])
            end
            _surfaceFlux[1] = C[base]
            _surfaceFlux[_nCells + 1] = C[base + _nCells * _strideCell - _strideNode]
        return nothing
    end

    # Numerical fluxes (c*)
    @inline function interfaceFluxUpwind!(_surfaceFlux::Vector{Float64}, C, idx::UnitRange{Int}, _strideNode::Int, _strideCell::Int, _nCells::Int, cIn::Float64, cOut::Union{Nothing, Float64}=nothing)
        base = first(idx)
        @inbounds for Cell in 2:_nCells
            if (v_in * rho_i[1] / rho_i[Cell]) >= 0
                _surfaceFlux[Cell] = (v_in * rho_i[1] ./ rho_[Cell]) * C[base + (Cell - 1) * _strideCell - _strideNode]
            else
                _surfaceFlux[Cell] = (v_in * rho_i[1] ./ rho_[Cell + 1]) * C[base + (Cell - 1) * _strideCell]
            end
        end
        _surfaceFlux[1] = v_in >= 0 ? (v_in) * cIn : (v_in * rho_i[1] ./ rho_[Cell]) * C[base]
        v_out = v_in * rho_i[1] / rho_i[_nCells + 1]
        if v_out >= 0
            _surfaceFlux[_nCells + 1] = (v_in * rho_i[1]) * C[base + _nCells * _strideCell - _strideNode]
        else
            if cOut === nothing
                _surfaceFlux[_nCells + 1] = (v_in * rho_i[1]) * C[base + _nCells * _strideCell - _strideNode]
            else
                _surfaceFlux[_nCells + 1] = (v_in * rho_i[1]) * cOut
            end
        end
        return nothing
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