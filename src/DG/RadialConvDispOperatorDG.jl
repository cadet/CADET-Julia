module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _invWeights, _polyDerM, _invMM, _MM00, _MM01, _nodes, v, d_rad, rho_i, cIn, c_star, g_star, Dg, _h, mul1)
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
        interfaceFluxCentral!(g_star, _h, idx, _strideNode, _strideCell, _nCells)

        # Numerical fluxes v c*
        interfaceFluxUpwind!(c_star, y, idx, _strideNode, _strideCell, _nCells, v, rho_i, cIn)

        # Surface Integral: B (v c*) - B_g (g*)
        surfaceIntegral!(Dc, _nNodes, _nCells, c_star, g_star, scale)

        # Volume Intergral: Sᵀ (v c) - S_g g
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _polyDerM, _MM00, _MM01, _h, mul1, rho_i, _deltarho, d_rad, _nodes, v)        
        
        # °c = M_ρ^{-1} c
        @inbounds for Cell in 1:_nCells
            @views Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes] = ( rho_i[Cell] .* _MM00 .+ (_deltarho / 2) .* _MM01 ) \ Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]
        end

        # °c = (2/Δρ) * [...]
        @. Dc = map * Dc

        return nothing
    end

    # compute strong derivative: - D c
    @inline function volumeIntegraly!(state, idx, stateDer::Vector{Float64}, _nCells::Int64, _nNodes::Int64, _polyDerM::Matrix{Float64}, mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes])) 
            broadcast!(-, @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]),@view(stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)]), mul1)
        end
        nothing
    end

    # Compute strong derivative: Dᵀ M(0,0) (v c) - S_g g
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _MM00::Matrix{Float64}, _MM01::Matrix{Float64}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Float64, _nodes::Vector{Float64}, v_in::Float64)
        base = first(idx)
        @inbounds @simd for Cell in 1:_nCells
            c_cell  = @view(y[base + (Cell - 1) * _nNodes : base - 1 + Cell * _nNodes])
            g_cell  = @view(g[1 + (Cell - 1) * _nNodes : Cell * _nNodes])
            Dc_cell = @view(Dc[1 + (Cell - 1) * _nNodes : Cell * _nNodes])
            #rho_nodes = @. rho_i[Cell] + (_deltarho/2) * (1.0 + _nodes)
            @. mul1 = (v_in * rho_i[1] / rho_i[Cell]) * c_cell

            # convective volume: Sᵀ ( (v c) )
            tmp_conv = transpose(_polyDerM) * (_MM00 * mul1)
            # dispersion volume: S_g g with S_g = Dᵀ * M_ρ * d_rad
            tmp_diff = transpose(_polyDerM) * ( (rho_i[Cell] .* _MM00 .+ (_deltarho/2) .* _MM01) * (d_rad .* g_cell) )
            
            # Accumulate
            @. Dc_cell += (tmp_conv - tmp_diff)
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

    # Exact Integration lifting of Surface  - [B (v c*) - B_g (g*)]
    @inline function surfaceIntegral!(Surf::Vector{Float64}, nNodes::Int, nCells::Int, c_star::Vector{Float64}, g_star::Vector{Float64}, scale::Vector{Float64})
        @inbounds for Cell in 1:nCells
            # accumulate: B (v c*) - B_g (g*)
            Surf[(Cell - 1) * nNodes + 1] -= (- c_star[Cell] - (- (scale[Cell] * g_star[Cell])))
            Surf[Cell * nNodes] -= (c_star[Cell + 1] - (scale[Cell + 1] * g_star[Cell + 1]))
        end
        return nothing
    end

    # Auxiliary fluxes c* = 0.5 (c_L + c_R)
    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C, idx::UnitRange{Int}, strideNode::Int, strideCell::Int, nCells::Int)
        base = first(idx)
        @inbounds for Cell in 2:nCells
                _surfaceFlux[Cell] = 0.5 * (C[base + (Cell - 1) * strideCell - strideNode] + C[base + (Cell - 1) * strideCell])
            end
            _surfaceFlux[1] = C[base]
            _surfaceFlux[nCells + 1] = C[base + nCells * strideCell - strideNode]
        return nothing
    end

    # Numerical fluxes (upwinding c*)
    @inline function interfaceFluxUpwind!(_surfaceFlux::Vector{Float64}, C, idx::UnitRange{Int}, _strideNode::Int, _strideCell::Int, nCells::Int, v_in::Float64, rho_i::Vector{Float64}, cIn::Float64, cOut::Union{Nothing, Float64}=nothing)
        base = first(idx)
        @inbounds for Cell in 2:nCells
            if (v_in * rho_i[1] / rho_i[Cell]) >= 0
                _surfaceFlux[Cell] = (v_in * rho_i[1] / rho_i[Cell]) * C[base + (Cell - 1) * _strideCell - _strideNode]
            else
                _surfaceFlux[Cell] = (v_in * rho_i[1] / rho_i[Cell]) * C[base + (Cell - 1) * _strideCell]
            end
        end
        _surfaceFlux[1] = v_in >= 0 ? v_in * cIn : v_in * C[base]
        v_out = v_in * rho_i[1] / rho_i[nCells + 1]
        if v_out >= 0
            _surfaceFlux[nCells + 1] = v_out * C[base + nCells * _strideCell - _strideNode]
        else
            _surfaceFlux[nCells + 1] = isnan(cOut) ? v_out * C[base + nCells * strideCell - strideNode] : v_out * cOut
        end
        return nothing
    end

    # Numerical fluxes (auxiliary central g*)
    @inline function interfaceFluxCentral!(_surfaceFlux::Vector{Float64}, G, idx::UnitRange{Int}, strideNode::Int, strideCell::Int, nCells::Int)
        base = first(idx)
        @inbounds @simd for Cell in 2:nCells
                _surfaceFlux[Cell] = 0.5 * (G[base + (Cell - 1) * strideCell - strideNode] + G[base + (Cell - 1) * strideCell])
            end
            _surfaceFlux[1] = G[base]
            _surfaceFlux[nCells + 1] = G[base + nCells * strideCell - strideNode]
        return nothing
    end

end