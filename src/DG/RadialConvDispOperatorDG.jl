module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _polyDerM, _invMM, _S, _rMM, _invrMM, _nodes, _weights, v, d_rad, rho_i, cIn, c_star, g_star, Dg, _h, mul1)
        fill!(Dg, 0.0)   # reset auxiliary buffer used to build g
        fill!(Dc, 0.0)   # reset residual accumulator for mobile phase

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

        # Surface Integral: B(v c*) - B_g(scale * g*)
        surfaceIntegral!(Dc, _nNodes, _nCells, idx, _strideNode, _strideCell, c_star, g_star, scale, v, rho_i)

        # Volume Integral: Sᵀ (v c) - S_g g, then apply M_ρ^{-1} and (2/Δρ) scaling
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _polyDerM, _S, _rMM, _invrMM, _h, mul1, rho_i, _deltarho, d_rad, _nodes, _weights, map, v)

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
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _S::Matrix{Float64}, _rMM::Vector{Matrix{Float64}}, _invrMM::Vector{Matrix{Float64}}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Union{Float64, Function}, _nodes::Vector{Float64}, _weights::Vector{Float64}, map::Float64, v_in::Float64)
        base = first(idx)
        @inbounds for Cell in 1:_nCells
            c_cell  = @view(y[base + (Cell - 1) * _nNodes : base - 1 + Cell * _nNodes])
            g_cell  = @view(g[1 + (Cell - 1) * _nNodes : Cell * _nNodes])
            Dc_cell = @view(Dc[1 + (Cell - 1) * _nNodes : Cell * _nNodes])
            rho_nodes = @. rho_i[Cell] + (_deltarho/2) * (1.0 + _nodes)
            v_nodes = v_in * rho_i[1] ./ rho_nodes
            @. mul1 = v_nodes * c_cell

            # convective volume: Sᵀ (v*c)
            tmp_conv = transpose(_S) * mul1

            # dispersion volume: S_g g with S_g = Dᵀ * M_ρ * d_rad(r)
            temp_vec_diff = _rMM[Cell] * (d_rad .* g_cell)
            tmp_diff = transpose(_polyDerM) * temp_vec_diff

            # Accumulate: Dc += (tmp_conv - tmp_diff)
            @. Dc_cell += (tmp_conv - tmp_diff)

            # Apply M_ρ^{-1} and scaling (2/Δρ) in one operation
            @views Dc_cell .= map .* (_invrMM[Cell] * Dc_cell)
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