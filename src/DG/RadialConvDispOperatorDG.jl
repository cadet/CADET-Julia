module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _polyDerM, _invMM, _MM01, _MM00, _rMM, _invrMM, _nodes, _weights, v, d_rad, rho_i, cIn, c_star, g_star, Dg, _h, mul1)
        fill!(Dg, 0.0)   # reset auxiliary buffer used to build g
        fill!(Dc, 0.0)   # reset residual accumulator for mobile phase
        map = 2.0 / _deltarho

        # Build strong derivative in ξ: Dg ← D * c
        volumeIntegraly!(y, idx, Dg, _nCells, _nNodes, _polyDerM, mul1)

        # Numerical fluxes c*
        interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells)

        # Lift face terms into g
        surfaceIntegraly!(Dg, y, idx, _strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg)

        # g = (2/Δρ) * [ M^{-1} B(c* - c) + D c ]
        @. _h = map * Dg

        # Auxiliary fluxes g*
        interfaceFlux!(g_star, _h, _strideNode, _strideCell, _nCells)

        # Volume Integral
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _polyDerM, _MM00, _rMM, _invrMM, _h, mul1, rho_i, _deltarho, d_rad, _nodes, _weights, v)

        # Surface Integral
        surfaceIntegral!(Dc, _nNodes, _nCells, c_star, g_star, d_rad, v, rho_i, _invrMM, cIn)

        @. Dc = map * Dc

        return nothing
    end

    # compute strong derivative: D c
    @inline function volumeIntegraly!(state, idx, stateDer::Vector{Float64}, _nCells::Int64, _nNodes::Int64, _polyDerM::Matrix{Float64}, mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[idx[1] + (Cell-1) * _nNodes : idx[1] - 1 + _nNodes + (Cell-1) * _nNodes]))
            broadcast!(+, @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]),@view(stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)]), mul1)
        end
        nothing
    end

    # Compute strong derivative
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _MM00::Matrix{Float64}, _rMM::Vector{Matrix{Float64}}, _invrMM::Vector{Matrix{Float64}}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Union{Float64, Function}, _nodes::Vector{Float64}, _weights::Vector{Float64}, v::Float64)
        base = first(idx)
        @inbounds for Cell in 1:_nCells
            broadcast!(+, @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), mul!(mul1, _invrMM[Cell], mul!(mul1, transpose(_polyDerM), mul!(mul1, _MM00, broadcast!(*, mul1, Ref(v), @view(y[base + (Cell - 1) * _nNodes : base - 1 + Cell * _nNodes]))))))
        end
        @inbounds for Cell in 1:_nCells
            broadcast!(-, @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), mul!(mul1, _invrMM[Cell], mul!(mul1, transpose(_polyDerM), mul!(mul1, _rMM[Cell], broadcast!(*, mul1, Ref(d_rad), @view(g[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]))))))
        end

        return nothing
    end

    # Exact-integration lifting of Surface contributions
    @inline function surfaceIntegraly!(stateDer,state,idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg)
        for Cell in 1:_nCells
            @inbounds @simd for Node in 1:_nNodes
                stateDer[1 + ( (Cell - 1) * strideCell_stateDer) + ( (Node - 1) * strideNode_stateDer)] +=
                    (_invMM[Node, 1])  * ( (state[idx[1] + ((Cell - 1) * strideCell_state)]) - (_surfaceFlux[Cell]) ) -
                    (_invMM[Node, end]) * ( (state[idx[1] + ((Cell - 1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (_surfaceFlux[Cell + 1]) )
            end
        end
        nothing
    end

    # Exact Integration lifting of Surface M_ρ{-1} [B(v c*) - B_g(g*)]
    @inline function surfaceIntegral!(Dc::Vector{Float64}, nNodes::Int, _nCells::Int, c_star::Vector{Float64}, g_star::Vector{Float64}, d_rad::Union{Float64, Function}, v::Float64, rho_i::Vector{Float64}, _invrMM::Vector{Matrix{Float64}}, cIn::Float64)
        surface_contrib = zeros(nNodes)

        @inbounds begin
            fill!(surface_contrib, 0.0)
            surface_contrib[1] = -(v * cIn - rho_i[1] * d_rad * g_star[1])
            surface_contrib[nNodes] = v * c_star[2] - rho_i[2] * d_rad * g_star[2]
            broadcast!(-, @view(Dc[1:nNodes]), @view(Dc[1:nNodes]), mul!(surface_contrib, _invrMM[1], surface_contrib))
        end

        @inbounds for Cell in 2:_nCells-1
            fill!(surface_contrib, 0.0)
            surface_contrib[1] = -(v * c_star[Cell] - rho_i[Cell] * d_rad * g_star[Cell])
            surface_contrib[nNodes] = v * c_star[Cell + 1] - rho_i[Cell + 1] * d_rad * g_star[Cell + 1]
            broadcast!(-, @view(Dc[(Cell - 1) * nNodes + 1 : Cell * nNodes]), @view(Dc[(Cell - 1) * nNodes + 1 : Cell * nNodes]), mul!(surface_contrib, _invrMM[Cell], surface_contrib))
        end

        @inbounds begin
            fill!(surface_contrib, 0.0)
            surface_contrib[1] = -(v * c_star[_nCells] - rho_i[_nCells] * d_rad * g_star[_nCells])
            surface_contrib[nNodes] = v * c_star[_nCells + 1]
            broadcast!(-, @view(Dc[(_nCells - 1) * nNodes + 1 : _nCells * nNodes]), @view(Dc[(_nCells - 1) * nNodes + 1 : _nCells * nNodes]), mul!(surface_contrib, _invrMM[_nCells], surface_contrib))
        end
        return nothing
    end

    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C, idx, _strideNode::Int64, _strideCell::Int64,_nCells::Int64)
        @inbounds for Cell in 2:_nCells
            _surfaceFlux[Cell] = 0.5 * ((C[idx[1] + (Cell - 1) * _strideCell - _strideNode]) + (C[idx[1] + (Cell - 1) * _strideCell]))
        end
        _surfaceFlux[1] = 0.5 * ((C[idx[1] ]) + (C[idx[1] ]))
        _surfaceFlux[_nCells+1] = 0.5 * ((C[idx[1] + _nCells * _strideCell - _strideNode]) + (C[idx[1] + _nCells * _strideCell - _strideNode]))
        return nothing
    end

    @inline function interfaceFlux!(_surfaceFlux::Vector{Float64}, G, _strideNode::Int, _strideCell::Int, _nCells::Int)
        @inbounds @simd for Cell in 2:_nCells
            _surfaceFlux[Cell] = 0.5 * (G[(Cell - 1) * _strideCell - _strideNode + 1] + G[(Cell - 1) * _strideCell + 1])
        end
        _surfaceFlux[1] = G[1]
        _surfaceFlux[_nCells + 1] = G[_nCells * _strideCell - _strideNode + 1]
        return nothing
    end

end
module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _polyDerM, _invMM, _MM01, _MM00, _rMM, _invrMM, _nodes, _weights, v, d_rad, rho_i, cIn, c_star, h_star, Dg, _h, mul1)

        fill!(Dg, 0.0)   # reset auxiliary buffer used to build g
        fill!(Dc, 0.0)   # reset residual accumulator for mobile phase
        map = 2.0 / _deltarho

        # Build strong derivative in ξ: Dg ← D * c
        volumeIntegraly!(y, idx, Dg, _nCells, _nNodes, _polyDerM, mul1)

        # Numerical fluxes c*
        interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells)

        # Lift face terms into g
        surfaceIntegraly!(Dg, y, idx, _strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg)

        # g = (2/Δρ) * [ M^{-1} B(c* - c) + D c ]
        @. _h = map * Dg

        # Compute complete numerical flux: h* = v*c* - rho_i*D_rad*g*
        interfaceFlux!(h_star, y, idx, d_rad, _h, v, _nCells, _nCells * _nNodes, _deltarho, cIn, _strideNode, _strideCell, 1, _nNodes, rho_i)

        # Volume Integral
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _polyDerM, _MM00, _rMM, _invrMM, _h, mul1, rho_i, _deltarho, d_rad, _nodes, _weights, v)

        # Surface Integral: Apply M_ρ^{-1} * B * (v*c* - rho_i*D*g*)
        surfaceIntegral!(Dc, 1, _nNodes, h_star, _nCells, _nNodes, _invrMM)

        @. Dc = map * Dc

        return nothing
    end

    # compute strong derivative: D c
    @inline function volumeIntegraly!(state, idx, stateDer::Vector{Float64}, _nCells::Int64, _nNodes::Int64, _polyDerM::Matrix{Float64}, mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[idx[1] + (Cell-1) * _nNodes : idx[1] - 1 + _nNodes + (Cell-1) * _nNodes]))
            broadcast!(+, @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]),@view(stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)]), mul1)
        end
        nothing
    end

    # Compute strong derivative
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _MM00::Matrix{Float64}, _rMM::Vector{Matrix{Float64}}, _invrMM::Vector{Matrix{Float64}}, S_g::Vector{Matrix{Float64}}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Union{Float64, Function}, _nodes::Vector{Float64}, _weights::Vector{Float64}, v::Float64)
        base = first(idx)
        @inbounds for Cell in 1:_nCells
            idx_range = (Cell - 1) * _nNodes + 1 : Cell * _nNodes
            broadcast!(+, @view(Dc[idx_range]), @view(Dc[idx_range]), mul!(mul1, _invrMM[Cell], mul!(mul1, transpose(_polyDerM), mul!(mul1, _MM00, broadcast!(*, mul1, Ref(v), @view(y[base + idx_range[1] - 1 : base + idx_range[end] - 1]))))))
            broadcast!(-, @view(Dc[idx_range]), @view(Dc[idx_range]), mul!(mul1, _invrMM[Cell], mul!(mul1, S_g[Cell], @view(g[idx_range]))))
        end
        nothing
    end

    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _MM00::Matrix{Float64}, _rMM::Vector{Matrix{Float64}}, _invrMM::Vector{Matrix{Float64}}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Union{Float64, Function}, _nodes::Vector{Float64}, _weights::Vector{Float64}, v::Float64)
        base = first(idx)
        @inbounds for Cell in 1:_nCells
            broadcast!(+, @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), mul!(mul1, _invrMM[Cell], mul!(mul1, transpose(_polyDerM), mul!(mul1, _MM00, broadcast!(*, mul1, Ref(v), @view(y[base + (Cell - 1) * _nNodes : base - 1 + Cell * _nNodes]))))))
        end
        @inbounds for Cell in 1:_nCells
            broadcast!(-, @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), mul!(mul1, _invrMM[Cell], mul!(mul1, transpose(_polyDerM), mul!(mul1, _rMM[Cell], broadcast!(*, mul1, Ref(d_rad), @view(g[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]))))))
        end

        return nothing
    end

    # Exact-integration lifting of Surface contributions
    @inline function surfaceIntegraly!(stateDer,state,idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg)
        for Cell in 1:_nCells
            @inbounds @simd for Node in 1:_nNodes
                stateDer[1 + ( (Cell - 1) * strideCell_stateDer) + ( (Node - 1) * strideNode_stateDer)] +=
                    (_invMM[Node, 1])  * ( (state[idx[1] + ((Cell - 1) * strideCell_state)]) - (_surfaceFlux[Cell]) ) -
                    (_invMM[Node, end]) * ( (state[idx[1] + ((Cell - 1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (_surfaceFlux[Cell + 1]) )
            end
        end
        nothing
    end
"""
    # Exact Integration lifting of Surface M_ρ{-1} [B(v c_internal - v c*) - B_g(0 - g*)]
    # Compute: M_ρ^{-1} * B * [(v * c_internal - v * c_numerical) - (0 - g*)]
    @inline function surfaceIntegral!(Dc::Vector{Float64}, y, idx, _strideNode::Int, _strideCell::Int, _polyDeg::Int, nNodes::Int, _nCells::Int, c_star::Vector{Float64}, g_star::Vector{Float64}, _invrMM::Vector{Matrix{Float64}}, v::Float64)
        base = first(idx)
        @inbounds for Cell in 1:_nCells
            # Get internal trace values at left and right faces of this cell
            c_left_internal = y[base + (Cell - 1) * _strideCell]  # first node of cell
            c_right_internal = y[base + (Cell - 1) * _strideCell + _polyDeg * _strideNode]  # last node of cell

            @simd for Node in 1:nNodes
                Dc[(Cell - 1) * nNodes + Node] +=
                    (_invrMM[Cell][Node, 1])   * (v * c_left_internal - c_star[Cell] + g_star[Cell]) -
                    (_invrMM[Cell][Node, end]) * (v * c_right_internal - c_star[Cell + 1] + g_star[Cell + 1])
            end
        end
        return nothing
    end
"""

    # Surface integral: Dc -= M_ρ^{-1} * B * flux = M_ρ^{-1} * (B_left * flux_left - B_right * flux_right)
    # where B_left = -1 at first node, B_right = +1 at last node
    @inline function surfaceIntegral!(stateDer, strideNode_stateDer, strideCell_stateDer, _surfaceFlux, _nCells, _nNodes, _invrMM)
        @inbounds for Cell in 1:_nCells
            @inbounds for Node in 1:_nNodes
                stateDer[1 + ((Cell-1) * strideCell_stateDer) + ((Node-1) * strideNode_stateDer)] -=
                    (_invrMM[Cell][Node, 1]) * (-_surfaceFlux[Cell]) +
                    (_invrMM[Cell][Node, end]) * (_surfaceFlux[Cell+1])
            end
        end
        return nothing
    end

    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C, idx, _strideNode::Int64, _strideCell::Int64,_nCells::Int64)
        @inbounds for Cell in 2:_nCells
            _surfaceFlux[Cell] = 0.5 * ((C[idx[1] + (Cell - 1) * _strideCell - _strideNode]) + (C[idx[1] + (Cell - 1) * _strideCell]))
        end
        _surfaceFlux[1] = 0.5 * ((C[idx[1] ]) + (C[idx[1] ]))
        _surfaceFlux[_nCells+1] = 0.5 * ((C[idx[1] + _nCells * _strideCell - _strideNode]) + (C[idx[1] + _nCells * _strideCell - _strideNode]))
        return nothing
    end

    # Compute numerical flux h* = v*c* - B_g*g*
    @inline function interfaceFlux!(_surfaceFlux, C, idx, d_rad, G, v, _nCells, _nPoints, _deltarho, cIn, _strideNode, _strideCell, _strideNode_g, _strideCell_g, rho_i)
        @inbounds for Cell in 2:_nCells
            _surfaceFlux[Cell] = v * 0.5 * ((C[idx[1] + (Cell-1) * _strideCell - _strideNode]) + (C[idx[1] + (Cell-1) * _strideCell])) - rho_i[Cell] * d_rad * 0.5 * ((G[1 + (Cell-1) * _strideCell_g - _strideNode_g]) + (G[1 + (Cell-1) * _strideCell_g]))
        end

        # Inlet boundary
        _surfaceFlux[1] = v * cIn

        # Outlet boundary
        _surfaceFlux[_nCells+1] = v * C[idx[1] + _nCells * _strideCell - _strideNode] - rho_i[_nCells+1] * d_rad * 0.5 * ((G[1 + _nCells * _strideCell_g - _strideNode_g]) + (-(G[_nPoints])))

        return nothing
    end

end