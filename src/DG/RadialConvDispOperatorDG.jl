module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _polyDerM, _invMM, _MM01, _MM00, _rMM, _invrMM, _nodes, _weights, v, d_rad, rho_i, cIn, c_star, g_star, Dg, _h, mul1)
        fill!(Dg, 0.0)   # reset auxiliary buffer used to build g
        fill!(Dc, 0.0)   # reset residual accumulator for mobile phase
        scale  = rho_i .* d_rad
        map = 2.0 / _deltarho

        # Build strong derivative in ξ: Dg ← D * c
        volumeIntegraly!(y, idx, Dg, _nCells, _nNodes, _polyDerM, mul1)

        # Numerical fluxes c*
        interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells, cIn)

        # Lift face terms into g and scale to physical derivative
        surfaceIntegraly!(Dg, y, idx, _strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg)

        # g = (2/Δρ) * [ M^{-1}B(c* - c) + D c ]
        @. _h = map * Dg

        # Auxiliary fluxes g*
        interfaceFluxCentral!(g_star, _h, _strideNode, _strideCell, _nCells)

        # Volume Integral: M_ρ^{-1} ( Dᵀ*M(0,0) (v c) - Dᵀ M_ρ (d_rad * g) ),
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _polyDerM, _MM00, _rMM, _invrMM, _h, mul1, rho_i, _deltarho, d_rad, _nodes, _weights, v)

        # Surface Integral: M_ρ^{-1} ( B(v c*) - B_g(scale * g*) )
        surfaceIntegral!(Dc, _nNodes, _nCells, c_star, g_star, scale, v, rho_i, _invrMM)

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

    # Compute strong derivative: M_ρ^{-1} (Dᵀ*M(0,0)*(v*c) - Dᵀ M_ρ (d_rad * g)),
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _MM00::Matrix{Float64}, _rMM::Vector{Matrix{Float64}}, _invrMM::Vector{Matrix{Float64}}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Union{Float64, Function}, _nodes::Vector{Float64}, _weights::Vector{Float64}, v::Float64)
        base = first(idx)
        # Convection term: M_ρ^{-1} * Dᵀ * M(0,0) * (v*c)
        @inbounds for Cell in 1:_nCells
            broadcast!(+, @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), mul!(mul1, _invrMM[Cell], mul!(mul1, transpose(_polyDerM), mul!(mul1, _MM00, broadcast!(*, mul1, Ref(v), @view(y[base + (Cell - 1) * _nNodes : base - 1 + Cell * _nNodes]))))))
        end
        # Diffusion term: -M_ρ^{-1} * Dᵀ * M_ρ * (d_rad * g)
        @inbounds for Cell in 1:_nCells
            broadcast!(-, @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), mul!(mul1, _invrMM[Cell], mul!(mul1, transpose(_polyDerM), mul!(mul1, _rMM[Cell], broadcast!(*, mul1, Ref(d_rad), @view(g[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]))))))
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
    @inline function surfaceIntegral!(Dc::Vector{Float64}, nNodes::Int, _nCells::Int, c_star::Vector{Float64}, g_star::Vector{Float64}, scale::Vector{Float64}, v::Float64, rho_i::Vector{Float64}, _invrMM::Vector{Matrix{Float64}})
        surface_contrib = zeros(nNodes)
        @inbounds for Cell in 1:_nCells
            fill!(surface_contrib, 0.0)
            surface_contrib[1] = -(v * c_star[Cell] - scale[Cell] * g_star[Cell])
            surface_contrib[nNodes] = v * c_star[Cell + 1] - scale[Cell + 1] * g_star[Cell + 1]
            broadcast!(+, @view(Dc[(Cell - 1) * nNodes + 1 : Cell * nNodes]), @view(Dc[(Cell - 1) * nNodes + 1 : Cell * nNodes]), mul!(surface_contrib, _invrMM[Cell], surface_contrib))
        end
        return nothing
    end

    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C,idx, _strideNode::Int64, _strideCell::Int64,_nCells::Int64, cIn::Float64)
        # Auxiliary flux: c* = 0.5 (c_l + c_r) for g - Determines the interfaces (because of lifting matrix, B) between the cells - hence the length is nCells +1

        # _surfaceFlux = zeros(_nCells+1)
        # calculate inner interface fluxes
        @inbounds for Cell in 2:_nCells
            _surfaceFlux[Cell] = 0.5 * ((C[idx[1] + (Cell-1) * _strideCell - _strideNode]) + (C[idx[1] + (Cell-1) * _strideCell]))
        end
        # calculate boundary interface fluxes
        _surfaceFlux[1] = cIn  # inlet boundary
        _surfaceFlux[_nCells+1] = C[idx[1] + _nCells * _strideCell - _strideNode]  # outlet boundary
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