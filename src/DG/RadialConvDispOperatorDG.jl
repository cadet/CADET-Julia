module RadialConvDispOperatorDG
    using LinearAlgebra

    # d_rad_i: Vector of dispersion coefficients at each cell interface (length nCells+1)
    #          For constant D, pass fill(D, nCells+1) or use the scalar version
    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _polyDerM, _invMM, _MM01, _MM00, _rMM, _invrMM, S_g, _nodes, _weights, v, d_rad_i, rho_i, cIn, c_star, g_star, Dg, g, mul1, mul2)
        fill!(Dg, 0.0)   # reset auxiliary buffer used to build g
        fill!(Dc, 0.0)   # reset residual accumulator for mobile phase
        map = 2.0 / _deltarho
        #v /= rho_i[1] / 4.0
        #v /= rho_i[1] / 3.0
        #v /= rho_i[1] / 2.0
        #v /= rho_i[1]
        #v *= map
        #v *= map / rho_i[1]
        v *= (rho_i[1] + rho_i[end]) / 2.0

        # Strong derivative in ξ: Dg ← D * c
        auxiliaryVolumeIntegral!(y, idx, _strideNode, _strideCell, Dg, _nCells, _nNodes, _polyDerM, mul1)

        # Numerical fluxes c* for auxiliary equation
        interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells)

        # Lift face terms into Dg: Dg += M^{-1} B (c* - c)
        auxiliarySurfaceIntegral!(Dg, y, idx, _strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg)

        # g = (2/Δρ) * Dg
        @. g = map * Dg

        computeNumericalFluxes!(c_star, g_star, y, idx, g, _strideNode, _strideCell, _nNodes, _nCells, v, d_rad_i, cIn)

        # Dc -= (2/Δρ) * M_ρ^{-1} * (B * v * c* - B_g * g*)
        surfaceIntegral!(Dc, c_star, g_star, _nCells, _nNodes, _invrMM, v, d_rad_i, rho_i, _deltarho)

        # Dc += (2/Δρ) * M_ρ^{-1} * (D^T * M^{(0,0)} * v * c - S_g * g)
        volumeIntegral!(Dc, y, idx, _strideNode, _strideCell, g, _nCells, _nNodes, _polyDerM, _MM00, _invrMM, S_g, _deltarho, v, rho_i, mul1, mul2)

        return nothing
    end

    # Auxiliary Central flux c*
    @inline function interfaceFluxAuxiliary!(c_star::Vector{Float64}, C, idx, _strideNode::Int64, _strideCell::Int64, _nCells::Int64)
        @inbounds for Cell in 2:_nCells
            c_star[Cell] = 0.5 * (C[idx[1] + (Cell - 1) * _strideCell - _strideNode] + C[idx[1] + (Cell - 1) * _strideCell])
        end
        c_star[1] = C[idx[1]]
        c_star[_nCells + 1] = C[idx[1] + _nCells * _strideCell - _strideNode]
        return nothing
    end

    # Compute strong derivative for auxiliary variable: Dg ← D * c
    @inline function auxiliaryVolumeIntegral!(state, idx, _strideNode::Int64, _strideCell::Int64, stateDer::Vector{Float64}, _nCells::Int64, _nNodes::Int64, _polyDerM::Matrix{Float64}, mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[idx[1] + (Cell-1) * _nNodes : idx[1] - 1 + _nNodes + (Cell-1) * _nNodes]))
            broadcast!(+, @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]), @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]), mul1)
        end
        nothing
    end

    # Lift face terms into auxiliary variable: Dg += M^{-1} B (c* - c)
    @inline function auxiliarySurfaceIntegral!(stateDer, state, idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer, _surfaceFlux, _nCells, _nNodes, _invMM, _polyDeg)
        for Cell in 1:_nCells
            @inbounds @simd for Node in 1:_nNodes
                stateDer[1 + ((Cell - 1) * strideCell_stateDer) + ((Node - 1) * strideNode_stateDer)] +=
                    (_invMM[Node, 1]) * ((state[idx[1] + ((Cell - 1) * strideCell_state)]) - (_surfaceFlux[Cell])) -
                    (_invMM[Node, end]) * ((state[idx[1] + ((Cell - 1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (_surfaceFlux[Cell + 1]))
            end
        end
        nothing
    end

    # Dc += (2/Δρ) * M_ρ^{-1} * (D^T * M^{(0,0)} * v * c - S_g * g)
    @inline function volumeIntegral!(Dc, y, idx, _strideNode::Int64, _strideCell::Int64, g, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _MM00::Matrix{Float64}, _invrMM::Vector{Matrix{Float64}}, S_g::Vector{Matrix{Float64}}, _deltarho::Float64, v, rho_i, mul1::Vector{Float64}, mul2::Vector{Float64})
        base = first(idx)
        # Convection term: M_ρ^{-1} * D^T * M^{(0,0)} * v * c
        @inbounds for Cell in 1:_nCells
            cell_idx = (Cell - 1) * _nNodes + 1 : Cell * _nNodes
            broadcast!(+, @view(Dc[cell_idx]), @view(Dc[cell_idx]), (2.0 / _deltarho) * v * (_invrMM[Cell] * (transpose(_polyDerM) * (_MM00 * @view(y[base + cell_idx[1] - 1 : base + cell_idx[end] - 1])))))
        end
        # Dispersion term: -M_ρ^{-1} * S_g * g
        @inbounds for Cell in 1:_nCells
            broadcast!(-, @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), @view(Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]), (2.0 / _deltarho) * (_invrMM[Cell] * (S_g[Cell] * @view(g[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]))))
        end
        return nothing
    end

    # d_rad_i: Vector of dispersion coefficients at each interface (length nCells+1)
    @inline function surfaceIntegral!(Dc, c_star::Vector{Float64}, g_star::Vector{Float64}, _nCells::Int, _nNodes::Int, _invrMM::Vector{Matrix{Float64}}, v::Float64, d_rad_i::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64)
        for Cell in 1:_nCells
            @inbounds @simd for Node in 1:_nNodes
                Dc[(Cell - 1) * _nNodes + Node] -= (2.0 / _deltarho) * (_invrMM[Cell][Node, 1] * (-v * c_star[Cell] + rho_i[Cell] * d_rad_i[Cell] * g_star[Cell]) + _invrMM[Cell][Node, _nNodes] * (v * c_star[Cell + 1] - rho_i[Cell + 1] * d_rad_i[Cell + 1] * g_star[Cell + 1]))
            end
        end
        return nothing
    end

    # Compute numerical fluxes c* and g*
    # d_rad_i: Vector of dispersion coefficients at each interface (length nCells+1)
    @inline function computeNumericalFluxes!(c_star, g_star, y, idx, g, _strideNode, _strideCell, _nNodes, _nCells, v, d_rad_i::Vector{Float64}, cIn)
        @inbounds for Cell in 2:_nCells
            c_star[Cell] = y[idx[1] + (Cell - 1) * _strideCell + ifelse(v >= 0.0, -_strideNode, 0)]
            g_star[Cell] = 0.5 * (g[(Cell - 1) * _nNodes] + g[(Cell - 1) * _nNodes + 1])
        end
        c_star[1] = cIn
        g_star[1] = v / d_rad_i[1] * (y[idx[1]] - cIn)

        c_star[_nCells + 1] = y[idx[1] + _nCells * _strideCell - _strideNode]
        g_star[_nCells + 1] = 0.0
        return nothing
    end

end