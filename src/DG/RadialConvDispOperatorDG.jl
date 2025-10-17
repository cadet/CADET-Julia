module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _invWeights, _polyDerM, _invMM, _MM00, _MM01, _nodes, v, d_rad, cIn, c_star, g_star, Dg, _h, mul1, rho_i, rho_ip1, radial_v, left_scale_vec, right_scale_vec)
        fill!(Dg, 0.0)   # reset auxiliary buffer used to build g
        fill!(Dc, 0.0)   # reset residual accumulator for mobile phase

        map = 2.0 / _deltarho

        # Build strong derivative in ξ: Dg ← -D * c
        volumeIntegraly!(y, idx, Dg, _nCells, _nNodes, _polyDerM, mul1)

        # Numerical fluxes c*
        interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells)

        # Lift face terms into g and scale to physical derivative
        surfaceIntegraly!(Dg, y, idx, _strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg, _invWeights)

        # g = (2/Δρ) * [ M^{-1}B(c* - c) - D c ]
        @. _h = map * Dg

        # Auxiliary face flux g*
        interfaceFluxCentral!(g_star, _h, 1, 1, _nNodes, _nCells)

        # Convective upwind flux c*
        interfaceFluxUpwind!(c_star, y, idx, _strideNode, _strideCell, _nCells, radial_v, cIn)

        # Surface term: B (v c*) - B_g (g*)
        surfaceIntegral!(Dc, _nNodes, _nCells, c_star, g_star, left_scale_vec, right_scale_vec)

        # Volume term : Sᵀ (v c) - S_g g
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _polyDerM, _MM00, _MM01, radial_v, _h, mul1, rho_i, _deltarho, d_rad, _nodes, v, rho_i[1])        
        
        @inbounds for Cell in 1:_nCells
            @views Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes] = ( rho_i[Cell] .* _MM00 .+ (_deltarho/2) .* _MM01 ) \ Dc[(Cell - 1) * _nNodes + 1 : Cell * _nNodes]
        end

        return nothing
    end

    # compute strong derivative - D c in ξ, precursor to g (Eq. 84b)
    @inline function volumeIntegraly!(state, idx, stateDer::Vector{Float64}, _nCells::Int64, _nNodes::Int64, _polyDerM::Matrix{Float64}, mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes])) 
            # @. stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)] -= mul1 #-D*c
            broadcast!(-, @view(stateDer[1 + (Cell-1) * _nNodes : _nNodes + (Cell-1) * _nNodes]),@view(stateDer[(1 + (Cell-1) * _nNodes) : (_nNodes + (Cell-1) * _nNodes)]), mul1)
        end
        nothing
    end

    # Compute strong derivative of the Volume term: Dᵀ M(0,0) (v c) - S_g g
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int, _nNodes::Int, _polyDerM::Matrix{Float64}, _MM00::Matrix{Float64}, _MM01::Matrix{Float64}, radial_v::Vector{Float64}, g::Vector{Float64}, mul1::Vector{Float64}, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Float64, _nodes::Vector{Float64}, v_in::Float64, rho_inner::Float64)        
        
        @inbounds for Cell in 1:_nCells
            c_cell  = @view y[idx[1] + (Cell - 1) * _nNodes : idx[1] - 1 + Cell * _nNodes]
            g_cell  = @view g[1 + (Cell - 1) * _nNodes : Cell * _nNodes]
            Dc_cell = @view Dc[1 + (Cell - 1) * _nNodes : Cell * _nNodes]
            
            # convective volume: Sᵀ (v ∘ c)  with Sᵀ = Dᵀ * M(0,0)
            # nodal radii in this cell: r(ξ) = ρ_i + (Δρ/2)*(1+ξ)
            r_nodes = @. rho_i[Cell] + (_deltarho/2) * (1.0 + _nodes)
            @. mul1 = (v_in * rho_inner) * (c_cell / r_nodes)   # elementwise

            # convective volume: Sᵀ ( (v ∘ c) )
            tmp_conv = transpose(_polyDerM) * (_MM00 * mul1)

            # diffusive volume: S_g g with S_g = Dᵀ * M_ρ * d_rad
            tmp_diff = transpose(_polyDerM) * ( (rho_i[Cell] .* _MM00 .+ (_deltarho/2) .* _MM01) * (d_rad .* g_cell) )
            
            # Accumulate volume
            @. Dc_cell += (tmp_conv - tmp_diff)
        end

        return nothing
    end

    # Exact-integration lifting of surface contributions
    @inline function surfaceIntegraly!(stateDer,state,idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg,_invWeights)
        for Cell in 1:_nCells
            @inbounds @simd for Node in 1:_nNodes
                stateDer[1 + ( (Cell - 1) * strideCell_stateDer) + ( (Node - 1) * strideNode_stateDer)] -=
                    (_invMM[Node, 1])  * ( (state[idx[1] + ((Cell - 1) * strideCell_state)]) - (_surfaceFlux[Cell]) ) -
                    (_invMM[Node, end]) * ( (state[idx[1] + ((Cell - 1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (_surfaceFlux[Cell + 1]) )
            end
        end
        nothing
    end

    # Exact-Integration lifting of Surface Contribution += B (v c*) - B_g (g*)
    @inline function surfaceIntegral!(Surf::Vector{Float64}, nNodes::Int, nCells::Int, c_star::Vector{Float64}, g_star::Vector{Float64}, left_scale::Vector{Float64}, right_scale::Vector{Float64})
        @inbounds for Cell in 1:nCells
            # endpoint indices for this cell
            iL = (Cell - 1) * nNodes + 1
            iR = Cell * nNodes
            # convective: +B (v c*)
            convL =  c_star[Cell]
            convR = - c_star[Cell + 1]
            # diffusive: - B_g (g*)
            diffL = - (left_scale[Cell]  * g_star[Cell])
            diffR = (right_scale[Cell] * g_star[Cell + 1])
            # accumulate
            Surf[iL] += (convL + diffL)
            Surf[iR] += (convR + diffR)
        end
        return nothing
    end

    # Auxiliary flux for g: c* = 0.5 (c_L + c_R)
    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C, idx::UnitRange{Int}, strideNode::Int, strideCell::Int, nCells::Int)
        base = first(idx)
        @inbounds for Cell in 2:nCells
                cL = C[base + (Cell - 1) * strideCell - strideNode]
                cR = C[base + (Cell - 1) * strideCell]
                _surfaceFlux[Cell] = 0.5 * (cL + cR)
            end
            # boundaries
            _surfaceFlux[1] = C[base]
            _surfaceFlux[nCells + 1] = C[base + nCells * strideCell - strideNode]
        return nothing
    end

    # Upwind numerical flux for c:
    @inline function interfaceFluxUpwind!(_surfaceFlux::Vector{Float64}, C, idx::UnitRange{Int}, strideNode::Int, strideCell::Int, nCells::Int, radial_v::Vector{Float64}, cIn::Float64, cOut::Float64 = NaN)
        base = first(idx)
        @inbounds begin
            # interior faces
            for cell in 2:nCells
                v  = radial_v[cell]
                cL = C[base + (cell - 1) * strideCell - strideNode]
                cR = C[base + (cell - 1) * strideCell]
                _surfaceFlux[cell] = v >= 0 ? v * cL : v * cR
            end
            # inlet
            _surfaceFlux[1] = radial_v[1] >= 0 ? radial_v[1] * cIn : radial_v[1] * C[base]
            # outlet
            if radial_v[nCells + 1] >= 0
                _surfaceFlux[nCells + 1] = radial_v[nCells + 1] * C[base + nCells * strideCell - strideNode]
            else
                _surfaceFlux[nCells + 1] = isnan(cOut) ? radial_v[nCells + 1] * C[base + nCells * strideCell - strideNode] : radial_v[nCells + 1] * cOut
            end
        end
        return nothing
    end

    # Central numerical flux of g (auxiliary trace g*)
    @inline function interfaceFluxCentral!(_surfaceFlux, SRC, idx_start, strideNode, strideCell, nCells)
        base = idx_start
        @inbounds begin
            for Cell in 2:nCells
                gL = SRC[base + (Cell - 1) * strideCell - strideNode]
                gR = SRC[base + (Cell - 1) * strideCell]
                _surfaceFlux[Cell] = 0.5 * (gL + gR)
            end
            _surfaceFlux[1]          = SRC[base]
            _surfaceFlux[nCells + 1] = SRC[base + nCells * strideCell - strideNode]
        end
        return nothing
    end

end

