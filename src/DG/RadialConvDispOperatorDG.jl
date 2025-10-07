module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements

    abstract type ExactInt end
    struct exact_integration <: ExactInt end


    @inline function compute_faces_v(v_in::Number, rho_i::Vector{Float64}, rho_ip1::Vector{Float64}, rho_inner::Float64, nCells::Int)
        faces_v = Vector{Float64}(undef, nCells + 1)
        # v(r) = v_in * (rho_inner / r)
        faces_v[1] = v_in * (rho_inner / rho_inner)
        @inbounds for cell in 2:nCells
            faces_v[cell] = v_in * (rho_inner / rho_i[cell])
        end
        faces_v[nCells + 1] = v_in * (rho_inner / rho_ip1[nCells])
        return faces_v
    end

    @inline function compute_faces_v(v::Function, rho_i::Vector{Float64}, rho_ip1::Vector{Float64}, rho_inner::Float64, nCells::Int)
        faces_v = Vector{Float64}(undef, nCells + 1)
        faces_v[1] = v(rho_inner)
        @inbounds for cell in 2:nCells
            faces_v[cell] = v(rho_i[cell])
        end
        faces_v[nCells + 1] = v(rho_ip1[nCells])
        return faces_v
    end

    # --- Precompute D at faces ---
    @inline function diff_at_faces(d_rad::Number, rho_i::Vector{Float64}, rho_ip1::Vector{Float64})
        D_left  = fill(d_rad, length(rho_i))
        D_right = fill(d_rad, length(rho_ip1))
        return D_left, D_right
    end

    @inline function diff_at_faces(d_rad::Function, rho_i::Vector{Float64}, rho_ip1::Vector{Float64})
        return d_rad.(rho_i), d_rad.(rho_ip1)
    end

    @inline apply_outlet!(c_star, g_star, ::Val{:neumann0}, outlet_value, nCells) = (g_star[nCells+1] = 0.0)
    @inline apply_outlet!(c_star, g_star, ::Val{:dirichlet}, outlet_value, nCells) = (c_star[nCells+1] = outlet_value)

    @inline apply_outlet!(c_star, g_star, ::Val{:do_nothing}, outlet_value, nCells) = nothing

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _polyDeg, _invWeights, _nodes::Vector{Float64}, _polyDerM, _invMM, _MM, invMrhoM::Vector{<:AbstractMatrix}, SgMatrix::Union{Nothing,Vector{<:AbstractMatrix}}, v, d_rad, cIn, c_star, g_star, Dg, _h, mul1, _exactInt, outlet_bc::Symbol, outlet_value::Float64, rho_i::Vector{Float64}, rho_ip1::Vector{Float64}, faces_v::Union{Nothing,Vector{Float64}}, left_scale_vec::Union{Nothing,Vector{Float64}}, right_scale_vec::Union{Nothing,Vector{Float64}})
        fill!(Dg,0.0)   # reset auxiliary buffer used to build g
        fill!(Dc,0.0)   # reset residual accumulator for mobile phase
        
        # --- Function barriers / constants ---
        rho_inner = rho_i[1]
        map = 2.0 / _deltarho
        # fill!(_h,0.0)

        # Step 1) Build strong derivative in ξ: Dg ← D c
        volumeIntegraly!(y,idx, Dg,_nCells,_nNodes,_polyDerM,mul1)

        # Step 2) Numerical fluxes c*: upwind based on faces_v sign
        interfaceFluxUpwind!(c_star, y, (idx isa UnitRange ? first(idx) : idx), _strideNode, _strideCell, _nCells, faces_v)
        c_star[1] = cIn

        # Step 3) Lift face terms into g and scale to physical derivative
        surfaceIntegraly!(Dg,y,idx,_strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg,_invWeights,_exactInt)
        Dg .*= map

        # Step 4) Diffusive face traces g*
        interfaceFluxCentral!(g_star, Dg, 1, 1, _nNodes, _nCells)

        # Step 5) Outlet BC
        apply_outlet!(c_star, g_star, Val(outlet_bc), outlet_value, _nCells)

        # --- Surface term: -(2/Δρ) M_ρ^{-1} [ Lift(v c*) + Lift(g*) ] ---
        surfaceIntegral!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _invMM, invMrhoM, _polyDeg, _invWeights, _exactInt, c_star, g_star; left_scale = left_scale_vec, right_scale = right_scale_vec, _h = _h, tmp = mul1, faces_v = faces_v)
        # --- Volume term: (2/Δρ) M_ρ^{-1} [ Sᵀ (v c) - S_g g ] ---
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _deltarho, _nodes, _polyDeg, _polyDerM, _invMM, _MM, invMrhoM, SgMatrix, v, Dg, d_rad, rho_inner, _h, mul1)

        nothing
    end

    # compute strong derivative D c in ξ (per cell), precursor to g
    @inline function volumeIntegraly!(state,idx, stateDer::Vector{Float64},_nCells::Int64,_nNodes::Int64,_polyDerM::Matrix{Float64},mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            c_cell = @view state[idx[1] + (Cell-1) * _nNodes : idx[1] - 1 + _nNodes + (Cell-1) * _nNodes]
            mul!(mul1, _polyDerM, c_cell)  # mul1 = D * c_cell
            sDer_cell = @view stateDer[1 + (Cell-1) * _nNodes : Cell * _nNodes]
            @inbounds @simd for n in 1:_nNodes
                sDer_cell[n] -= mul1[n]
            end
        end
        nothing
    end

    # Volume term : (2/Δρ) M_ρ^{-1} [ Sᵀ (v c) - S_g g ]
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int64, _nNodes::Int64, _deltarho, _nodes::Vector{Float64}, _polyDeg::Int64, _polyDerM::Matrix{Float64}, _invMM::Matrix{Float64}, _MM::Matrix{Float64}, invMrhoM::Vector{<:AbstractMatrix}, SgMatrix, v, g, d_rad, rho_inner, _h, mul1::Vector{Float64})
        is_v_fun = v isa Function
        local_drad = (d_rad isa Function) ? d_rad : (r -> d_rad)
        map = 2.0 / _deltarho
        # local _MM = DGElements.MMatrix(_nodes, _polyDeg)
        @inbounds for Cell in 1:_nCells
            y_cell  = @view y[idx[1] + (Cell-1)*_nNodes : idx[1] - 1 + _nNodes + (Cell-1)*_nNodes]
            g_cell  = @view g[1 + (Cell-1)*_nNodes      : Cell*_nNodes]
            Dh_cell = @view Dc[1 + (Cell-1)*_nNodes      : Cell*_nNodes]
            tmpD    = @view _h[1 + (Cell-1)*_nNodes      : Cell*_nNodes]

            # tmpD ← v(ρ) .* c
            if is_v_fun
                rho0 = rho_inner + (Cell-1) * _deltarho
                @inbounds @simd for n in 1:_nNodes
                    r = rho0 + (_nodes[n] + 1.0) * (_deltarho/2)
                    tmpD[n] = v(r) * y_cell[n]
                end
            else
                # Treat scalar v as v_in and apply RFC profile v(r) = v_in * rho_inner / r
                rho0 = rho_inner + (Cell-1) * _deltarho
                v_in = v
                @inbounds @simd for n in 1:_nNodes
                    r = rho0 + (_nodes[n] + 1.0) * (_deltarho/2)
                    tmpD[n] = (v_in * rho_inner / r) * y_cell[n]
                end
            end
            # Correct SBP form: Sᵀ (v c) = Dᵀ * ( M * (v c) )
            # Build reference mass matrix M once per call (size depends only on nodes/polyDeg)
            mul!(mul1, _MM, tmpD)                                     # mul1 = M * (v c)
            mul!(tmpD, transpose(_polyDerM), mul1)                     # tmpD = Dᵀ * (M (v c)) = Sᵀ (v c)

            # Diffusion operator: mul1 = S_g * g_cell  (ρ D_rad weighted stiffness)
            let rho0 = rho_inner + (Cell-1) * _deltarho
                if SgMatrix !== nothing
                    mul!(mul1, SgMatrix[Cell], g_cell)
                else
                    Sg_cell = DGElements.weighted_stiff_Matrix(_nodes, _polyDeg, rho0, _deltarho, local_drad)
                    mul!(mul1, Sg_cell, g_cell)
                end
            end

            @inbounds @simd for n in 1:_nNodes
                tmpD[n] -= mul1[n]
            end
            mul!(mul1, invMrhoM[Cell], tmpD)                          # mul1 = M_ρ^{-1} ( ... )
            @inbounds @simd for n in 1:_nNodes
                Dh_cell[n] += map * mul1[n]
            end
        end
        nothing
    end

    # Lift face fluxes into interiors and accumulate in Dc:
    # Dc += -(2/Δρ) M_ρ^{-1} [ Lift(v c*) + Lift(g* with left/right scaling) ].
    @inline function surfaceIntegral!(Dc, y, idx, _strideNode::Int64, _strideCell::Int64, _nNodes::Int64, _nCells::Int64, _deltarho, _invMM::Matrix{Float64}, invMrhoM::Vector{<:AbstractMatrix}, _polyDeg, _invWeights, _exactInt, c_star, g_star; left_scale::Union{Float64,AbstractVector}, right_scale::Union{Float64,AbstractVector}, faces_v::Vector{Float64}, _h = nothing, tmp = nothing)
        if _h === nothing
            _h = similar(Dc)
        end
        if tmp === nothing
            tmp = similar(Dc)
        end
        fill!(_h, 0.0)

        # Fused lifting of (v*c*) and (g* with left/right scaling) in a single pass
        @inbounds for Cell in 1:_nCells
            lscale = (left_scale isa AbstractVector)  ? left_scale[Cell]  : left_scale
            rscale = (right_scale isa AbstractVector) ? right_scale[Cell] : right_scale

            vcL = faces_v[Cell]   * c_star[Cell]
            vcR = faces_v[Cell+1] * c_star[Cell+1]
            gL  = g_star[Cell]
            gR  = g_star[Cell+1]

            # Fill only first and last nodes; invMM will spread to all nodes.
            for Node in 1:_nNodes
                base = 1 + (Cell-1)*_nNodes + (Node-1)
                acc = 0.0
                if Node == 1
                    # +(vcL) - lscale*gL
                    acc +=  vcL - lscale * gL
                elseif Node == _nNodes
                    # -(vcR) + rscale*gR
                    acc += -vcR + rscale * gR
                end
                _h[base] += acc
            end
        end

        # Apply M^{-1}, then M_ρ^{-1}, and scale by -(2/Δρ)
        @inbounds for Cell in 1:_nCells
            h_cell  = @view _h[1 + (Cell-1)*_nNodes : Cell*_nNodes]
            Dh_cell = @view Dc[1 + (Cell-1)*_nNodes : Cell*_nNodes]
            mul!(tmp, _invMM, h_cell)
            mul!(h_cell, invMrhoM[Cell], tmp)
            @. Dh_cell += -(2.0 / _deltarho) * h_cell
        end

        nothing
    end

    # Exact-integration lifting of surface contributions
    # uses inverse-mass first/last columns to spread face fluxes to all nodes.
    @inline function surfaceIntegraly!(stateDer,state,idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg,_invWeights,_exactInt::exact_integration; left_scale=1.0, right_scale=1.0)
        for Cell in 1:_nCells
            lscale = (left_scale isa AbstractVector)  ? left_scale[Cell]  : left_scale
            rscale = (right_scale isa AbstractVector) ? right_scale[Cell] : right_scale
            @inbounds @simd for Node in 1:_nNodes
                stateDer[1 + ( (Cell-1) * strideCell_stateDer) + ( (Node-1) * strideNode_stateDer)] -=
                    (_invMM[Node, 1])  * ( (state[idx[1] + ((Cell-1) * strideCell_state)]) - (lscale*_surfaceFlux[Cell]) ) -
                    (_invMM[Node, end]) * ( (state[idx[1] + ((Cell-1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (rscale*_surfaceFlux[Cell+1]) )
            end
        end
        nothing
    end

        # Always-upwind numerical flux for advection (requires faces_v)
    @inline function interfaceFluxUpwind!(_surfaceFlux::Vector{Float64}, SRC, idx_start::Int, strideNode::Int, strideCell::Int, nCells::Int, faces_v::Vector{Float64})
        @inbounds begin
            @simd for Cell in 2:nCells
                vface = faces_v[Cell]
                cL = SRC[idx_start + (Cell-1) * strideCell - strideNode]
                cR = SRC[idx_start + (Cell-1) * strideCell]
                _surfaceFlux[Cell] = (vface >= 0) ? cL : cR
            end
            # copy boundary values
            _surfaceFlux[1] = SRC[idx_start]
            _surfaceFlux[nCells + 1] = SRC[idx_start + nCells * strideCell - strideNode]
        end
        return nothing
    end

    # Central numerical flux (used for diffusion auxiliary trace g*)
    @inline function interfaceFluxCentral!(_surfaceFlux::Vector{Float64}, SRC, idx_start::Int, strideNode::Int, strideCell::Int, nCells::Int)
        @inbounds begin
            @simd for Cell in 2:nCells
                l = SRC[idx_start + (Cell-1) * strideCell - strideNode]
                r = SRC[idx_start + (Cell-1) * strideCell]
                _surfaceFlux[Cell] = 0.5 * (l + r)
            end
            _surfaceFlux[1] = SRC[idx_start]
            _surfaceFlux[nCells + 1] = SRC[idx_start + nCells * strideCell - strideNode]
        end
        return nothing
    end

end