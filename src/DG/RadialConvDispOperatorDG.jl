module RadialConvDispOperatorDG
    using LinearAlgebra
    using CADETJulia: DGElements
    # Governing PDE:
    #   ∂c/∂t = -(v/ρ) ∂c/∂ρ + (1/ρ) ∂/∂ρ ( ρ D_rad ∂c/∂ρ ) + ((1-ε_c)/ε_c) * (3 k_f / R_p) * ( c^(p)|_{ρ=R_p} - c )
    #
    # Mapping each cell (ρ_i, ρ_{i+1}) → ξ ∈ (-1,1) with ρ(ξ)=ρ_i + (ξ+1)Δρ_i/2 and introducing the auxiliary g := (2/Δρ_i) ∂c/∂ξ
    # The semi-discrete system (elementwise):
    #
    #   ċ = (2/Δρ_i) M_ρ^{-1} [ Sᵀ (v c) - S_g g - ( B (v c*) - B_g g* ) ] - ((1-ε_c)/ε_c) * (3/R_p) * M_ρ^{-1} M_K ( c - c^(p) )
    #   g  = (2/Δρ_i) [ M^{-1} B (c* - c) - D c ]
    #
    # where
    #   M_ρ   = (Δρ_i/2) M^(0,1) + ρ_i M^(0,0)        (ρ-weighted mass matrix),
    #   M     = M^(0,0)                               (mass matrix),
    #   Sᵀ    = Dᵀ M                                  (SBP identity),
    #   S_g   ≔ ρ D_rad-weighted stiffness (quadrature)
    #   B, B_g lift numerical face fluxes
    #   c*, g* are numerical traces

    abstract type ExactInt end
    # A tag type selecting exact integration on surfaces (“exact LGL”).
    struct exact_integration <: ExactInt end

    # Assemble Dc = advection + diffusion + film exchange.
    #
    # Inputs:
    #   y[idx]      nodal c in current radial component block
    #   Dg          workspace for auxiliary g ≈ (2/Δρ) ∂c/∂ξ after lifting
    #   invMrhoM    per-cell M_ρ^{-1}
    #   SgMatrix    ρ D_rad stiffness
    #   MKMatrix    film quadrature
    #   v           scalar v or function r ↦ v(r)
    #   d_rad       scalar D_rad or function r ↦ D_rad(r)
    #   c_star,g_star  face traces (c*, g*)
    #   outlet_bc   :neumann0 → g*|_right=0 (zero gradient), :dirichlet → c*|_right=const

    # --- Precompute face velocities (two specialized methods) ---
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

    # --- Outlet BC  ---
    @inline apply_outlet!(c_star, g_star, ::Val{:neumann0}, outlet_value, nCells) = (g_star[nCells+1] = 0.0)
    @inline apply_outlet!(c_star, g_star, ::Val{:dirichlet}, outlet_value, nCells) = (c_star[nCells+1] = outlet_value)
    # Fallback: do nothing
    @inline apply_outlet!(c_star, g_star, ::Val{:do_nothing}, outlet_value, nCells) = nothing

    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nPoints, _nNodes, _nCells, _deltarho, _polyDeg, _invWeights, _nodes::Vector{Float64}, _polyDerM, _invMM, invMrhoM::Vector{<:AbstractMatrix}, SgMatrix::Vector{<:AbstractMatrix}, MKMatrix, v, d_rad, cIn, c_star, g_star, Dg, _h, mul1, _exactInt, Rp, kf, cp, outlet_bc::Symbol, outlet_value::Float64, rho_i::Vector{Float64}, rho_ip1::Vector{Float64}, rho_inner; faces_v::Union{Nothing,Vector{Float64}}=nothing, left_scale_vec::Union{Nothing,Vector{Float64}}=nothing, right_scale_vec::Union{Nothing,Vector{Float64}}=nothing)
        # Reset auxiliary arrays: Dg builds g; Dc accumulates RHS
        fill!(Dg,0.0)   # reset auxiliary buffer used to build g
        fill!(Dc,0.0)   # reset residual accumulator for mobile phase
        # --- Function barriers / constants ---
        map = 2.0 / _deltarho
        use_faces_v = faces_v !== nothing
        # fill!(_h,0.0)
        # --- Velocity at faces ---
        if faces_v === nothing
            faces_v = compute_faces_v(v, rho_i, rho_ip1, rho_inner, _nCells)
        end

        # Step 1) Build strong derivative in ξ: Dg ← D c
        volumeIntegraly!(y,idx, Dg,_nCells,_nNodes,_polyDerM,mul1)

        # Step 2) Numerical fluxes c*: central or upwind based on faces_v sign
        interfaceFluxAuxiliary!(c_star, y, (idx isa UnitRange ? first(idx) : idx), _strideNode, _strideCell, _nCells; faces_v=faces_v)
        c_star[1] = cIn

        # Step 3) Lift face terms into g and scale to physical derivative
        surfaceIntegraly!(Dg,y,idx,_strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg,_invWeights,_exactInt)
        Dg .*= map

        # Step 4) Diffusive face traces g*
        interfaceFluxAuxiliary!(g_star, Dg, 1, 1, _nNodes, _nCells)

        # Step 5) Outlet BC
        apply_outlet!(c_star, g_star, Val(outlet_bc), outlet_value, _nCells)

        # --- Left/right diffusion scales ---
        if left_scale_vec === nothing || right_scale_vec === nothing
            D_left, D_right = diff_at_faces(d_rad, rho_i, rho_ip1)
            left_scale_vec  = @. - rho_i  * D_left
            right_scale_vec = @. + rho_ip1 * D_right
        end

        # --- Surface term: -(2/Δρ) M_ρ^{-1} [ Lift(v c*) + Lift(g*) ] ---
        surfaceIntegral!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _invMM, invMrhoM, _polyDeg, _invWeights, _exactInt,
                         c_star, g_star; left_scale = left_scale_vec, right_scale = right_scale_vec, v = 0.0, _h = _h, tmp = mul1, faces_v = faces_v)

        # --- Volume term: (2/Δρ) M_ρ^{-1} [ Sᵀ (v c) - S_g g ] ---
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _deltarho, _nodes, _polyDeg, _polyDerM, _invMM, invMrhoM, SgMatrix, v, Dg, d_rad, rho_inner, _h, mul1)

        # --- Film exchange term ---
        filmvolumeIntegral!(Dc, y, cp, idx, _nCells, _nNodes, _nodes, _deltarho, rho_inner, invMrhoM, MKMatrix, Rp, kf, _h, mul1)
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
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int64, _nNodes::Int64, _deltarho, _nodes::Vector{Float64}, _polyDeg::Int64, _polyDerM::Matrix{Float64}, _invMM::Matrix{Float64}, invMrhoM::Vector{<:AbstractMatrix}, SgMatrix, v, g, d_rad, rho_inner, _h, mul1::Vector{Float64})
        # (idea 8) Function barriers / constants
        is_v_fun = v isa Function
        local_drad = (d_rad isa Function) ? d_rad : (r -> d_rad)
        map = 2.0 / _deltarho
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
            local _MM = DGElements.MMatrix(_nodes, _polyDeg)
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

    # Film Diffusion term:
    # Dc += (3/R_p) * M_ρ^{-1} * MK * ( c^(p) - c )
    # MK = diag( w_m * ρ(ξ_m) * k_f(ρ(ξ_m)) ) using LGL quadrature points ξ_m.
    @inline function filmvolumeIntegral!(Dc, y, cp, idx, _nCells::Int64, _nNodes::Int64, _nodes::Vector{Float64}, _deltarho, rho_inner, invMrhoM::Vector{<:AbstractMatrix}, MKMatrix, Rp, kf, _h, mul1::Vector{Float64})
        pref = (3.0 / Rp)
        if pref == 0.0
            return  # usual LRMP/GRM path: Rp=Inf ⇒ no film source added here
        end
        @inbounds for Cell in 1:_nCells
            y_cell  = @view y[idx[1] + (Cell-1)*_nNodes : idx[1]-1 + _nNodes + (Cell-1)*_nNodes]
            cp_cell = @view cp[idx[1] + (Cell-1)*_nNodes : idx[1]-1 + _nNodes + (Cell-1)*_nNodes]
            Dh_cell = @view Dc[1 + (Cell-1)*_nNodes      : Cell*_nNodes]

            @. mul1 = cp_cell - y_cell                               # (c^(p) - c) at nodes

            # Build and apply MK_cell = diag( w(ξ_m) * ρ(ξ_m) * k_f(ρ(ξ_m)) )
            let rho0 = rho_inner + (Cell-1) * _deltarho
                if MKMatrix !== nothing
                    mul!(_h, MKMatrix[Cell], mul1)
                else
                    local_kf = (kf isa Function) ? (hatrho -> hatrho .* kf.(hatrho)) : (hatrho -> hatrho .* kf)
                    MK_cell = DGElements.weightedQuadrature(_nodes, rho0, _deltarho, local_kf)
                    mul!(_h, MK_cell, mul1)
                end
            end
            @. Dh_cell += pref * mul1
        end
        nothing
    end

    # Unified interfaceFluxAuxiliary! for both advection (c*) and diffusion (g*) traces
    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, SRC, idx_start::Int, strideNode::Int, strideCell::Int, nCells::Int; mode::Symbol=:auto, faces_v::Union{Nothing,Vector{Float64}}=nothing)
        central = (mode === :central) || (mode === :auto && faces_v === nothing)
        @inbounds begin
            @simd for Cell in 2:nCells
                if central
                    l = SRC[idx_start + (Cell-1) * strideCell - strideNode]
                    r = SRC[idx_start + (Cell-1) * strideCell]
                    _surfaceFlux[Cell] = 0.5 * (l + r)
                else
                    vface = faces_v[Cell]
                    cL = SRC[idx_start + (Cell-1) * strideCell - strideNode]
                    cR = SRC[idx_start + (Cell-1) * strideCell]
                    _surfaceFlux[Cell] = (vface >= 0) ? cL : cR
                end
            end
            _surfaceFlux[1] = SRC[idx_start]
            _surfaceFlux[nCells + 1] = SRC[idx_start + nCells * strideCell - strideNode]
        end
        return nothing
    end

    # Lift face fluxes into interiors and accumulate in Dc:
    # Dc += -(2/Δρ) M_ρ^{-1} [ Lift(v c*) + Lift(g* with left/right scaling) ].
    @inline function surfaceIntegral!(Dc, y, idx, _strideNode::Int64, _strideCell::Int64, _nNodes::Int64, _nCells::Int64, _deltarho, _invMM::Matrix{Float64}, invMrhoM::Vector{<:AbstractMatrix}, _polyDeg, _invWeights, _exactInt, c_star, g_star; left_scale= -1.0, right_scale= +1.0, v=0.0, _h = nothing, tmp = nothing, faces_v::Union{Nothing,Vector{Float64}}=nothing)
        if _h === nothing
            _h = similar(Dc)
        end
        if tmp === nothing
            tmp = similar(Dc)
        end
        fill!(_h, 0.0)

        # Fused lifting of (v*c*) and (g* with left/right scaling) in a single pass
        @inbounds for Cell in 1:_nCells
            # Resolve scalar/Vector scales for this cell
            lscale = (left_scale isa AbstractVector)  ? left_scale[Cell]  : left_scale
            rscale = (right_scale isa AbstractVector) ? right_scale[Cell] : right_scale

            # Pre-fetch boundary values for this cell (face indices Cell and Cell+1)
            # vc_left/right: either from faces_v (multiply here), or compute using scalar v
            vcL = (faces_v === nothing) ? (v * c_star[Cell])     : (faces_v[Cell]   * c_star[Cell])
            vcR = (faces_v === nothing) ? (v * c_star[Cell+1])   : (faces_v[Cell+1] * c_star[Cell+1])
            gL  = g_star[Cell]
            gR  = g_star[Cell+1]

            # Assemble raw face lifts using canonical basis vectors; only first/last nodes get direct face contributions here.
            @inbounds for Node in 1:_nNodes
                base = 1 + (Cell-1)*_nNodes + (Node-1)
                # Start with zero; only first/last nodes get direct face contributions here.
                acc = 0.0
                if Node == 1
                    # Left face lift: -( -vcL ) and -( -lscale*gL )
                    acc -= ( -vcL )
                    acc -= ( -lscale * gL )
                elseif Node == _nNodes
                    # Right face lift: -( -vcR ) and -( -rscale*gR )
                    acc -= ( -vcR )
                    acc -= ( -rscale * gR )
                end
                _h[base] += acc
            end
        end

        # Apply M^{-1} and then M_ρ^{-1} cellwise, and scale with -(2/Δρ)
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


end