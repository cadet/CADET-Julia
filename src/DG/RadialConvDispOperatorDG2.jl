module RadialConvDispOperatorDG
    using LinearAlgebra
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
    #   S     = Dᵀ M                                  (SBP identity),
    #   S_g   ≔ ρ D_rad-weighted stiffness (built by quadrature),
    #   B, B_g lift numerical face fluxes
    #   c*, g* are numerical traces
    # ------------------------------------------------------------------------

    abstract type ExactInt end
    # A tag type selecting exact integration on surfaces (matches PDF’s “exact LGL”).
    struct exact_integration <: ExactInt end

    # ------------------------------------------------------------------------
    # Assemble Dc = advection + diffusion + film exchange.
    #
    # Inputs:
    #   y[idx]      nodal c in current radial component block
    #   Dg          workspace for auxiliary g ≈ (2/Δρ) ∂c/∂ξ after lifting
    #   invMrhoM    per-cell M_ρ^{-1}
    #   SgMatrix    optional cached ρ D_rad stiffness
    #   MKMatrix    optional cached film quadrature
    #   v           scalar v or function r ↦ v(r)
    #   d_rad       scalar D_rad or function r ↦ D_rad(r)
    #   c_star,g_star  face traces (c*, g*)
    #   outlet_bc   :neumann0 → g*|_right=0 (zero gradient), :dirichlet → c*|_right=const
    # ------------------------------------------------------------------------
    @inline function radialresidualImpl!(Dc, y, idx, _strideNode, _strideCell, _nPoints, _nNodes, _nCells, _deltarho, _polyDeg, _invWeights, _nodes::Vector{Float64}, _polyDerM, _invMM, invMrhoM::Vector{<:AbstractMatrix}, SgMatrix::Vector{<:AbstractMatrix}, MKMatrix, advection_flux::Symbol, v, d_rad, cIn, c_star, g_star, Dg, _h, mul1, _exactInt, eps_c, Rp, kf, cp, outlet_bc::Symbol, outlet_value::Float64, rho_i::Vector{Float64}, rho_ip1::Vector{Float64}, rho_inner)
        # Reset auxiliary arrays: Dg builds g; Dc accumulates RHS
        fill!(Dg,0.0)   # reset auxiliary buffer used to build g
        fill!(Dc,0.0)   # reset residual accumulator for mobile phase

        # Step 1) Build strong derivative in ξ: Dg ← D c (pre-SBP)
        volumeIntegraly!(y,idx, Dg,_nCells,_nNodes,_polyDerM,mul1)

        # Step 2) Numerical traces c*. If v(r) known per face, enable upwind.
        if v isa Function
            faces_v = Vector{Float64}(undef, _nCells+1)
            faces_v[1] = v(rho_inner)
            @inbounds for Cell in 2:_nCells
                faces_v[Cell] = v(rho_i[Cell])
            end
            faces_v[_nCells+1] = v(rho_ip1[_nCells])
            interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells; advection_flux=advection_flux, faces_v=faces_v)
        else
            interfaceFluxAuxiliary!(c_star, y, idx, _strideNode, _strideCell, _nCells; advection_flux=advection_flux)
        end

        c_star[1] = cIn  # enforce Dirichlet inlet at ρ_in (left boundary of global domain)

        # Step 3) SBP lifting of face terms into g: Dg ← Dg + Lift( faces of c )
        surfaceIntegraly!(Dg,y,idx,_strideNode, _strideCell, 1, _nNodes, c_star, _nCells, _nNodes, _invMM, _polyDeg,_invWeights,_exactInt)

        @. Dg *= (2.0 / _deltarho)  # map ∂/∂ξ → (2/Δρ) ∂/∂ρ, so g = (2/Δρ) ∂c/∂ξ

        # Step 4) Build numerical traces g* (central), for diffusive flux ρ D_rad g
        interfaceFluxAuxiliary!(g_star, Dg, 1, _nNodes, _nCells)

        # Step 5) Outlet boundary condition: Neumann(0) or Dirichlet
        if outlet_bc === :neumann0
            g_star[_nCells + 1] = 0.0
        elseif outlet_bc === :dirichlet
            c_star[_nCells + 1] = outlet_value
        end

        # Step 6) Diffusive face scalings:  (B_g uses ± ρ D at faces)
        D_left  = (d_rad isa Function) ? d_rad.(rho_i)  : fill(d_rad, _nCells)
        D_right = (d_rad isa Function) ? d_rad.(rho_ip1) : fill(d_rad, _nCells)
        left_scale_vec  = @. - rho_i  * D_left   # left face:  -ρ_i D(ρ_i)
        right_scale_vec = @. + rho_ip1 * D_right # right face: +ρ_{i+1} D(ρ_{i+1})

        # Step 7) Surface term:
        # Dc += -(2/Δρ) M_ρ^{-1} [ B (v c*) + B_g g* ]  (faces_vc holds v_face * c*)
        if v isa Function
            faces_v = Vector{Float64}(undef, _nCells+1)
            faces_v[1] = v(rho_inner)
            @inbounds for Cell in 2:_nCells
                faces_v[Cell] = v(rho_i[Cell])
            end
            faces_v[_nCells+1] = v(rho_ip1[_nCells])
            faces_vc = similar(faces_v)
            @inbounds for i in 1:_nCells+1
                faces_vc[i] = faces_v[i] * c_star[i]
            end
            surfaceIntegral!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _invMM, invMrhoM, _polyDeg, _invWeights, _exactInt, c_star, g_star, left_scale = left_scale_vec, right_scale = right_scale_vec, v = 0.0, _h = _h, tmp = mul1, faces_vc = faces_vc)
        else
            surfaceIntegral!(Dc, y, idx, _strideNode, _strideCell, _nNodes, _nCells, _deltarho, _invMM, invMrhoM, _polyDeg, _invWeights, _exactInt, c_star, g_star, left_scale = left_scale_vec, right_scale = right_scale_vec, v = v, _h = _h, tmp = mul1)
        end

        # Step 8) Volume term:
        # Dc += (2/Δρ) M_ρ^{-1} [ Sᵀ (v c) - S_g g ]
        volumeIntegral!(Dc, y, idx, _nCells, _nNodes, _deltarho, _nodes, _polyDeg, _polyDerM, _invMM, invMrhoM, SgMatrix, v, Dg, d_rad, rho_inner, _h, mul1)

        # Step 9) Film exchange term:
        # Dc += ((1-ε_c)/ε_c) * (3/R_p) M_ρ^{-1} MK ( c^(p) - c )
        # Note: callers often neutralize this via Rp=Inf (pref=0) and/or ε_c=1.
        filmvolumeIntegral!(Dc, y, cp, idx, _nCells, _nNodes, _nodes, _deltarho, rho_inner, invMrhoM, MKMatrix, eps_c, Rp, kf, _h, mul1)

        # Dc now holds the semi-discrete RHS for the mobile phase
        nothing
    end

    # ------------------------------------------------------------------------
    # compute strong derivative D c in ξ (per cell), precursor to g
    @inline function volumeIntegraly!(state,idx, stateDer::Vector{Float64},_nCells::Int64,_nNodes::Int64,_polyDerM::Matrix{Float64},mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            mul!(mul1, _polyDerM, @view(state[idx[1] + (Cell-1) * _nNodes : idx[1] - 1 + _nNodes + (Cell-1) * _nNodes]))  # mul1 = D * c_cell
            broadcast!((a,b)->a+b, @view(stateDer[1 + (Cell-1) * _nNodes : Cell * _nNodes]), @view(stateDer[1 + (Cell-1) * _nNodes : Cell * _nNodes]), mul1)  # accumulate into g
        end
        nothing
    end

    # ------------------------------------------------------------------------
    # Volume contribution to Eq. (14a): (2/Δρ) M_ρ^{-1} [ Sᵀ (v c) - S_g g ]
    @inline function volumeIntegral!(Dc, y, idx, _nCells::Int64, _nNodes::Int64, _deltarho, _nodes::Vector{Float64}, _polyDeg::Int64, _polyDerM::Matrix{Float64}, _invMM::Matrix{Float64}, invMrhoM::Vector{<:AbstractMatrix}, SgMatrix, v, g, d_rad, rho_inner, _h, mul1::Vector{Float64})
        @inbounds for Cell in 1:_nCells
            y_cell  = @view y[idx[1] + (Cell-1)*_nNodes : idx[1] - 1 + _nNodes + (Cell-1)*_nNodes]
            g_cell  = @view g[1 + (Cell-1)*_nNodes      : Cell*_nNodes]
            Dh_cell = @view Dc[1 + (Cell-1)*_nNodes      : Cell*_nNodes]
            tmpD    = @view _h[1 + (Cell-1)*_nNodes      : Cell*_nNodes]

            # tmpD ← v(ρ) .* c (sample v at nodes if functional)
            if v isa Function
                rho_i = rho_inner + (Cell-1) * _deltarho
                @inbounds for n in 1:_nNodes
                    tmpD[n] = v(rho_i + (_nodes[n] + 1.0) * (_deltarho/2)) * y_cell[n]
                end
            else
                @. tmpD = v * y_cell
            end
            ldiv!(mul1, _invMM, tmpD)                                 # mul1 = M^{-1} * (v c)
            mul!(tmpD, transpose(_polyDerM), mul1)                     # tmpD = Dᵀ (M^{-1} (v c)) = Sᵀ (v c)

            # Diffusion operator: mul1 = S_g * g_cell  (ρ D_rad weighted stiffness)
            if d_rad isa Function
                rho_i = rho_inner + (Cell-1) * _deltarho
                Sg_cell = DGElements.weighted_stiff_Matrix(_nodes, _polyDeg, rho_i, _deltarho, d_rad)
                mul!(mul1, Sg_cell, g_cell)
            elseif d_rad isa Number
                rho_i = rho_inner + (Cell-1) * _deltarho
                Sg_cell = DGElements.weighted_stiff_Matrix(_nodes, _polyDeg, rho_i, _deltarho, _ -> d_rad)
                mul!(mul1, Sg_cell, g_cell)
            elseif SgMatrix !== nothing
                mul!(mul1, SgMatrix[Cell], g_cell)
            else
                rho_i = rho_inner + (Cell-1) * _deltarho
                if d_rad isa Function
                    Sg_cell = DGElements.weighted_stiff_Matrix(_nodes, _polyDeg, rho_i, _deltarho, d_rad)
                else
                    Sg_cell = DGElements.weighted_stiff_Matrix(_nodes, _polyDeg, rho_i, _deltarho, _ -> d_rad)
                end
                mul!(mul1, Sg_cell, g_cell)
            end

            @. tmpD -= mul1                                           # tmpD = Sᵀ (v c) - S_g g
            mul!(mul1, invMrhoM[Cell], tmpD)                          # mul1 = M_ρ^{-1} ( ... )
            @. Dh_cell += (2.0 / _deltarho) * mul1                    # scale by (2/Δρ) (mapping)
        end
        nothing
    end

    # ------------------------------------------------------------------------
    # Film exchange contribution:
    # Dc += ((1-ε_c)/ε_c) * (3/R_p) * M_ρ^{-1} * MK * ( c^(p) - c )
    # MK = diag( w_m * ρ(ξ_m) * k_f(ρ(ξ_m)) ) using LGL quadrature points ξ_m.
    @inline function filmvolumeIntegral!(Dc, y, cp, idx, _nCells::Int64, _nNodes::Int64, _nodes::Vector{Float64}, _deltarho, rho_inner, invMrhoM::Vector{<:AbstractMatrix}, MKMatrix, eps_c, Rp, kf, _h, mul1::Vector{Float64})
        pref = (3.0 / Rp)
        if pref == 0.0
            return  # usual LRMP/GRM path: Rp=Inf ⇒ no film source added here
        end
        @inbounds for Cell in 1:_nCells
            y_cell  = @view y[idx[1] + (Cell-1)*_nNodes : idx[1]-1 + _nNodes + (Cell-1)*_nNodes]
            cp_cell = @view cp[idx[1] + (Cell-1)*_nNodes : idx[1]-1 + _nNodes + (Cell-1)*_nNodes]
            Dh_cell = @view Dc[1 + (Cell-1)*_nNodes      : Cell*_nNodes]

            @. mul1 = cp_cell - y_cell                               # (c^(p) - c) at nodes

            # Build MK_cell = diag( w(ξ_m) * ρ(ξ_m) * k_f(ρ(ξ_m)) )
            if kf isa Function
                rho_i = rho_inner + (Cell-1) * _deltarho
                MK_cell = DGElements.weightedQuadrature(_nodes, rho_i, _deltarho, hatrho -> hatrho .* kf.(hatrho))
                mul!(_h, MK_cell, mul1)
            elseif kf isa Number
                rho_i = rho_inner + (Cell-1) * _deltarho
                MK_cell = DGElements.weightedQuadrature(_nodes, rho_i, _deltarho, hatrho -> hatrho .* kf)
                mul!(_h, MK_cell, mul1)
            elseif MKMatrix !== nothing
                mul!(_h, MKMatrix[Cell], mul1)
            else
                rho_i = rho_inner + (Cell-1) * _deltarho
                if kf isa Function
                    MK_cell = DGElements.weightedQuadrature(_nodes, rho_i, _deltarho, hatrho -> hatrho .* kf.(hatrho))
                    mul!(_h, MK_cell, mul1)
                elseif kf isa Number
                    MK_cell = DGElements.weightedQuadrature(_nodes, rho_i, _deltarho, hatrho -> hatrho .* kf)
                    mul!(_h, MK_cell, mul1)
                else
                    # no film term if MK not provided and no kf information
                end
            end

            # Apply ε_c scaling and M_ρ^{-1}:
            # ((1-ε_c)/ε_c) * (3/R_p) * M_ρ^{-1} * MK * (c^(p) - c)
            @. _h = ((1 - eps_c) / eps_c) * _h
            mul!(mul1, invMrhoM[Cell], _h)
            @. Dh_cell += pref * mul1
        end
        nothing
    end

    # ------------------------------------------------------------------------
    # Face traces c* for advection:
    #   central: average; upwind: pick by sign of face velocity.
    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, C, idx, _strideNode::Int64, _strideCell::Int64, _nCells::Int64; advection_flux::Symbol=:central, faces_v::Union{Nothing,Vector{Float64}}=nothing)
        if advection_flux == :central || faces_v === nothing
            @inbounds for Cell in 2:_nCells
                _surfaceFlux[Cell] = 0.5 * ((C[idx[1] + (Cell-1) * _strideCell - _strideNode]) + (C[idx[1] + (Cell-1) * _strideCell]))
            end
            _surfaceFlux[1] = 0.5 * ((C[idx[1] ]) + (C[idx[1] ]))
            _surfaceFlux[_nCells+1] = 0.5 * ((C[idx[1] + _nCells * _strideCell - _strideNode]) + (C[idx[1] + _nCells * _strideCell - _strideNode]))
        else
            @inbounds for Cell in 2:_nCells
                vface = faces_v[Cell]
                cL = C[idx[1] + (Cell-1) * _strideCell - _strideNode]
                cR = C[idx[1] + (Cell-1) * _strideCell]
                _surfaceFlux[Cell] = (vface >= 0) ? cL : cR
            end
            vface_left = faces_v[1]
            cL = C[idx[1]]
            cR = C[idx[1]]
            _surfaceFlux[1] = (vface_left >= 0) ? cL : cR
            vface_right = faces_v[end]
            cL = C[idx[1] + _nCells * _strideCell - _strideNode]
            cR = C[idx[1] + _nCells * _strideCell - _strideNode]
            _surfaceFlux[_nCells+1] = (vface_right >= 0) ? cL : cR
        end
        nothing
    end

    # Face traces g* for diffusion (central), used in B_g g*
    @inline function interfaceFluxAuxiliary!(_surfaceFlux::Vector{Float64}, G::Vector{Float64}, _strideNode_g::Int64, _strideCell_g::Int64, _nCells::Int64)
        @inbounds for Cell in 2:_nCells
            _surfaceFlux[Cell] = 0.5 * (G[1 + (Cell-1) * _strideCell_g - _strideNode_g] + G[1 + (Cell-1) * _strideCell_g])
        end
        _surfaceFlux[1] = 0.5 * (G[1] + G[1])
        _surfaceFlux[_nCells + 1] = 0.5 * (G[1 + _nCells * _strideCell_g - _strideNode_g] + G[1 + _nCells * _strideCell_g - _strideNode_g])
        nothing
    end

    # ------------------------------------------------------------------------
    # Lift face fluxes into interiors and accumulate in Dc:
    # Dc += -(2/Δρ) M_ρ^{-1} [ Lift(v c*) + Lift(g* with left/right scaling) ].
    @inline function surfaceIntegral!(Dc, y, idx, _strideNode::Int64, _strideCell::Int64, _nNodes::Int64, _nCells::Int64, _deltarho, _invMM::Matrix{Float64}, invMrhoM::Vector{<:AbstractMatrix}, _polyDeg, _invWeights, _exactInt, c_star, g_star; left_scale= -1.0, right_scale= +1.0, v=0.0, _h = nothing, tmp = nothing, faces_vc::Union{Nothing,Vector{Float64}}=nothing)
        if _h === nothing
            _h = similar(Dc)
        end
        if tmp === nothing
            tmp = similar(Dc)
        end

        fill!(_h, 0.0)  # scratch buffer for lifted face contributions per node

        if faces_vc === nothing
            # build v*c* on faces using scalar v (constant)
            local_fvc = Vector{Float64}(undef, _nCells+1)
            @inbounds for i in 1:(_nCells+1)
                local_fvc[i] = v * c_star[i]
            end
            surfaceIntegraly!(_h, y, idx, _strideNode, _strideCell, 1, _nNodes, local_fvc, _nCells, _nNodes, _invMM, _polyDeg, _invWeights, _exactInt, use_star_only=true)
            surfaceIntegraly!(_h, y, idx, _strideNode, _strideCell, 1, _nNodes, g_star, _nCells, _nNodes, _invMM, _polyDeg, _invWeights, _exactInt, left_scale=left_scale, right_scale=right_scale, use_star_only=true)
        else
            # faces_vc already equals v_face * c*
            surfaceIntegraly!(_h, y, idx, _strideNode, _strideCell, 1, _nNodes, faces_vc, _nCells, _nNodes, _invMM, _polyDeg, _invWeights, _exactInt, use_star_only=true)
            surfaceIntegraly!(_h, y, idx, _strideNode, _strideCell, 1, _nNodes, g_star, _nCells, _nNodes, _invMM, _polyDeg, _invWeights, _exactInt, left_scale=left_scale, right_scale=right_scale, use_star_only=true)
        end
        @inbounds for Cell in 1:_nCells
            h_cell  = @view _h[1 + (Cell-1)*_nNodes : Cell*_nNodes]
            Dh_cell = @view Dc[1 + (Cell-1)*_nNodes : Cell*_nNodes]
            ldiv!(tmp, _invMM, h_cell)
            mul!(h_cell, invMrhoM[Cell], tmp)
            @. Dh_cell += -(2.0 / _deltarho) * h_cell
        end

        nothing
    end

    # Exact-integration lifting of surface contributions
    # uses inverse-mass first/last columns to spread face fluxes to all nodes (SBP lifting).
    @inline function surfaceIntegraly!(stateDer,state,idx, strideNode_state, strideCell_state, strideNode_stateDer, strideCell_stateDer,_surfaceFlux,_nCells,_nNodes,_invMM, _polyDeg,_invWeights,_exactInt::exact_integration; left_scale=1.0, right_scale=1.0, use_star_only::Bool=false)
        @inbounds for Cell in 1:_nCells
            lscale = (left_scale isa AbstractVector)  ? left_scale[Cell]  : left_scale
            rscale = (right_scale isa AbstractVector) ? right_scale[Cell] : right_scale
            @inbounds for Node in 1:_nNodes
                if !use_star_only
                    stateDer[1 + ( (Cell-1) * strideCell_stateDer) + ( (Node-1) * strideNode_stateDer)] -= (_invMM[Node, 1]) * ((state[idx[1] + ((Cell-1) * strideCell_state)]) - (lscale*_surfaceFlux[Cell])) - (_invMM[Node, end]) * ((state[idx[1] + ((Cell-1) * strideCell_state) + (_polyDeg * strideNode_state)]) - (rscale*_surfaceFlux[Cell+1]))  # (Left/Right face contribution via inverse-mass first/last column)
                else
                    stateDer[1 + ( (Cell-1) * strideCell_stateDer) + ( (Node-1) * strideNode_stateDer)] -= (_invMM[Node, 1]) * (- lscale*_surfaceFlux[Cell]) - (_invMM[Node, end]) * (- rscale*_surfaceFlux[Cell+1])  # (Left/Right face contribution via inverse-mass first/last column)
                end
            end
        end
        nothing
    end
end