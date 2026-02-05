"""
    radialGetGBlock(cellIdx, _nNodes, _polyDerM, _nCells, _invMM, _deltarho, _invWeights, _polyDeg)
# Returns
The G block matrix of size (nNodes × nNodes+2).
"""
function radialGetGBlock(cellIdx, _nNodes, _polyDerM, _nCells, _invMM, _deltarho, _invWeights, _polyDeg)
    gBlock = zeros(_nNodes, _nNodes + 2)
    @. @views gBlock[:, 2:_nNodes+1] = _polyDerM

    if cellIdx != 1 && cellIdx != _nCells # Inner cell (Eq. A.1)
        @. @views gBlock[:, 1] -= 0.5 * _invMM[:, 1]
        @. @views gBlock[:, 2] += 0.5 * _invMM[:, 1]
        @. @views gBlock[:, _nNodes + 1] -= 0.5 * _invMM[:, _nNodes]
        @. @views gBlock[:, _nNodes + 2] += 0.5 * _invMM[:, _nNodes]
    elseif cellIdx == 1 # Left boundary cell
        if cellIdx == _nCells # Single cell: both boundary conditions, g = (2/Δρ) D c
            return gBlock * 2 / _deltarho
        end
        # Left BC: c+ = c- → no left correction; only right correction
        @. @views gBlock[:, _nNodes + 1] -= 0.5 * _invMM[:, _nNodes]
        @. @views gBlock[:, _nNodes + 2] += 0.5 * _invMM[:, _nNodes]
    elseif cellIdx == _nCells # Right boundary cell
        # Right BC: c+ = c- → no right correction; only left correction
        @. @views gBlock[:, 1] -= 0.5 * _invMM[:, 1]
        @. @views gBlock[:, 2] += 0.5 * _invMM[:, 1]
    elseif cellIdx == 0 || cellIdx == _nCells + 1 # Ghost cells
        gBlock .= 0
    end
    gBlock .*= 2 / _deltarho

    return gBlock
end


"""
    auxBlockGstar(cellIdx, leftG, middleG, rightG, _nNodes, _nCells)
# Returns
The assembled G* block matrix of size (nNodes × 3*nNodes+2).
"""
function radialAuxBlockGstar(cellIdx, leftG, middleG, rightG, _nNodes, _nCells)
    gStarDC = zeros(_nNodes, 3 * _nNodes + 2)

    if cellIdx != 1
        @. @views gStarDC[1, _nNodes+1:_nNodes+2+_nNodes] .+= middleG[1, 1:_nNodes+2]
        @. @views gStarDC[1, 1:_nNodes+2] .+= leftG[_nNodes, 1:_nNodes+2]
    end
    if cellIdx != _nCells
        @. @views gStarDC[_nNodes, _nNodes+1:_nNodes+2+_nNodes] .+= middleG[_nNodes, 1:_nNodes+2]
        @. @views gStarDC[_nNodes, 1+2*_nNodes : 2*_nNodes+_nNodes+2] .+= rightG[1, 1:_nNodes+2]
    end

    gStarDC .*= 0.5

    return gStarDC
end


"""
    radialDGjacobianConvBlock(cellIdx, _nNodes, _polyDerM, _MM00, _invrMM_cell, _deltarho, _nCells)

# Returns
The convection block matrix (without velocity factor v).
"""
function radialDGjacobianConvBlock(cellIdx, _nNodes, _polyDerM, _MM00, _invrMM_cell, _deltarho, _nCells)
    convBlock = zeros(_nNodes, _nNodes + 1)

    # Weak form volume integral: M_ρ^{-1} * D^T * M^{(0,0)}
    DtM00 = transpose(_polyDerM) * _MM00
    @views convBlock[:, 2:_nNodes+1] = _invrMM_cell * DtM00

    # Surface integral (upwind forward flow): c*_L = c_{-1}, c*_R = c_N
    # Left face: M_ρ^{-1}[:,1] * (c_{-1})  →  contributes to j=-1 column
    @views convBlock[:, 1] .+= _invrMM_cell[:, 1]

    # Right face: -M_ρ^{-1}[:,end] * (c_N)  →  contributes to j=N column
    @views convBlock[:, _nNodes+1] .-= _invrMM_cell[:, _nNodes]

    convBlock .*= 2.0 / _deltarho

    return convBlock
end


"""
    radialDGjacobianDispBlock(cellIdx, _nNodes, _polyDerM, _invMM, _invrMM_cell, _rMM_cell, _rho_i_left, _rho_i_right, _deltarho, _invWeights, _polyDeg, _nCells)

# Returns
The dispersion block matrix (without d_rad factor).
"""
function radialDGjacobianDispBlock(cellIdx, _nNodes, _polyDerM, _invMM, _invrMM_cell, _rMM_cell, _rho_i_left, _rho_i_right, _deltarho, _invWeights, _polyDeg, _nCells)

    # 1. Compute G blocks (auxiliary Jacobian, cell-independent)
    gBlock = radialGetGBlock(cellIdx, _nNodes, _polyDerM, _nCells, _invMM, _deltarho, _invWeights, _polyDeg)
    gBlockLeft = radialGetGBlock(cellIdx - 1, _nNodes, _polyDerM, _nCells, _invMM, _deltarho, _invWeights, _polyDeg)
    gBlockRight = radialGetGBlock(cellIdx + 1, _nNodes, _polyDerM, _nCells, _invMM, _deltarho, _invWeights, _polyDeg)

    # 2. Compute G* (numerical flux Jacobian of auxiliary equation)
    gStarDC = radialAuxBlockGstar(cellIdx, gBlockLeft, gBlock, gBlockRight, _nNodes, _nCells)

    # 3. S_g_norm = D^T * M_ρ (dispersion stiffness without d_rad)
    S_g_norm = transpose(_polyDerM) * _rMM_cell

    # 4. Assemble dispersion block
    dispBlock = zeros(_nNodes, 3 * _nNodes + 2)

    # Volume term: -M_ρ^{-1} * S_g_norm * G̅  (on middle columns: nNodes+1 to 2*nNodes+2)
    dispBlock[:, _nNodes + 1 : _nNodes + 2 + _nNodes] .= _invrMM_cell * (-S_g_norm * gBlock)

    # Surface term: M_ρ^{-1} * B_g_norm * G̅*
    # B_g_norm acts on rows of gStarDC:
    #   row 1 (left face):  scaled by -rho_i_left   (from B_g[1,1] = -rho_left)
    #   row N (right face): scaled by +rho_i_right  (from B_g[N,N] = +rho_right)
    # But gStarDC is (nNodes × 3*nNodes+2) with non-zero rows only at 1 and nNodes.
    bg_gstar = zeros(_nNodes, 3 * _nNodes + 2)
    @views bg_gstar[1, :] .= -_rho_i_left .* gStarDC[1, :]
    @views bg_gstar[_nNodes, :] .= _rho_i_right .* gStarDC[_nNodes, :]

    dispBlock .+= _invrMM_cell * bg_gstar

    # Scale by 2/Δρ
    dispBlock .*= 2.0 / _deltarho

    return dispBlock
end


"""
    radialAddConvJacobian!(jacobian, convBlocks, _nCells, _nNodes, v, compstride)
"""
function radialAddConvJacobian!(jacobian, convBlocks, _nCells, _nNodes, v, compstride)
    # First cell: no j=-1 dependency (inlet), use columns 2:end
    @. @views jacobian[1 + compstride : _nNodes + compstride, 1 + compstride : _nNodes + compstride] += v * convBlocks[1][:, 2:end]

    # Remaining cells: use all columns (j=-1, 0, ..., N)
    for cell in 2:_nCells
        @. @views jacobian[1 + (cell-1)*_nNodes + compstride : (cell-1)*_nNodes + _nNodes + compstride, (cell-1)*_nNodes + compstride : (cell-1)*_nNodes + _nNodes + compstride] += v * convBlocks[cell]
    end
    nothing
end


"""
    radialAddDispJacobian!(jacobian, dispBlocks, _nCells, _nNodes, d_rad, compstride)
"""
function radialAddDispJacobian!(jacobian, dispBlocks, _nCells, _nNodes, d_rad, compstride)
    nTotal = _nCells * _nNodes

    for cell in 1:_nCells
        # Row range in global Jacobian
        row_start = 1 + (cell-1)*_nNodes + compstride
        row_end = cell*_nNodes + compstride

        # dispBlock column k maps to global node g = (cell-2)*nNodes + k - 1
        # Valid range: g ∈ [1, nTotal]
        # So k ∈ [2 - (cell-2)*nNodes, nTotal + 1 - (cell-2)*nNodes]
        k_start = max(1, 2 - (cell-2)*_nNodes)
        k_end = min(3*_nNodes+2, nTotal + 1 - (cell-2)*_nNodes)

        if k_start > k_end
            continue
        end

        g_start = (cell-2)*_nNodes + k_start - 1  # global node (1-based)
        g_end = (cell-2)*_nNodes + k_end - 1

        @views jacobian[row_start:row_end, g_start + compstride : g_end + compstride] .+= d_rad .* dispBlocks[cell][:, k_start:k_end]
    end
    nothing
end


"""
    ConvDispJacobian(model::rLRM, v, p)
"""
function ConvDispJacobian(model::rLRM, v, p)
    nNodes = model.ConvDispOpInstance.nNodes
    nCells = model.nCells
    nComp = model.nComp
    deltarho = model.ConvDispOpInstance.deltarho
    polyDerM = model.ConvDispOpInstance.polyDerM
    invMM = model.ConvDispOpInstance.invMM
    MM00 = model.ConvDispOpInstance.MM00
    invrMM = model.ConvDispOpInstance.invrMM
    rMM = model.ConvDispOpInstance.rMM
    rho_i = model.ConvDispOpInstance.rho_i
    invWeights = model.ConvDispOpInstance.invWeights
    polyDeg = model.polyDeg

    ConvDispJac = zeros(nComp * model.bindStride + model.adsStride, nComp * model.bindStride + model.adsStride)

    # Precompute per-cell convection blocks
    convBlocks = Vector{Matrix{Float64}}(undef, nCells)
    for cell in 1:nCells
        convBlocks[cell] = radialDGjacobianConvBlock(cell, nNodes, polyDerM, MM00, invrMM[cell], deltarho, nCells)
    end

    # Precompute per-cell dispersion blocks (without d_rad factor)
    dispBlocks = Vector{Matrix{Float64}}(undef, nCells)
    for cell in 1:nCells
        dispBlocks[cell] = radialDGjacobianDispBlock(cell, nNodes, polyDerM, invMM, invrMM[cell], rMM[cell], rho_i[cell], rho_i[cell+1], deltarho, invWeights, polyDeg, nCells)
    end

    # Inlet Danckwerts BC correction: g*[1] = v/d_rad * (c_0 - c_in)
    # Surface integral gives -(2/Δρ) * invrMM[1][:,1] * rho_i[1] * d_rad * (v/d_rad) * c_0
    # The d_rad cancels, yielding a v-scaled boundary correction on column 1 of cell 1.
    inlet_corr = -(2.0 / deltarho) * rho_i[1] .* @view(invrMM[1][:, 1])

    # Assemble for all components
    for i in 1:nComp
        compstride = (i-1) * nNodes * nCells

        # Dispersion (scaled by per-component d_rad)
        d_rad_i = isa(model.d_rad, Vector) ? model.d_rad[i] : model.d_rad
        if isa(d_rad_i, Function)
            # For variable d_rad, use a representative value (first interface)
            d_rad_i = d_rad_i(rho_i[1])
        end
        radialAddDispJacobian!(ConvDispJac, dispBlocks, nCells, nNodes, d_rad_i, compstride)

        # Convection (scaled by velocity v)
        radialAddConvJacobian!(ConvDispJac, convBlocks, nCells, nNodes, v, compstride)

        # Inlet g* boundary correction (scales with velocity, not d_rad)
        @views ConvDispJac[1+compstride:nNodes+compstride, 1+compstride] .+= v .* inlet_corr
    end

    return ConvDispJac
end


"""
    ConvDispJacobian(model::rLRMP, v, p)
Film diffusion Jacobian blocks:
    ∂RHS_c/∂c   += -Fc * (3/Rp) * M_ρ^{-1} * M_K  (diagonal block)
    ∂RHS_c/∂c_p += +Fc * (3/Rp) * M_ρ^{-1} * M_K  (off-diagonal block)
    ∂RHS_cp/∂c  += +(3/Rp/eps_p) * M_ρ^{-1} * M_K  (off-diagonal block)
    ∂RHS_cp/∂c_p += -(3/Rp/eps_p) * M_ρ^{-1} * M_K  (diagonal block)
"""
function ConvDispJacobian(model::rLRMP, v, p)
    nNodes = model.ConvDispOpInstance.nNodes
    nCells = model.nCells
    nComp = model.nComp
    deltarho = model.ConvDispOpInstance.deltarho
    polyDerM = model.ConvDispOpInstance.polyDerM
    invMM = model.ConvDispOpInstance.invMM
    MM00 = model.ConvDispOpInstance.MM00
    invrMM = model.ConvDispOpInstance.invrMM
    rMM = model.ConvDispOpInstance.rMM
    rho_i = model.ConvDispOpInstance.rho_i
    invWeights = model.ConvDispOpInstance.invWeights
    polyDeg = model.polyDeg

    ConvDispJac = zeros(nComp * model.bindStride + model.adsStride, nComp * model.bindStride + model.adsStride)

    # Precompute per-cell convection blocks
    convBlocks = Vector{Matrix{Float64}}(undef, nCells)
    for cell in 1:nCells
        convBlocks[cell] = radialDGjacobianConvBlock(cell, nNodes, polyDerM, MM00, invrMM[cell], deltarho, nCells)
    end

    # Precompute per-cell dispersion blocks (without d_rad factor)
    dispBlocks = Vector{Matrix{Float64}}(undef, nCells)
    for cell in 1:nCells
        dispBlocks[cell] = radialDGjacobianDispBlock(cell, nNodes, polyDerM, invMM, invrMM[cell], rMM[cell], rho_i[cell], rho_i[cell+1], deltarho, invWeights, polyDeg, nCells)
    end

    # Precompute per-cell film diffusion coupling: M_ρ^{-1} * M_K
    # This is needed for the film diffusion terms in the Jacobian
    filmDiffBlocks = Vector{Matrix{Float64}}(undef, nCells)
    for cell in 1:nCells
        filmDiffBlocks[cell] = invrMM[cell] * model.FilmDiffOpInstance.M_K[cell]
    end

    jacobiStride = nNodes * nCells

    # Inlet Danckwerts BC correction (same as rLRM)
    inlet_corr = - (2.0 / deltarho) * rho_i[1] .* @view(invrMM[1][:, 1])

    # Assemble for all components
    for i in 1:nComp
        compstride = (i-1) * jacobiStride

        # Dispersion (scaled by per-component d_rad)
        d_rad_i = isa(model.d_rad, Vector) ? model.d_rad[i] : model.d_rad
        if isa(d_rad_i, Function)
            d_rad_i = d_rad_i(rho_i[1])
        end
        radialAddDispJacobian!(ConvDispJac, dispBlocks, nCells, nNodes, d_rad_i, compstride)

        # Convection (scaled by velocity v)
        radialAddConvJacobian!(ConvDispJac, convBlocks, nCells, nNodes, v, compstride)

        # Inlet g* boundary correction (scales with velocity, not d_rad)
        @views ConvDispJac[1+compstride:nNodes+compstride, 1+compstride] .+= v .* inlet_corr

        # Film diffusion coupling terms (per-cell assembly)
        # kf coefficient for this component
        kf_i = isa(model.kf, Vector) ? model.kf[i] : model.kf
        if isa(kf_i, Function)
            kf_i = kf_i(rho_i[1])
        end

        # The film diffusion coefficient Q = 3/Rp
        Q = model.FilmDiffOpInstance.Q  # = 3/Rp

        for cell in 1:nCells
            cell_start = 1 + (cell-1)*nNodes + compstride
            cell_end = cell*nNodes + compstride
            pore_start = 1 + jacobiStride*nComp + (cell-1)*nNodes + compstride
            pore_end = jacobiStride*nComp + cell*nNodes + compstride

            # For the kf scaling: M_K is pre-computed with kf embedded, but if kf differs per component,
            # we need to scale. Since FilmDiffOpInstance uses kf[1], scale by kf_i/kf_base.
            kf_base = isa(model.kf, Vector) ? model.kf[1] : model.kf
            if isa(kf_base, Function)
                kf_base = kf_base(rho_i[1])
            end
            kf_scale = kf_i / kf_base

            # Film diffusion block for this cell: Q * kf_scale * filmDiffBlocks[cell]
            filmBlock = Q * kf_scale * filmDiffBlocks[cell]

            # ∂(dc/dt)/∂c += -Fc * filmBlock
            @views ConvDispJac[cell_start:cell_end, cell_start:cell_end] .-= model.Fc .* filmBlock

            # ∂(dc/dt)/∂c_p += +Fc * filmBlock
            @views ConvDispJac[cell_start:cell_end, pore_start:pore_end] .+= model.Fc .* filmBlock

            # ∂(dcp/dt)/∂c += +(1/eps_p) * filmBlock
            @views ConvDispJac[pore_start:pore_end, cell_start:cell_end] .+= filmBlock ./ model.eps_p

            # ∂(dcp/dt)/∂c_p += -(1/eps_p) * filmBlock
            @views ConvDispJac[pore_start:pore_end, pore_start:pore_end] .-= filmBlock ./ model.eps_p
        end
    end

    return ConvDispJac
end
