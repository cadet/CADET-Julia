using LinearAlgebra
using Printf
using Plots
using CADETJulia

# --------------------------- utilities ---------------------------------------
# map DG nodes to physical radius in a cell
hatrho(nodes, ρi, Δρ) = @. ρi + (nodes + 1) * (Δρ/2)

# build face velocities for constant volumetric flow (v ∝ 1/r)
function face_velocities(op; v_in::Float64=1.0, ρ_in::Float64=op.rho_i[1])
    rv = Vector{Float64}(undef, op.nCells + 1)
    rv[1] = v_in
    @inbounds for f in 2:op.nCells
        rv[f] = v_in * (ρ_in / op.rho_i[f])
    end
    rv[op.nCells + 1] = v_in * (ρ_in / op.rho_ip1[op.nCells])
    return rv
end

# build full vector of physical radii at all DG nodes
function radii_vector(op)
    r = similar(op.Dc)
    @inbounds for cell in 1:op.nCells
        ρnodes = hatrho(op.nodes, op.rho_i[cell], op.deltarho)
        r[(cell-1)*op.nNodes+1 : cell*op.nNodes] .= ρnodes
    end
    return r
end

# ---------------- residual assembly for c(r)=r, ----------------
function assemble_residual_c_eq_r(polyDeg::Int, nCells::Int; ρ_in::Float64=1.0, ρ_out::Float64=2.0, v_in::Float64=1.0, d_rad::Float64=0.0, c_in::Float64=ρ_in)
    op = CADETJulia.RadialConvDispOp(polyDeg, nCells, ρ_in, ρ_out)

    # state y: c(r) = r
    y = similar(op.Dc); fill!(y, 0.0)
    @inbounds for cell in 1:op.nCells
        ρnodes = hatrho(op.nodes, op.rho_i[cell], op.deltarho)
        y[(cell-1)*op.nNodes+1 : cell*op.nNodes] .= ρnodes
    end

    radial_v = face_velocities(op; v_in=v_in, ρ_in=ρ_in)
    left_scale_vec  = op.rho_i   .* d_rad
    right_scale_vec = op.rho_ip1 .* d_rad

    # scratch aliases
    Dc     = op.Dc
    fill!(Dc, 0.0)
    c_star = op.c_star
    g_star = op.g_star
    Dg     = op.Dg
    h      = op.h
    mul1   = op.mul1
    idx    = 1:op.nPoints

    cIn = c_in
    CADETJulia.RadialConvDispOperatorDG.radialresidualImpl!(Dc, y, idx, op.strideNode, op.strideCell, op.nNodes, op.nCells, op.deltarho, polyDeg, op.invWeights, op.polyDerM, op.invMM, op.MM00, op.MM01, op.nodes, v_in, d_rad, cIn, c_star, g_star, Dg, h, mul1, op.rho_i, op.rho_ip1, radial_v, left_scale_vec, right_scale_vec)

    return op, y, Dc
end

# ---------------------------- EOC ----------------------------
function residual_err_inf(p::Int, nC::Int; ρ_in::Float64=0.05, ρ_out::Float64=0.10, v_in::Float64=2/60, d_rad::Float64=1e-6, c_in::Float64=ρ_in)
    _, _, Dc = assemble_residual_c_eq_r(p, nC; ρ_in=ρ_in, ρ_out=ρ_out, v_in=v_in, d_rad=d_rad, c_in=c_in)
    return maximum(abs.(Dc))
end

# EOC for halving h (log2 of consecutive error ratios)
function eoc_from_errs(errs::AbstractVector{<:Real})
    if length(errs) < 2
        return Float64[]
    end
    e = similar(collect(float.(errs[1:end-1])))
    @inbounds for k in 1:length(e)
        e[k] = log(errs[k]/errs[k+1]) / log(2)
    end
    return e
end

# h-refinement study: nCells doubles each step
function eoc_h(p::Int; nCells_list = [4,8,16,32], ρ_in=0.05, ρ_out=0.10, v_in=2/60, d_rad=1e-6, c_in=ρ_in)
    errs = [residual_err_inf(p, nC; ρ_in=ρ_in, ρ_out=ρ_out, v_in=v_in, d_rad=d_rad, c_in=c_in) for nC in nCells_list]
    e = eoc_from_errs(errs)
    @info "h-refinement: errs" nCells_list errs
    @info "h-refinement: EOC" e
    return (; nCells_list, errs, eoc=e)
end

# p-refinement study: increase polynomial degree at fixed mesh
function eoc_p(nC::Int; p_list = [1,2,3,4,5], ρ_in=0.05, ρ_out=0.10, v_in=2/60, d_rad=1e-6, c_in=ρ_in)
    errs = [residual_err_inf(p, nC; ρ_in=ρ_in, ρ_out=ρ_out, v_in=v_in, d_rad=d_rad, c_in=c_in) for p in p_list]
    @info "p-refinement: errs" p_list errs
    return (; p_list, errs)
end

# ---------------------------- EOC plotting helpers ----------------------------
function plot_eoc_h(p::Int; nCells_list=[4,8,16,32], ρ_in=0.05, ρ_out=0.10, v_in=2/60, d_rad=1e-6, c_in=ρ_in)
    data = eoc_h(p; nCells_list=nCells_list, ρ_in=ρ_in, ρ_out=ρ_out, v_in=v_in, d_rad=d_rad, c_in=c_in)
    hs = (ρ_out - ρ_in) ./ data.nCells_list
    plt = plot(hs, data.errs; xscale=:log10, yscale=:log10, marker=:circle, xlabel="h = Δρ", ylabel="‖Residual‖∞", title="h-refinement residual (p=$(p), d_rad=$(d_rad), c_in=$(c_in))")
    png = joinpath(@__DIR__, "eoc_h_residual_p$(p).png")
    savefig(plt, png); println("Saved h-refinement plot → ", png)
    display(plt)
    return (; h=hs, errs=data.errs, eoc=data.eoc)
end

function plot_eoc_p(nC::Int; p_list=[1,2,3,4,5], ρ_in=0.05, ρ_out=0.10, v_in=2/60, d_rad=1e-6, c_in=ρ_in)
    data = eoc_p(nC; p_list=p_list, ρ_in=ρ_in, ρ_out=ρ_out, v_in=v_in, d_rad=d_rad, c_in=c_in)
    plt = plot(data.p_list, data.errs; yscale=:log10, marker=:square, xlabel="p (polyDeg)", ylabel="‖Residual‖∞", title="p-refinement residual (nCells=$(nC), d_rad=$(d_rad), c_in=$(c_in))")
    png = joinpath(@__DIR__, "eoc_p_residual_n$(nC).png")
    savefig(plt, png); println("Saved p-refinement plot → ", png)
    display(plt)
    return data
end


# ---------------------------- main (script/REPL) ------------------------------
function __main__()
    eoc_h(3; nCells_list=[16,32,64,128,256], ρ_in=0.05, ρ_out=0.1, v_in=2/60, d_rad=1e-6, c_in=0.05)
    eoc_p(32; p_list=[5,20,40,50], ρ_in=0.05, ρ_out=0.1, v_in=2/60, d_rad=1e-6, c_in=0.05)
    plot_eoc_h(3; nCells_list=[16,32,64,128,256], ρ_in=0.05, ρ_out=0.1, v_in=2/60, d_rad=1e-6, c_in=0.05)
    plot_eoc_p(32; p_list=[5,20,40,50], ρ_in=0.05, ρ_out=0.1, v_in=2/60, d_rad=1e-6, c_in=0.05)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__) || isinteractive()
    __main__()
end