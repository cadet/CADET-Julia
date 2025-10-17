using LinearAlgebra
using Statistics
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

# ---------------- residual assembly for c(r)=r, convection-only ---------------
function assemble_residual_c_eq_r(polyDeg::Int, nCells::Int; ρ_in::Float64=1.0, ρ_out::Float64=2.0, v_in::Float64=1.0)
    op = CADETJulia.RadialConvDispOp(polyDeg, nCells, ρ_in, ρ_out)

    # state y: c(r) = r
    y = similar(op.Dc); fill!(y, 0.0)
    @inbounds for cell in 1:op.nCells
        ρnodes = hatrho(op.nodes, op.rho_i[cell], op.deltarho)
        y[(cell-1)*op.nNodes+1 : cell*op.nNodes] .= ρnodes
    end

    d_rad = 0.0
    radial_v = face_velocities(op; v_in=v_in, ρ_in=ρ_in)
    left_scale_vec  = op.rho_i   .* d_rad
    right_scale_vec = op.rho_ip1 .* d_rad

    # scratch aliases
    Dc     = op.Dc;     fill!(Dc, 0.0)
    c_star = op.c_star
    g_star = op.g_star
    Dg     = op.Dg
    h      = op.h
    mul1   = op.mul1
    idx    = 1:op.nPoints

    cIn = ρ_in

    CADETJulia.RadialConvDispOperatorDG.radialresidualImpl!(
        Dc, y, idx,
        op.strideNode, op.strideCell,
        op.nNodes, op.nCells, op.deltarho,
        polyDeg, op.invWeights, op.polyDerM, op.invMM, op.MM00, op.MM01,
        op.nodes, v_in, d_rad, cIn,
        c_star, g_star, Dg, h, mul1,
        op.rho_i, op.rho_ip1,
        radial_v, left_scale_vec, right_scale_vec,
    )

    return op, y, Dc
end

# ---------------------------- comparison + plot -------------------------------
function compare_residual_c_eq_r(polyDeg::Int=3, nCells::Int=16; ρ_in::Float64=1.0, ρ_out::Float64=2.0, v_in::Float64=1.0)
    op, y, Dc = assemble_residual_c_eq_r(polyDeg, nCells; ρ_in=ρ_in, ρ_out=ρ_out, v_in=v_in)
    r = radii_vector(op)

    # analytic residual: (1/r) d_r (r v c) with v = v_in*ρ_in/r, c=r  ⇒ v_in*ρ_in/r
    expected = zeros(length(r))
    err_vec  = Dc .- expected

    # print diagnostics
    @printf("\n=== Residual check for c(r)=r (p=%d, nCells=%d) ===\n", polyDeg, nCells)
    @printf("‖Dc‖∞      = %.6e\n", maximum(abs.(Dc)))
    @printf("‖expected‖∞= %.6e\n", maximum(abs.(expected)))
    @printf("‖error‖∞   = %.6e\n", maximum(abs.(err_vec)))
    @printf("‖error‖₂    = %.6e\n", norm(err_vec))
    @printf("mean(error)= %.6e\n", mean(err_vec))
    println("Dc[1:8]       = ", collect(Dc[1:min(end,8)]))
    println("expected[1:8] = ", collect(expected[1:min(end,8)]))
    println("error[1:8]    = ", collect(err_vec[1:min(end,8)]))

    # plot residual Dc and analytic (expected) only
    p = plot(xlabel="r", ylabel="value", title="Residual Dc vs expected (vc volume term)")
    plot!(p, r, Dc; label="Residual Dc")
    plot!(p, r, expected; label="Analytic v_in*ρ_in/r", linestyle=:dash)

    p2 = plot(r, err_vec; label="error = Dc - analytic", xlabel="r", ylabel="error",
              title="Residual error for c(r)=r")

    # save
    png1 = joinpath(@__DIR__, "compare_residual_p$(polyDeg)_n$(nCells).png")
    png2 = joinpath(@__DIR__, "compare_residual_error_p$(polyDeg)_n$(nCells).png")
    savefig(p, png1);  println("Saved plot → ", png1)
    savefig(p2, png2); println("Saved plot → ", png2)

    display(p); display(p2)
    return (; r, Dc, expected, err_vec)
end

# ---------------------------- main (script/REPL) ------------------------------
function __main__()
    compare_residual_c_eq_r(3, 16; ρ_in=1.0, ρ_out=2.0, v_in=1.0)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__) || isinteractive()
    __main__()
end