using CADETJulia
using CADETJulia: DGElements, RadialConvDispOperatorDG
using Plots

function mms_compare_plot(; p=4, nCells=8, ri=0.01, ro=0.05, v_in=2/60, D=1.0e-4)
    # Manufactured solution and derivatives
    A, B = 1.0, 0.0
    c_fun(r)  = A * (r - ri)^2 * (ro - r)^2 + B
    dcdr(r)   = 2A * ((r - ri)*(ro - r)^2 - (r - ri)^2*(ro - r))
    d2cdr2(r) = 2A * ((ro - r)^2 - 4*(r - ri)*(ro - r) + (r - ri)^2)

    v_of_r(r) = v_in * (ri / r)
    F_fun(r) = -(v_of_r(r)/r) * dcdr(r) + D*( d2cdr2(r) + (1/r)*dcdr(r) )

    # DG geometry / matrices
    nodes, invw = DGElements.lglnodes(p)
    nNodes  = p + 1
    deltarho = (ro - ri) / nCells
    rho_i  = [ri + (k-1)*deltarho for k in 1:nCells]
    rho_ip1 = [ri + k*deltarho for k in 1:nCells]
    nPoints = nNodes * nCells

    invWeights = invw
    invMM      = DGElements.invMMatrix(nodes, p)
    MM         = DGElements.MMatrix(nodes, p)
    polyDerM   = DGElements.derivativeMatrix(p, nodes)
    invMrhoM   = [DGElements.invMrhoMatrix(nodes, p, deltarho, rho_i[c]) for c in 1:nCells]
    SgMatrix   = [DGElements.weighted_stiff_Matrix(nodes, p, rho_i[c], deltarho, _ -> D) for c in 1:nCells]

    # map reference nodes to radii (cell-by-cell, LGL order)
    r_all = Float64[]
    for k in 1:nCells
        ρ0 = rho_i[k]
        for ξ in nodes
            push!(r_all, ρ0 + (ξ + 1.0) * (deltarho/2))
        end
    end

    # buffers
    c_star = zeros(nCells+1)
    g_star = zeros(nCells+1)
    Dc_num = zeros(nPoints)
    Dg     = zeros(nPoints)
    hbuf   = zeros(nPoints)
    mul1   = zeros(nNodes)

    y = [c_fun(r) for r in r_all]
    cp_dummy = copy(y) # film term disabled below

    F_anal = [F_fun(r) for r in r_all]
    idx = 1:nPoints

    # Face velocities (RFC: v(r) = v_in * ri / r at faces)
    faces_v = RadialConvDispOperatorDG.compute_faces_v(v_in, rho_i, rho_ip1, ri, nCells)

    # Left/right diffusion scale factors ρ D at faces (for each face)
    D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(D, rho_i, rho_ip1)
    left_scale_vec  = @. rho_i  * D_left
    right_scale_vec = @. rho_ip1 * D_right

    # Call DG residual builder (Neumann-0 at outlet, central advective flux)
    RadialConvDispOperatorDG.radialresidualImpl!(
        Dc_num, y, idx,
        1, nNodes,
        nNodes, nCells, deltarho, p,
        invWeights, nodes, polyDerM, invMM, MM, invMrhoM,
        SgMatrix,
        v_in, D,
        c_fun(ri),
        c_star, g_star, Dg, hbuf, mul1,
        RadialConvDispOperatorDG.exact_integration(),
        :neumann0, 0.0,
        rho_i, rho_ip1,
        faces_v,
        left_scale_vec,
        right_scale_vec
    )

    abs_err = maximum(abs.(Dc_num .- F_anal))
    rel_err = abs_err / (maximum(abs.(F_anal)) + eps())

    # helper: test tolerance as in the unit tests
    tol_for_p(p::Int) = p == 1 ? 0.05 : p == 2 ? 0.05 : p == 3 ? 0.025 : p == 4 ? 1/60 : p == 5 ? 0.0125 : 0.01
    tol = tol_for_p(p)
    verdict = rel_err < tol ? "PASS (good)" : "FAIL (bad)"
    @info "MMS comparison" abs_err rel_err tol verdict

    # === Make plots ===
    plt_num_vs_anal = plot(r_all, F_anal, lw=2, label="Analytical RHS",
                           xlabel="r [m]", ylabel="RHS",
                           title="Radial MMS: numerical vs analytical (p=$(p), nCells=$(nCells))")
    plot!(plt_num_vs_anal, r_all, Dc_num, lw=2, ls=:dash, label="DG residual (numerical)")

    plt_error = plot(r_all, Dc_num .- F_anal, lw=2, label="Pointwise error",
                     xlabel="r [m]", ylabel="Error",
                     title="Pointwise error (p=$(p), nCells=$(nCells))")

    plt_combined = plot(plt_num_vs_anal, plt_error, layout=(2,1), size=(900, 800),
                        title="Radial MMS comparison & error (rel_err=$(round(rel_err, sigdigits=3)), tol=$(tol): $(verdict))")

    # Display individually and combined
    display(plt_num_vs_anal)
    display(plt_error)
    display(plt_combined)

    # Save all three
    savefig(plt_num_vs_anal, "radial_mms_num_vs_anal_p$(p)_n$(nCells).png")
    savefig(plt_error,       "radial_mms_error_p$(p)_n$(nCells).png")
    savefig(plt_combined,    "radial_mms_combined_p$(p)_n$(nCells).png")

    return (r=r_all, F=F_anal, DG=Dc_num, abs_err=abs_err, rel_err=rel_err, tol=tol, verdict=verdict)
end

# Example:
res = mms_compare_plot(p=4, nCells=8)
@info "Rel error verdict" rel_err=res.rel_err tol=res.tol verdict=res.verdict