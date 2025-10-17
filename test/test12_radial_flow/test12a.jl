using Test, CADETJulia
using .CADETJulia.RadialConvDispOperatorDG
using LinearAlgebra
using Plots

function run_radial_mms_test(; polyDeg=4, nCells=8, rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4)
    # @info "Running radial DG test" polyDeg nCells rin rout v_in D

    # --- Build operator using the same internals as the solver ---
    Rad = CADETJulia.RadialConvDispOp(polyDeg, nCells, rin, rout; d_rad_const=D)
    nNodes  = Rad.nNodes
    nPoints = Rad.nPoints

    nodes      = Rad.nodes
    polyDerM   = Rad.polyDerM
    invMM      = Rad.invMM
    MM         = Rad.MM
    invWeights = Rad.invWeights
    invMrhoM   = Rad.invMrhoM
    SgMatrix   = Rad.SgMatrix
    deltarho   = Rad.deltarho
    rho_i      = Rad.rho_i
    rho_ip1    = Rad.rho_ip1

    # Face data (ALWAYS upwind for advection), diffusion face scaling = ρ·D
    faces_v = RadialConvDispOperatorDG.compute_faces_v(v_in, rho_i, rho_ip1, rin, nCells)
    D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(D, rho_i, rho_ip1)
    left_scale_vec  = @. rho_i  * D_left
    right_scale_vec = @. rho_ip1 * D_right

    # --- Manufactured polynomial c(r) with ∂c/∂r = 0 at both ends ---
    A, B = 1.0, 0.0
    c_fun(r)  = A * (r - rin)^2 * (rout - r)^2 + B
    dcdr(r)   = 2*A * ((r - rin)*(rout - r)^2 - (r - rin)^2*(rout - r))
    d2cdr2(r) = 2*A * ((rout - r)^2 - 4*(r - rin)*(rout - r) + (r - rin)^2)

    v_of_r(r) = v_in * (rin / r)
    # Conservative advection in cylindrical coords: -(1/r)∂(r v c)/∂r.
    # For RFC v(r)=v_in*rin/r this reduces to -v(r) * dc/dr.
    Lc(r) = -v_of_r(r) * dcdr(r) + D * (d2cdr2(r) + (1/r) * dcdr(r))

    # --- State on DG nodes ---
    y = zeros(Float64, nPoints)
    radii = similar(y)
    for cell in 1:nCells
        rL = rin + (cell-1) * deltarho
        for n in 1:nNodes
            r = rL + (nodes[n] + 1.0) * (deltarho/2)
            idx = (cell-1)*nNodes + n
            radii[idx] = r
            y[idx] = c_fun(r)
        end
    end

    # --- Apply DG operator ---
    Dc = zeros(Float64, nPoints)
    c_star = zeros(Float64, nCells + 1)
    g_star = zeros(Float64, nCells + 1)
    Dg = zeros(Float64, nPoints)
    h = zeros(Float64, nPoints)
    mul1 = zeros(Float64, nNodes)

    cIn = c_fun(rin)

    RadialConvDispOperatorDG.radialresidualImpl!(Dc, y, 1:nPoints, 1, nNodes, nNodes, nCells, deltarho, polyDeg, invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix, v_in, D, cIn, c_star, g_star, Dg, h, mul1, RadialConvDispOperatorDG.exact_integration(), rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec)

    # Analytic RHS on the same nodes
    L_exact = @. Lc(radii)

    # Robust relative L2 (guard for tiny reference):
    num = norm(Dc - L_exact)
    den = max(norm(L_exact), 1.0)   # scale to 1 if reference is tiny, to avoid bogus huge rel_err
    rel_err = num / den
    abs_err = maximum(abs.(Dc - L_exact))
    # @info "comparison" abs_err rel_err tol = 1e-2 verdict = (rel_err < 1e-2 ? "PASS" : "FAIL")

    @test rel_err < 1e-2
    #println("Test passed. rel. L2 error = ", rel_err, ", max|err| = ", abs_err)
    return rel_err
end

# Run on include
#run_radial_mms_test()

# --- Helper: compute MMS error without printing ---
function _radial_mms_error(polyDeg, nCells; rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4)
    Rad = CADETJulia.RadialConvDispOp(polyDeg, nCells, rin, rout; d_rad_const=D)
    nNodes  = Rad.nNodes
    nPoints = Rad.nPoints

    nodes      = Rad.nodes
    polyDerM   = Rad.polyDerM
    invMM      = Rad.invMM
    MM         = Rad.MM
    invWeights = Rad.invWeights
    invMrhoM   = Rad.invMrhoM
    SgMatrix   = Rad.SgMatrix
    deltarho   = Rad.deltarho
    rho_i      = Rad.rho_i
    rho_ip1    = Rad.rho_ip1

    # Face data (RFC profile and ρ·D face scaling)
    faces_v = RadialConvDispOperatorDG.compute_faces_v(v_in, rho_i, rho_ip1, rin, nCells)
    D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(D, rho_i, rho_ip1)
    left_scale_vec  = @. rho_i  * D_left
    right_scale_vec = @. rho_ip1 * D_right

    # Manufactured solution with Neumann(0) at both ends
    A, B = 1.0, 0.0
    c_fun(r)  = A * (r - rin)^2 * (rout - r)^2 + B
    dcdr(r)   = 2*A * ((r - rin)*(rout - r)^2 - (r - rin)^2*(rout - r))
    d2cdr2(r) = 2*A * ((rout - r)^2 - 4*(r - rin)*(rout - r) + (r - rin)^2)

    v_of_r(r) = v_in * (rin / r)
    Lc(r) = -v_of_r(r) * dcdr(r) + D * (d2cdr2(r) + (1/r) * dcdr(r))

    # State values on DG nodes + their radii
    y = zeros(Float64, nPoints)
    radii = similar(y)
    for cell in 1:nCells
        rL = rin + (cell-1) * deltarho
        for n in 1:nNodes
            r = rL + (nodes[n] + 1.0) * (deltarho/2)
            idx = (cell-1)*nNodes + n
            radii[idx] = r
            y[idx] = c_fun(r)
        end
    end

    # Apply DG operator
    Dc = zeros(Float64, nPoints)
    c_star = zeros(Float64, nCells + 1)
    g_star = zeros(Float64, nCells + 1)
    Dg = zeros(Float64, nPoints)
    h = zeros(Float64, nPoints)
    mul1 = zeros(Float64, nNodes)

    cIn = c_fun(rin)

    RadialConvDispOperatorDG.radialresidualImpl!(Dc, y, 1:nPoints, 1, nNodes, nNodes, nCells, deltarho, polyDeg, invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix, v_in, D, cIn, c_star, g_star, Dg, h, mul1, RadialConvDispOperatorDG.exact_integration(), rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec)

    L_exact = @. Lc(radii)
    num = norm(Dc - L_exact)
    den = max(norm(L_exact), 1.0)
    rel_err = num / den
    abs_err = maximum(abs.(Dc - L_exact))
    return rel_err, abs_err
end

# --- Utility: least-squares slope on log-log ---
function _loglog_slope(x::AbstractVector, y::AbstractVector)
    lx = log.(x); ly = log.(y)
    n = length(x)
    sx = sum(lx); sy = sum(ly)
    sxx = sum(lx .* lx); sxy = sum(lx .* ly)
    (n*sxy - sx*sy) / (n*sxx - sx*sx)
end

# --- p-refinement EOC: vary polynomial degree at fixed mesh ---
function run_radial_mms_eoc_p(; p_list = 1:7, nCells=8, rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4)
    #println("\nEOC (p-refinement) — MMS for radial DG (nCells=$(nCells))")
    errs = Float64[]; dof_scale = Float64[]  # use (p+1) as resolution metric
    #@printf("%6s  %8s  %12s  %12s\n", "p", "nNodes", "rel_L2_err", "abs_max_err")
    #@printf("%s\n", repeat('-', 46))
    for p in p_list
        rel, ab = _radial_mms_error(p, nCells; rin=rin, rout=rout, v_in=v_in, D=D)
        push!(errs, rel); push!(dof_scale, p+1)
        #@printf("%6d  %8d  %12.4e  %12.4e\n", p, (p+1)*nCells, rel, ab)
    end
    slope = _loglog_slope(dof_scale, errs)
    #@printf("Estimated slope on log(err) vs log(p+1): %.3f\n", slope)
    return (; p_list = collect(p_list), errs, slope)
end

# --- h-refinement EOC: vary nCells at fixed polynomial degree ---
function run_radial_mms_eoc_h(; nCells_list = [4, 6, 8, 12, 16, 24], polyDeg=4, rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4)
    #println("\nEOC (h-refinement) — MMS for radial DG (polyDeg=$(polyDeg))")
    errs = Float64[]; hvals = Float64[]  # characteristic cell size Δr
    #@printf("%8s  %8s  %12s  %12s\n", "nCells", "(p+1)N", "rel_L2_err", "abs_max_err")
    #@printf("%s\n", repeat('-', 48))
    for nc in nCells_list
        rel, ab = _radial_mms_error(polyDeg, nc; rin=rin, rout=rout, v_in=v_in, D=D)
        push!(errs, rel); push!(hvals, (rout-rin)/nc)
        #@printf("%8d  %8d  %12.4e  %12.4e\n", nc, (polyDeg+1)*nc, rel, ab)
    end
    slope = _loglog_slope(1.0 ./ hvals, errs)  # slope wrt resolution ~ 1/h
    #@printf("Estimated slope on log(err) vs log(1/h): %.3f\n", slope)
    return (; nCells_list = collect(nCells_list), errs, slope)
end

# --- Snapshot for plotting: radii r, numeric residual Dc, analytic RHS L_exact ---
function _radial_mms_snapshot(polyDeg, nCells; rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4)
    Rad = CADETJulia.RadialConvDispOp(polyDeg, nCells, rin, rout; d_rad_const=D)
    nNodes  = Rad.nNodes
    nPoints = Rad.nPoints

    nodes      = Rad.nodes
    polyDerM   = Rad.polyDerM
    invMM      = Rad.invMM
    MM         = Rad.MM
    invWeights = Rad.invWeights
    invMrhoM   = Rad.invMrhoM
    SgMatrix   = Rad.SgMatrix
    deltarho   = Rad.deltarho
    rho_i      = Rad.rho_i
    rho_ip1    = Rad.rho_ip1

    faces_v = RadialConvDispOperatorDG.compute_faces_v(v_in, rho_i, rho_ip1, rin, nCells)
    D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(D, rho_i, rho_ip1)
    left_scale_vec  = @. rho_i  * D_left
    right_scale_vec = @. rho_ip1 * D_right

    # Manufactured solution (same as MMS test)
    A, B = 1.0, 0.0
    c_fun(r)  = A * (r - rin)^2 * (rout - r)^2 + B
    dcdr(r)   = 2*A * ((r - rin)*(rout - r)^2 - (r - rin)^2*(rout - r))
    d2cdr2(r) = 2*A * ((rout - r)^2 - 4*(r - rin)*(rout - r) + (r - rin)^2)
    v_of_r(r) = v_in * (rin / r)
    Lc(r) = -v_of_r(r) * dcdr(r) + D * (d2cdr2(r) + (1/r) * dcdr(r))

    # Discretize manufactured c on DG nodes
    y = zeros(Float64, nPoints)
    radii = similar(y)
    for cell in 1:nCells
        rL = rin + (cell-1) * deltarho
        for n in 1:nNodes
            r = rL + (nodes[n] + 1.0) * (deltarho/2)
            idx = (cell-1)*nNodes + n
            radii[idx] = r
            y[idx] = c_fun(r)
        end
    end

    # Apply DG operator to get numeric residual Dc
    Dc = zeros(Float64, nPoints)
    c_star = zeros(Float64, nCells + 1)
    g_star = zeros(Float64, nCells + 1)
    Dg = zeros(Float64, nPoints)
    h = zeros(Float64, nPoints)
    mul1 = zeros(Float64, nNodes)
    cIn = c_fun(rin)

    RadialConvDispOperatorDG.radialresidualImpl!(Dc, y, 1:nPoints, 1, nNodes, nNodes, nCells, deltarho, polyDeg, invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix, v_in, D, cIn, c_star, g_star, Dg, h, mul1, RadialConvDispOperatorDG.exact_integration(), rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec)

    L_exact = @. Lc(radii)
    return radii, Dc, L_exact
end

# --- Transient MMS test (space+time): u_t = L_d(u) + S(t), exact u(r,t) = e^{α t} φ(r)
function run_radial_mms_time_test(; polyDeg=4, nCells=8, rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4,
                                   α=1.0, T=0.1, reltol=1e-8, abstol=1e-10)
    # Build operator cache
    Rad = CADETJulia.RadialConvDispOp(polyDeg, nCells, rin, rout; d_rad_const=D)
    nNodes  = Rad.nNodes
    nPoints = Rad.nPoints
    nodes      = Rad.nodes
    polyDerM   = Rad.polyDerM
    invMM      = Rad.invMM
    MM         = Rad.MM
    invWeights = Rad.invWeights
    invMrhoM   = Rad.invMrhoM
    SgMatrix   = Rad.SgMatrix
    deltarho   = Rad.deltarho
    rho_i      = Rad.rho_i
    rho_ip1    = Rad.rho_ip1

    faces_v = RadialConvDispOperatorDG.compute_faces_v(v_in, rho_i, rho_ip1, rin, nCells)
    D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(D, rho_i, rho_ip1)
    left_scale_vec  = @. rho_i  * D_left
    right_scale_vec = @. rho_ip1 * D_right

    # Manufactured spatial profile φ(r)
    A, B = 1.0, 0.0
    φ(r)  = A * (r - rin)^2 * (rout - r)^2 + B

    # Node radii and initial condition y0 = φ(r)
    radii = zeros(Float64, nPoints)
    y0 = similar(radii)
    for cell in 1:nCells
        rL = rin + (cell-1) * deltarho
        for n in 1:nNodes
            r = rL + (nodes[n] + 1.0) * (deltarho/2)
            idx = (cell-1)*nNodes + n
            radii[idx] = r
            y0[idx] = φ(r)
        end
    end

    # Precompute discrete L_d(φ) to build source ensuring exact solution u_exact = e^{α t} φ(r)
    Lφ = zeros(Float64, nPoints)
    c_star = zeros(Float64, nCells + 1)
    g_star = zeros(Float64, nCells + 1)
    Dg = zeros(Float64, nPoints)
    h = zeros(Float64, nPoints)
    mul1 = zeros(Float64, nNodes)
    cIn = φ(rin)
    RadialConvDispOperatorDG.radialresidualImpl!(Lφ, y0, 1:nPoints, 1, nNodes, nNodes, nCells, deltarho, polyDeg, invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix, v_in, D, cIn, c_star, g_star, Dg, h, mul1, RadialConvDispOperatorDG.exact_integration(), rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec)

    # Source term S0 = α φ − L_d(φ); full S(t) = e^{α t} * S0
    S0 = α .* y0 .- Lφ

    # ODE: du/dt = L_d(u) + e^{α t} * S0
    function rhs!(du, u, p, t)
        # du ← L_d(u)
        RadialConvDispOperatorDG.radialresidualImpl!(du, u, 1:nPoints, 1, nNodes, nNodes, nCells, deltarho, polyDeg, invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix, v_in, D, cIn, c_star, g_star, Dg, h, mul1, RadialConvDispOperatorDG.exact_integration(), rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec)
        # add source e^{α t} * S0
        @. du += exp(α * t) * S0
        return nothing
    end

    prob = ODEProblem(rhs!, y0, (0.0, T))
    sol = solve(prob, QNDF(autodiff=AutoFiniteDiff()); reltol=reltol, abstol=abstol, save_everystep=false)

    uT_exact = @. exp(α * T) * y0
    num = norm(sol.u[end] .- uT_exact)
    den = max(norm(uT_exact), 1.0)
    rel_err = num / den

    @test rel_err < 1e-6
    return rel_err
end

run_radial_mms_time_test()

# --- Utility: time series at one representative node (for plotting c vs t) ---
function radial_mms_time_series(; polyDeg=4, nCells=8, rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4,
                                 α=1.0, T=0.5, Nt=200)
    # Build operator cache (same as in run_radial_mms_time_test)
    Rad = CADETJulia.RadialConvDispOp(polyDeg, nCells, rin, rout; d_rad_const=D)
    nNodes  = Rad.nNodes
    nPoints = Rad.nPoints
    nodes      = Rad.nodes
    polyDerM   = Rad.polyDerM
    invMM      = Rad.invMM
    MM         = Rad.MM
    invWeights = Rad.invWeights
    invMrhoM   = Rad.invMrhoM
    SgMatrix   = Rad.SgMatrix
    deltarho   = Rad.deltarho
    rho_i      = Rad.rho_i
    rho_ip1    = Rad.rho_ip1

    faces_v = RadialConvDispOperatorDG.compute_faces_v(v_in, rho_i, rho_ip1, rin, nCells)
    D_left, D_right = RadialConvDispOperatorDG.diff_at_faces(D, rho_i, rho_ip1)
    left_scale_vec  = @. rho_i  * D_left
    right_scale_vec = @. rho_ip1 * D_right

    # Manufactured φ(r)
    φ(r) = (r - rin)^2 * (rout - r)^2

    # Initial condition y0 = φ(r) at DG nodes
    y0 = zeros(Float64, nPoints)
    radii = similar(y0)
    for cell in 1:nCells
        rL = rin + (cell-1) * deltarho
        for n in 1:nNodes
            r = rL + (nodes[n] + 1.0) * (deltarho/2)
            idx = (cell-1)*nNodes + n
            radii[idx] = r
            y0[idx] = φ(r)
        end
    end

    # Discrete L_d(φ)
    Lφ = zeros(Float64, nPoints)
    c_star = zeros(Float64, nCells + 1)
    g_star = zeros(Float64, nCells + 1)
    Dg = zeros(Float64, nPoints)
    h = zeros(Float64, nPoints)
    mul1 = zeros(Float64, nNodes)
    cIn = φ(rin)
    RadialConvDispOperatorDG.radialresidualImpl!(Lφ, y0, 1:nPoints, 1, nNodes, nNodes, nCells, deltarho, polyDeg, invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix, v_in, D, cIn, c_star, g_star, Dg, h, mul1, RadialConvDispOperatorDG.exact_integration(), rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec)

    # Source S(t) = e^{α t} * (α φ − L_d(φ))
    S0 = α .* y0 .- Lφ

    function rhs!(du, u, p, t)
        RadialConvDispOperatorDG.radialresidualImpl!(du, u, 1:nPoints, 1, nNodes, nNodes, nCells, deltarho, polyDeg, invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix, v_in, D, cIn, c_star, g_star, Dg, h, mul1, RadialConvDispOperatorDG.exact_integration(), rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec)
        @. du += exp(α * t) * S0
        return nothing
    end

    tspan = (0.0, T)
    saveat = range(tspan[1], tspan[2], length=Nt)
    prob = ODEProblem(rhs!, y0, tspan)
    sol = solve(prob, QNDF(autodiff=AutoFiniteDiff()); reltol=1e-9, abstol=1e-11, saveat=saveat)

    # Pick a representative node: middle cell, middle node
    cell_mid = cld(nCells, 2)
    node_mid = cld(nNodes, 2)
    idx_mid = (cell_mid-1)*nNodes + node_mid

    tvec = sol.t
    c_num = [u[idx_mid] for u in sol.u]
    c_exact = @. exp(α * tvec) * y0[idx_mid]
    return tvec, c_num, c_exact
end

# --- Always show quick plots when this file runs ---
function _radial_show_plots()
    # Time series at a representative node
    t, cnum, cex = radial_mms_time_series()
    p1 = plot!(t, cnum, label="numeric", xlabel="t", ylabel="c(r*, t)")
    plot!(p1, t, cex, label="exact", linestyle=:dash)
    #savefig(p1, joinpath(@__DIR__, "c_vs_t.png"))
    display(p1)

    # Spatial snapshot: DG residual vs analytical RHS + |error|
    r, Dc_snap, Lc_vals = _radial_mms_snapshot(4, 8; rin=0.01, rout=0.05, v_in=2/60, D=1e-4)
    p2 = plot!(r, Lc_vals, lw=2, label="Analytical RHS F(r)")
    plot!(p2, r, Dc_snap, lw=2, label="DG residual")
    plot!(p2, r, abs.(Dc_snap .- Lc_vals), lw=2, label="|Error|")
    xlabel!(p2, "r [m]"); ylabel!(p2, "Magnitude")
    title!(p2, "Radial DG vs Analytical Residual")
    savefig(p2, joinpath(@__DIR__, "residual_vs_exact.png"))
    display(p2)
end

# Run plots immediately when this file is included/run
_radial_show_plots()

#run_radial_mms_test()
#if get(ENV, "RUN_EOC", "0") == "1"
#    run_radial_mms_eoc_p(p_list = 1:6, nCells = 8)
#    run_radial_mms_eoc_h(nCells_list = [4, 6, 8, 12, 16], polyDeg = 4)
#    res_p = run_radial_mms_eoc_p(p_list = 1:6, nCells = 8)
#    res_h = run_radial_mms_eoc_h(nCells_list = [4,6,8,12,16], polyDeg = 4)
#end

# --- p-refinement ---
#    plot(res_p.p_list .+ 1, res_p.errs, xscale = :log10, yscale = :log10, marker = :o, label = "p-refinement (fixed nCells)", xlabel = "Polynomial degree + 1", ylabel = "Relative L2 error", title = "EOC - p refinement")
#    annotate!(3, 1e-2, text("Slope ≈ $(round(res_p.slope, digits=2))", 10))

    # --- h-refinement ---
#    plot(1.0 ./ ((0.05 - 0.01) ./ res_h.nCells_list), res_h.errs, xscale = :log10, yscale = :log10, marker = :diamond, label = "h-refinement (fixed p=4)", xlabel = "1/h (mesh resolution)", ylabel = "Relative L2 error", title = "EOC - h refinement")
#    annotate!(30, 1e-2, text("Slope ≈ $(round(res_h.slope, digits=2))", 10))


#rel, abs_err = _radial_mms_error(4, 8; rin=0.01, rout=0.05)
#r = range(0.01, 0.05, length=8*(4+1))
#c_fun(r) = (r-0.01)^2 * (0.05-r)^2
#v_of_r(r) = (2/60) * 0.01 / r
#D = 1e-4
#dcdr(r) = 2*((r-0.01)*(0.05-r)^2 - (r-0.01)^2*(0.05-r))
#d2cdr2(r) = 2*((0.05-r)^2 - 4*(r-0.01)*(0.05-r) + (r-0.01)^2)
#Lc(r) = -v_of_r(r)*dcdr(r) + D*(d2cdr2(r) + (1/r)*dcdr(r))
#Lc_vals = Lc.(r)

# --- Spatial comparison plot (numeric residual vs analytical RHS) ---
#r, Dc_snap, Lc_vals = _radial_mms_snapshot(4, 8; rin=0.01, rout=0.05, v_in=2/60, D=1e-4)

#plot(r, Lc_vals, lw=2, label="Analytical RHS F(r)")
#plot!(r, Dc_snap, lw=2, label="DG residual")
#plot!(r, abs.(Dc_snap .- Lc_vals), lw=2, label="|Error|")
#xlabel!("r [m]")
#ylabel!("Magnitude")
#title!("Radial DG vs Analytical Residual")