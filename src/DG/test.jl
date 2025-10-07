using CADETJulia
using .CADETJulia.RadialConvDispOperatorDG
using LinearAlgebra
using Test
using Printf

function run_radial_mms_test(; polyDeg=4, nCells=8, rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4)
    @info "Running radial DG test" polyDeg nCells rin rout v_in D

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

    RadialConvDispOperatorDG.radialresidualImpl!(
        Dc, y, 1:nPoints,
        1, nNodes, nNodes, nCells, deltarho, polyDeg,
        invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix,
        v_in, D, cIn, c_star, g_star, Dg, h, mul1,
        RadialConvDispOperatorDG.exact_integration(),
        :neumann0, 0.0,
        rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec
    )

    # Analytic RHS on the same nodes
    L_exact = @. Lc(radii)

    # Robust relative L2 (guard for tiny reference):
    num = norm(Dc - L_exact)
    den = max(norm(L_exact), 1.0)   # scale to 1 if reference is tiny, to avoid bogus huge rel_err
    rel_err = num / den
    abs_err = maximum(abs.(Dc - L_exact))
    @info "comparison" abs_err rel_err tol = 1/60 verdict = (rel_err < 1/60 ? "PASS" : "FAIL")

    @test rel_err < 1/60
    println("Test passed. rel. L2 error = ", rel_err, ", max|err| = ", abs_err)
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

    RadialConvDispOperatorDG.radialresidualImpl!(
        Dc, y, 1:nPoints,
        1, nNodes, nNodes, nCells, deltarho, polyDeg,
        invWeights, nodes, polyDerM, invMM, MM, invMrhoM, SgMatrix,
        v_in, D, cIn, c_star, g_star, Dg, h, mul1,
        RadialConvDispOperatorDG.exact_integration(),
        :neumann0, 0.0,
        rho_i, rho_ip1, faces_v, left_scale_vec, right_scale_vec
    )

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
    println("\nEOC (p-refinement) — MMS for radial DG (nCells=$(nCells))")
    errs = Float64[]; dof_scale = Float64[]  # use (p+1) as resolution metric
    @printf("%6s  %8s  %12s  %12s\n", "p", "nNodes", "rel_L2_err", "abs_max_err")
    @printf("%s\n", repeat('-', 46))
    for p in p_list
        rel, ab = _radial_mms_error(p, nCells; rin=rin, rout=rout, v_in=v_in, D=D)
        push!(errs, rel); push!(dof_scale, p+1)
        @printf("%6d  %8d  %12.4e  %12.4e\n", p, (p+1)*nCells, rel, ab)
    end
    slope = _loglog_slope(dof_scale, errs)
    @printf("Estimated slope on log(err) vs log(p+1): %.3f\n", slope)
    return (; p_list = collect(p_list), errs, slope)
end

# --- h-refinement EOC: vary nCells at fixed polynomial degree ---
function run_radial_mms_eoc_h(; nCells_list = [4, 6, 8, 12, 16, 24], polyDeg=4, rin=0.01, rout=0.05, v_in=2/60, D=1.0e-4)
    println("\nEOC (h-refinement) — MMS for radial DG (polyDeg=$(polyDeg))")
    errs = Float64[]; hvals = Float64[]  # characteristic cell size Δr
    @printf("%8s  %8s  %12s  %12s\n", "nCells", "(p+1)N", "rel_L2_err", "abs_max_err")
    @printf("%s\n", repeat('-', 48))
    for nc in nCells_list
        rel, ab = _radial_mms_error(polyDeg, nc; rin=rin, rout=rout, v_in=v_in, D=D)
        push!(errs, rel); push!(hvals, (rout-rin)/nc)
        @printf("%8d  %8d  %12.4e  %12.4e\n", nc, (polyDeg+1)*nc, rel, ab)
    end
    slope = _loglog_slope(1.0 ./ hvals, errs)  # slope wrt resolution ~ 1/h
    @printf("Estimated slope on log(err) vs log(1/h): %.3f\n", slope)
    return (; nCells_list = collect(nCells_list), errs, slope)
end

# Run on include
run_radial_mms_test()
# Optional: quick EOC sweeps (comment out for CI if desired)
run_radial_mms_eoc_p(p_list = 1:6, nCells = 8)
run_radial_mms_eoc_h(nCells_list = [4, 6, 8, 12, 16], polyDeg = 4)