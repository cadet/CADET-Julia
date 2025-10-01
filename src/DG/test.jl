using CADETJulia
using CADETJulia: DGElements, RadialConvDispOperatorDG
using Plots

function mms_compare_plot(; p=4, nCells=8, ri=0.01, ro=0.05, v_in=1e-4, D=1.0e-8)
    # Manufactured solution and derivatives
    A, B = 2.0, 0.3
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

    # Call DG residual builder (Neumann-0 at outlet, central advective flux)
    RadialConvDispOperatorDG.radialresidualImpl!(
        Dc_num, y, idx,
        1, nNodes, nPoints, nNodes, nCells, deltarho, p,
        invWeights, nodes, polyDerM, invMM, invMrhoM,
        SgMatrix, nothing,                             # diffusion matrix cached; no film MK
        v_in, D,                                       # scalar v_in; constant D
        c_fun(ri),                                     # Dirichlet inlet = exact
        c_star, g_star, Dg, hbuf, mul1,
        RadialConvDispOperatorDG.exact_integration(),
        Inf, 0.0,                                      # Rp=∞, kf=0 → no film exchange
        cp_dummy,
        :neumann0, 0.0,                                # zero-gradient outlet
        rho_i, rho_ip1, ri
    )

    abs_err = maximum(abs.(Dc_num .- F_anal))
    rel_err = abs_err / (maximum(abs.(F_anal)) + eps())
    @info "abs_err = $abs_err, rel_err = $rel_err"

    plt1 = plot(r_all, F_anal, lw=2, label="Analytical RHS",
                xlabel="r [m]", ylabel="RHS",
                title="Radial MMS comparison (p=$p, nCells=$nCells)")
    plot!(plt1, r_all, Dc_num, lw=2, ls=:dash, label="DG residual (numerical)")

    plt2 = plot(r_all, Dc_num .- F_anal, lw=2, label="Pointwise error",
                xlabel="r [m]", ylabel="Error", title="Pointwise error")

    display(plot(plt1, plt2, layout=(2,1)))
    return (r=r_all, F=F_anal, DG=Dc_num, abs_err=abs_err, rel_err=rel_err)
end

# Example:
mms_compare_plot(p=4, nCells=8)