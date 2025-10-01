using Test, CADETJulia, CSV, DataFrames
using CADETJulia: DGElements, RadialConvDispOperatorDG

function run_radial_mms(; polyDeg_list = 1:6, nCells_list  = (2, 4, 8), ri = 0.01, ro = 0.03, v_in = 1.0, D = 1.0e-4)

    # Manufactured c(r) and derivatives
    A = 2.0
    B = 0.3
    c_fun(r)  = A * (r - ri)^2 * (ro - r)^2 + B
    # First derivative and second derivative (derived explicitly)
    # Let s = r - ri, t = ro - r. For c = A*s^2*t^2:
    # dc/dr = 2A*( s*t^2 - s^2*t )
    # d2c/dr2 = 2A*( t^2 - 4*s*t + s^2 )
    dcdr(r) = 2A * ((r - ri)*(ro - r)^2 - (r - ri)^2*(ro - r))
    d2cdr2(r) = 2A * ((ro - r)^2 - 4*(r - ri)*(ro - r) + (r - ri)^2)

    # Analytical RHS F(c) at radius r:
    v_of_r(r) = v_in * (ri / r)  # matches code's internal rule for scalar v
    F_fun(r) = -(v_of_r(r)/r) * dcdr(r) + D*( d2cdr2(r) + (1/r)*dcdr(r) )

    # Helper: map LGL nodes in each cell to physical radii
    function cell_nodes(polyDeg, nCells, ri, ro)
        nodes, invw = DGElements.lglnodes(polyDeg)  # ξ ∈ [-1,1]
        deltarho = (ro - ri) / nCells
        rho_i  = [ri + (k-1)*deltarho for k in 1:nCells]
        # Return a flat vector of radii at all nodal points, cell by cell
        r_all = Float64[]
        for k in 1:nCells
            rho_left = rho_i[k]
            for ξ in nodes
                r = rho_left + (ξ + 1.0) * (deltarho / 2.0)
                push!(r_all, r)
            end
        end
        return nodes, r_all, deltarho, rho_i
    end

    # Set up fixed “DG operator” buffers like radialresidualImpl! expects
    errs = Dict{Tuple{Int,Int}, Float64}()

    for nCells in nCells_list, p in polyDeg_list
        nNodes  = p + 1
        nPoints = nNodes * nCells

        nodes, r_all, deltarho, rho_i = cell_nodes(p, nCells, ri, ro)
        invWeights = DGElements.lglnodes(p)[2]
        invMM      = DGElements.invMMatrix(nodes, p)
        polyDerM   = DGElements.derivativeMatrix(p, nodes)
        invMrhoM   = [DGElements.invMrhoMatrix(nodes, p, deltarho, rho_i[c]) for c in 1:nCells]
        SgMatrix   = [DGElements.weighted_stiff_Matrix(nodes, p, rho_i[c], deltarho, _ -> D) for c in 1:nCells]

        # work buffers
        c_star = zeros(nCells+1)
        g_star = zeros(nCells+1)
        Dc_num = zeros(nPoints)
        Dg     = zeros(nPoints)
        hbuf   = zeros(nPoints)
        mul1   = zeros(nNodes)

        # State vector y for one “component block” (what radialresidualImpl! slices via idx)
        y = [c_fun(r) for r in r_all]

        # cp, film, etc. are all neutralized by Rp=Inf, kf=0
        cp_dummy = copy(y)

        # Build analytical RHS at nodes
        F_anal = [F_fun(r) for r in r_all]

        # Call DG residual builder
        idx = 1:nPoints
        RadialConvDispOperatorDG.radialresidualImpl!(
            Dc_num,              # OUT
            y, idx,              # state + slice idx
            1,                   # strideNode
            nNodes,              # strideCell
            nPoints, nNodes, nCells, deltarho, p,
            invWeights, nodes, polyDerM, invMM, invMrhoM,
            SgMatrix, nothing,   # S_g cached, MKMatrix = nothing
            v_in,                # scalar v_in ⇒ code maps to v(r)=v_in*ri/r
            D,                   # constant diffusion
            c_fun(ri),           # inlet Dirichlet = exact c at ri
            c_star, g_star, Dg, hbuf, mul1,
            RadialConvDispOperatorDG.exact_integration(),
            Inf,                 # Rp = ∞  ⇒ film term disabled (prefactor 0)
            0.0,                 # kf = 0  ⇒ no film exchange
            cp_dummy,            # cp placeholder (unused when film off)
            :neumann0, 0.0,      # outlet BC: zero gradient (matches c'(ro)=0)
            rho_i, [ri + k*deltarho for k in 1:nCells], ri
        )

        # Compute error norms
        abs_err = maximum(abs.(Dc_num .- F_anal))
        rel_err = abs_err / (maximum(abs.(F_anal)) + eps())

        @info "MMS radial test" nCells p abs_err rel_err
        errs[(nCells,p)] = rel_err

        # Basic assertions (tolerances loosen for low p / coarse grids)
        tol = 5e-2 / max(1, p-1)  # heuristic: tighter with higher p
        @test rel_err < tol
    end

    return errs
end

@testset "Radial DGSEM MMS (convection-diffusion, const D, central flux)" begin
    errs = run_radial_mms()
    @info "Relative errors by (nCells, polyDeg)" errs
end