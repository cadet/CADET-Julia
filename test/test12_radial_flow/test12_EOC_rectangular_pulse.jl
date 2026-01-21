#!/usr/bin/env julia
# EOC (Experimental Order of Convergence) Test for Radial Flow
# Tests spatial convergence using a rectangular pulse injection

using Printf
using LinearAlgebra
using OrdinaryDiffEqBDF
using SciMLBase

# Include DGElements module directly
include(joinpath(@__DIR__, "..", "..", "src", "DG", "DGElements.jl"))
using .DGElements

# Load the radial convection-dispersion operator
include(joinpath(@__DIR__, "..", "..", "src", "DG", "RadialConvDispOperatorDG.jl"))
using .RadialConvDispOperatorDG

# ==================== PHYSICAL PARAMETERS ====================
rin = 0.1
rout = 1.1
v = 1.0e-2       # Interstitial velocity [m/s]
D = 1.0e-4       # Dispersion coefficient [m²/s]

# ==================== RECTANGULAR PULSE PARAMETERS ====================
t_start = 10.0       # Start time of rectangular pulse [s]
t_end = 50.0         # End time of rectangular pulse [s]
c_max = 1.0          # Maximum concentration

function inlet_rectangular(t)
    if t >= t_start && t <= t_end
        return c_max
    else
        return 0.0
    end
end

# ==================== SIMULATION PARAMETERS ====================
t_final = 150.0
abstol = 1e-12
reltol = 1e-12

# ==================== EOC STUDY CONFIGURATION ====================
polyDeg_list = [1, 2, 3, 4, 5]
nCells_list = [1, 2, 4, 8, 16, 32, 64, 128, 256]

# ==================== PRINT HEADER ====================
println("EOC Test - Rectangular Pulse (Radial Flow)")
println("Parameters: rin=$rin, rout=$rout, v=$v, D=$D, t_final=$t_final")
println("Pulse: t_start=$t_start, t_end=$t_end")
println()
@printf("%-8s %-8s %-14s %-14s %-10s %-10s\n", "polyDeg", "nCells", "L2_err", "Linf_err", "EOC_L2", "EOC_Linf")
println("-" ^ 70)

# ==================== LOOP OVER POLYNOMIAL DEGREES ====================
for polyDeg in polyDeg_list
    n_levels = length(nCells_list)
    solutions = []
    rho_nodes_list = []

    # Solve at each refinement level
    for (level, nCells) in enumerate(nCells_list)
        nNodes = polyDeg + 1
        nPoints = nCells * nNodes
        deltarho = (rout - rin) / nCells

        # Setup DG discretization
        nodes, invWeights = DGElements.cglnodes(polyDeg)
        polyDerM = DGElements.derivativeMatrix(polyDeg, nodes)
        invMM = DGElements.invMMatrix(nodes, polyDeg, 0, 0)
        MM00 = DGElements.MMatrix(nodes, polyDeg, 0, 0)
        MM01 = DGElements.MMatrix(nodes, polyDeg, 0, 1)

        # Cell interfaces and node positions
        rho_i = range(rin, rout, length=nCells+1) |> collect
        rho_nodes = zeros(nPoints)
        for Cell in 1:nCells
            for node in 1:nNodes
                idx = (Cell-1) * nNodes + node
                rho_nodes[idx] = rho_i[Cell] + (deltarho/2) * (1 + nodes[node])
            end
        end

        # Weighted mass matrices
        rMM = DGElements.weightedMMatrix(nodes, polyDeg, nCells, rho_i, deltarho)
        invrMM = DGElements.invweightedMMatrix(nodes, polyDeg, nCells, rho_i, deltarho)
        S_g = DGElements.dispMMatrix(nodes, polyDeg, rho_i, deltarho, D, polyDerM, rMM)

        # Pre-allocate workspace
        idx_comp = [1]
        c_star = zeros(nCells + 1)
        g_star = zeros(nCells + 1)
        Dg = zeros(nPoints)
        g = zeros(nPoints)
        mul1 = zeros(nNodes)
        mul2 = zeros(nNodes)

        # ODE right-hand side
        function ode_rhs!(dc_dt, c, p, t)
            cIn = inlet_rectangular(t)
            RadialConvDispOperatorDG.radialresidualImpl!(
                dc_dt, c, idx_comp,
                1, nNodes, nNodes, nCells,
                deltarho, polyDeg,
                polyDerM, invMM, MM01, MM00,
                rMM, invrMM, S_g,
                nodes, invWeights,
                v, D, rho_i, cIn,
                c_star, g_star, Dg, g, mul1, mul2
            )
            return nothing
        end

        # Solve ODE
        c_initial = zeros(nPoints)
        tspan = (0.0, t_final)
        prob = ODEProblem(ode_rhs!, c_initial, tspan)
        sol = solve(prob, QNDF(autodiff=false), abstol=abstol, reltol=reltol,
                    save_everystep=false, save_start=false)

        push!(solutions, sol.u[end])
        push!(rho_nodes_list, rho_nodes)
    end

    # Compute EOC using finest mesh as reference
    c_ref = solutions[end]
    rho_ref = rho_nodes_list[end]
    nCells_ref = nCells_list[end]
    nNodes_ref = polyDeg + 1
    deltarho_ref = (rout - rin) / nCells_ref

    # High-order Lagrange interpolation
    function lagrange_interp(rho_target, c_ref, rho_ref_nodes, cell_start, nNodes)
        result = 0.0
        for j in 1:nNodes
            Lj = 1.0
            for k in 1:nNodes
                if k != j
                    Lj *= (rho_target - rho_ref_nodes[cell_start + k - 1]) /
                          (rho_ref_nodes[cell_start + j - 1] - rho_ref_nodes[cell_start + k - 1])
                end
            end
            result += c_ref[cell_start + j - 1] * Lj
        end
        return result
    end

    errors_L2 = Float64[]
    errors_Linf = Float64[]
    h_list = Float64[]

    for level in 1:(n_levels-1)
        nCells = nCells_list[level]
        h = (rout - rin) / nCells
        c_coarse = solutions[level]
        rho_coarse = rho_nodes_list[level]

        c_ref_interp = zeros(length(rho_coarse))
        for i in 1:length(rho_coarse)
            rho_target = rho_coarse[i]
            cell_ref = min(nCells_ref, max(1, Int(floor((rho_target - rin) / deltarho_ref)) + 1))
            if rho_target <= rin
                cell_ref = 1
            elseif rho_target >= rout
                cell_ref = nCells_ref
            end
            cell_start = (cell_ref - 1) * nNodes_ref + 1
            c_ref_interp[i] = lagrange_interp(rho_target, c_ref, rho_ref, cell_start, nNodes_ref)
        end

        diff = c_coarse .- c_ref_interp
        L2_error = sqrt(sum(diff.^2) / length(diff))
        Linf_error = maximum(abs.(diff))

        push!(errors_L2, L2_error)
        push!(errors_Linf, Linf_error)
        push!(h_list, h)
    end

    # Print results for this polynomial degree
    for i in 1:(n_levels-1)
        nCells = nCells_list[i]
        if i == 1
            @printf("%-8d %-8d %.6e  %.6e  %-10s %-10s\n",
                    polyDeg, nCells, errors_L2[i], errors_Linf[i], "--", "--")
        else
            eoc_L2 = log(errors_L2[i-1] / errors_L2[i]) / log(h_list[i-1] / h_list[i])
            eoc_Linf = log(errors_Linf[i-1] / errors_Linf[i]) / log(h_list[i-1] / h_list[i])
            @printf("%-8d %-8d %.6e  %.6e  %-10.2f %-10.2f\n",
                    polyDeg, nCells, errors_L2[i], errors_Linf[i], eoc_L2, eoc_Linf)
        end
    end
end

println("-" ^ 70)
println("Note: Rectangular pulse has discontinuities - EOC may be reduced for high p")
