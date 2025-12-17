#!/usr/bin/env julia
# EOC (Experimental Order of Convergence) Test for Radial Flow
# Tests spatial convergence using a smooth Gaussian pulse injection
# Compares solutions at different mesh refinements

using Printf
using LinearAlgebra
using OrdinaryDiffEqBDF
using SciMLBase

# Add the CADETJulia module path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))
using CADETJulia: DGElements

# Load the radial convection-dispersion operator
include(joinpath(@__DIR__, "..", "..", "src", "DG", "RadialConvDispOperatorDG.jl"))
using .RadialConvDispOperatorDG

println("=" ^ 80)
println("EOC TEST - GAUSSIAN PULSE INJECTION")
println("Spatial Convergence Study for Radial Flow")
println("=" ^ 80)
println()

# ==================== PHYSICAL PARAMETERS ====================
rin = 0.1
rout = 1.0
v = 0.01       # Interstitial velocity [m/s]
D = 1.0e-5     # Dispersion coefficient [m²/s]

println("Physical Parameters:")
println("  Inner radius: $rin m")
println("  Outer radius: $rout m")
println("  Velocity: $v m/s")
println("  Dispersion: $D m²/s")
println()

# ==================== GAUSSIAN PULSE PARAMETERS ====================
t_center = 5.0      # Center time of Gaussian pulse [s]
t_width = 2.0       # Width (standard deviation) of Gaussian [s]
c_max = 1.0         # Maximum concentration

function inlet_gaussian(t)
    """Smooth Gaussian pulse inlet boundary condition."""
    return c_max * exp(-((t - t_center) / t_width)^2)
end

println("Inlet Boundary Condition:")
println("  Gaussian pulse: c(t) = $c_max * exp(-((t - $t_center) / $t_width)²)")
println()

# ==================== SIMULATION PARAMETERS ====================
t_final = 100.0
abstol = 1e-10
reltol = 1e-8

println("Simulation Parameters:")
println("  Final time: $t_final s")
println("  ODE tolerances: abstol=$abstol, reltol=$reltol")
println()

# ==================== EOC STUDY CONFIGURATION ====================
# Loop over polynomial degrees
polyDeg_list = [1, 2, 3, 4, 5]

# Mesh refinement levels: nCells = 4, 8, 16, 32
nCells_list = [1, 2, 4, 8, 16, 32, 64, 128]

println("EOC Configuration:")
println("  Polynomial degrees: ", polyDeg_list)
println("  Mesh levels: ", nCells_list)
println("  Convergence type: h-refinement")
println()

# ==================== LOOP OVER POLYNOMIAL DEGREES ====================
for polyDeg in polyDeg_list
    println("=" ^ 80)
    println("POLYNOMIAL DEGREE: $polyDeg")
    println("=" ^ 80)
    println()

    n_levels = length(nCells_list)

    # Storage for solutions at each level
    solutions = []
    rho_nodes_list = []
    nPoints_list = []
    compute_times = []  # Store computation time for each level

    # ==================== SOLVE AT EACH REFINEMENT LEVEL ====================
    println("=" ^ 80)
    println("COMPUTING SOLUTIONS AT EACH REFINEMENT LEVEL")
    println("=" ^ 80)
    println()

    for (level, nCells) in enumerate(nCells_list)
        println("-" ^ 80)
        println("Level $level: nCells = $nCells")
        println("-" ^ 80)

        nNodes = polyDeg + 1
        nPoints = nCells * nNodes
        deltarho = (rout - rin) / nCells

        println("  Total DOFs: $nPoints")
        println("  Cell width: $deltarho m")

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
        rMM, invrMM = DGElements.weightedMMatrix(nodes, polyDeg, rho_i, deltarho)

        # Dispersion stiffness matrix
        S_g = DGElements.dispMMatrix(nodes, polyDeg, rho_i, deltarho, D, polyDerM, rMM)

        # Pre-allocate workspace
        idx_comp = [1]
        c_star = zeros(nCells + 1)
        h_star = zeros(nCells + 1)
        Dg = zeros(nPoints)
        g_phys = zeros(nPoints)
        mul1 = zeros(nNodes)

        # ODE right-hand side
        function ode_rhs!(dc_dt, c, p, t)
            cIn = inlet_gaussian(t)

            RadialConvDispOperatorDG.radialresidualImpl!(
                dc_dt, c, idx_comp,
                1, nNodes, nNodes, nCells,
                deltarho, polyDeg,
                polyDerM, invMM, MM01, MM00,
                rMM, invrMM, S_g,
                nodes, invWeights,
                v, D, rho_i, cIn,
                c_star, h_star, Dg, g_phys, mul1
            )

            return nothing
        end

        # Initial condition: zero everywhere
        c_initial = zeros(nPoints)

        # Solve ODE
        println("  Solving...")
        tspan = (0.0, t_final)
        prob = ODEProblem(ode_rhs!, c_initial, tspan)

        # Time the solve
        time_start = time()
        sol = solve(prob, QNDF(autodiff=false), abstol=abstol, reltol=reltol,
                    save_everystep=false, save_start=false)
        time_elapsed = time() - time_start

        println("  Solution complete!")
        println("    Compute time: $(round(time_elapsed, digits=3)) s")
        println("    Function evaluations: $(sol.destats.nf)")
        println("    Accepted steps: $(sol.destats.naccept)")
        println()

        # Store solution
        push!(solutions, sol.u[end])
        push!(rho_nodes_list, rho_nodes)
        push!(nPoints_list, nPoints)
        push!(compute_times, time_elapsed)
    end

    # ==================== COMPUTE EOC ====================
    println("=" ^ 80)
    println("COMPUTING EXPERIMENTAL ORDER OF CONVERGENCE")
    println("=" ^ 80)
    println()

    println("Using finest mesh (level $n_levels) as reference solution")
    println()

    # Use finest mesh as reference
    c_ref = solutions[end]
    rho_ref = rho_nodes_list[end]

    # Compute errors relative to finest solution
    errors_L2 = Float64[]
    errors_Linf = Float64[]
    h_list = Float64[]

    for level in 1:(n_levels-1)
        nCells = nCells_list[level]
        h = (rout - rin) / nCells  # Characteristic mesh size

        c_coarse = solutions[level]
        rho_coarse = rho_nodes_list[level]

        # Interpolate reference solution to coarse grid points
        c_ref_interp = zeros(length(rho_coarse))
        for i in 1:length(rho_coarse)
            # Find reference solution at this position via linear interpolation
            rho_target = rho_coarse[i]

            # Find bracketing indices in reference grid
            idx_right = findfirst(r -> r >= rho_target, rho_ref)
            if idx_right === nothing
                c_ref_interp[i] = c_ref[end]
            elseif idx_right == 1
                c_ref_interp[i] = c_ref[1]
            else
                idx_left = idx_right - 1
                # Linear interpolation
                alpha = (rho_target - rho_ref[idx_left]) / (rho_ref[idx_right] - rho_ref[idx_left])
                c_ref_interp[i] = (1 - alpha) * c_ref[idx_left] + alpha * c_ref[idx_right]
            end
        end

        # Compute errors
        diff = c_coarse .- c_ref_interp
        L2_error = sqrt(sum(diff.^2) / length(diff))
        Linf_error = maximum(abs.(diff))

        push!(errors_L2, L2_error)
        push!(errors_Linf, Linf_error)
        push!(h_list, h)

        println("Level $level (nCells = $nCells, h = $(round(h, digits=6))):")
        @printf("  L2 error:   %.6e\n", L2_error)
        @printf("  L∞ error:   %.6e\n", Linf_error)
        println()
    end

    # ==================== COMPUTE CONVERGENCE RATES ====================
    println("-" ^ 80)
    println("CONVERGENCE RATES")
    println("-" ^ 80)
    println()

    @printf("%-8s %-12s %-14s %-14s %-10s %-10s\n",
            "Level", "nCells", "L2 error", "L∞ error", "EOC(L2)", "EOC(L∞)")
    println("-" ^ 80)

    for i in 1:(n_levels-1)
        nCells = nCells_list[i]

        if i == 1
            @printf("%-8d %-12d %.6e  %.6e  %-10s %-10s\n",
                    i, nCells, errors_L2[i], errors_Linf[i], "--", "--")
        else
            # Compute convergence rate: log(e1/e2) / log(h1/h2)
            eoc_L2 = log(errors_L2[i-1] / errors_L2[i]) / log(h_list[i-1] / h_list[i])
            eoc_Linf = log(errors_Linf[i-1] / errors_Linf[i]) / log(h_list[i-1] / h_list[i])

            @printf("%-8d %-12d %.6e  %.6e  %-10.2f %-10.2f\n",
                    i, nCells, errors_L2[i], errors_Linf[i], eoc_L2, eoc_Linf)
        end
    end

    println()

    # ==================== SAVE RESULTS ====================
    output_file = "test12_EOC_gaussian_pulse_p$(polyDeg).csv"
    open(output_file, "w") do io
        println(io, "# EOC Test - Gaussian Pulse Injection in Radial Flow")
        println(io, "# rin=$rin m, rout=$rout m, v=$v m/s, D=$D m²/s")
        println(io, "# polyDeg=$polyDeg (fixed), t_final=$t_final s")
        println(io, "# Gaussian pulse: center=$t_center s, width=$t_width s")
        println(io, "# level,nCells,h,L2_error,Linf_error,EOC_L2,EOC_Linf,compute_time")

        for i in 1:(n_levels-1)
            nCells = nCells_list[i]
            h = h_list[i]
            comp_time = compute_times[i]

            if i == 1
                @printf(io, "%d,%d,%.6e,%.6e,%.6e,,,%.6f\n",
                        i, nCells, h, errors_L2[i], errors_Linf[i], comp_time)
            else
                eoc_L2 = log(errors_L2[i-1] / errors_L2[i]) / log(h_list[i-1] / h_list[i])
                eoc_Linf = log(errors_Linf[i-1] / errors_Linf[i]) / log(h_list[i-1] / h_list[i])

                @printf(io, "%d,%d,%.6e,%.6e,%.6e,%.2f,%.2f,%.6f\n",
                        i, nCells, h, errors_L2[i], errors_Linf[i], eoc_L2, eoc_Linf, comp_time)
            end
        end
    end

    println("Results saved to: $output_file")
    println()

    # ==================== SUMMARY ====================
    println("=" ^ 80)
    println("SUMMARY")
    println("=" ^ 80)
    println()

    # Compute average EOC (excluding first level)
    if length(errors_L2) >= 2
        avg_eoc_L2 = 0.0
        avg_eoc_Linf = 0.0
        for i in 2:(n_levels-1)
            eoc_L2 = log(errors_L2[i-1] / errors_L2[i]) / log(h_list[i-1] / h_list[i])
            eoc_Linf = log(errors_Linf[i-1] / errors_Linf[i]) / log(h_list[i-1] / h_list[i])
            avg_eoc_L2 += eoc_L2
            avg_eoc_Linf += eoc_Linf
        end
        avg_eoc_L2 /= (n_levels - 2)
        avg_eoc_Linf /= (n_levels - 2)

        @printf("Average EOC (L2):   %.2f\n", avg_eoc_L2)
        @printf("Average EOC (L∞):   %.2f\n", avg_eoc_Linf)
        println()
    end

    println("Expected convergence rate:")
    println("  For DG with polynomial degree p=$polyDeg:")
    println("  - Optimal rate: p+1 = $(polyDeg+1)")
    println("  - Minimum acceptable: p = $polyDeg")
    println()

    # Check if we achieved expected convergence
    if length(errors_L2) >= 2
        final_eoc_L2 = log(errors_L2[end-1] / errors_L2[end]) / log(h_list[end-1] / h_list[end])
        if final_eoc_L2 >= polyDeg
            println("✓ PASS: Achieved expected convergence rate (EOC ≥ $polyDeg)")
        else
            println("⚠ WARNING: Convergence rate lower than expected")
        end
    end

    println()
end  # End of polyDeg loop

println()
println("=" ^ 80)
println("EOC TEST COMPLETE - ALL POLYNOMIAL DEGREES")
println("=" ^ 80)
