#!/usr/bin/env julia
# Detailed Analysis: Compare radial vs axial flow for multiple R ratios
# Shows detailed concentration profiles and statistics

using Printf
using LinearAlgebra
using OrdinaryDiffEqBDF
using SciMLBase

push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))
using CADETJulia: DGElements

include(joinpath(@__DIR__, "..", "..", "src", "DG", "RadialConvDispOperatorDG.jl"))
using .RadialConvDispOperatorDG

include(joinpath(@__DIR__, "..", "..", "src", "DG", "ConvDispOperatorDG.jl"))
using .ConvDispOperatorDG

println("=" ^ 80)
println("DETAILED RADIAL VS AXIAL COMPARISON - MULTIPLE R RATIOS")
println("=" ^ 80)
println()

# ==================== PARAMETERS ====================
L = 0.1            # Column length [m]
r_in_base = 0.1    # Base inner radius [m]
v_in = 1.0e-2      # Inlet velocity for radial [m/s]
D = 1.0e-5         # Dispersion coefficient [m²/s]

# Pulse parameters
t_inject_start = 0.0
t_inject_end = 5.0
c_inject = 1.0

function inlet_pulse(t)
    if t_inject_start <= t <= t_inject_end
        return c_inject
    else
        return 0.0
    end
end

# Discretization
polyDeg = 4
nCells = 8
nNodes = polyDeg + 1
nPoints = nCells * nNodes

abstol = 1e-12
reltol = 1e-10

# Test these R ratios
radius_ratios = [5.0, 4.0, 3.0, 2.5, 2.25, 2.0, 1.9, 1.8, 1.75, 1.5, 1.1]

println("Testing R values: ", radius_ratios)
println()
println("Parameters:")
@printf("  L = %.2f m, v_in = %.4f m/s, D = %.2e m²/s\n", L, v_in, D)
@printf("  p=%d, nCells=%d, DOFs=%d\n", polyDeg, nCells, nPoints)
println()

# Storage
all_results = []

for R in radius_ratios
    println("=" ^ 80)
    @printf("ANALYZING R = %.2f\n", R)
    println("=" ^ 80)
    println()

    # Calculate radii
    r_in = r_in_base / (R - 1)
    r_out = R * r_in
    L_actual = r_out - r_in

    @printf("Geometry: r_in=%.4f m, r_out=%.4f m, L=%.4f m\n", r_in, r_out, L_actual)

    # Residence time
    tau_rad = (L_actual / v_in) * log(R)
    v_axial = L_actual / tau_rad

    @printf("Residence time: τ_rad = %.4f s\n", tau_rad)
    @printf("Matched axial velocity: v_axial = %.6f m/s\n", v_axial)
    @printf("Radial inlet velocity: v_in = %.6f m/s\n", v_in)
    println()

    t_final = max(100.0, 4 * tau_rad)

    # ==================== SOLVE AXIAL FLOW ====================
    deltaZ = L_actual / nCells

    nodes_ax, invWeights_ax = DGElements.lglnodes(polyDeg)
    polyDerM_ax = DGElements.derivativeMatrix(polyDeg, nodes_ax)
    invMM_ax = DGElements.invMMatrix(nodes_ax, polyDeg, 0, 0)

    z_nodes = zeros(nPoints)
    for cell in 1:nCells
        for node in 1:nNodes
            idx = (cell - 1) * nNodes + node
            z_nodes[idx] = (cell - 1) * deltaZ + (deltaZ / 2) * (1 + nodes_ax[node])
        end
    end

    idx_ax = [1]
    c_star_ax = zeros(nCells + 1)
    h_star_ax = zeros(nCells + 1)
    Dc_ax = zeros(nPoints)
    _h_ax = zeros(nPoints)
    mul1_ax = zeros(nNodes)

    function ode_rhs_axial!(dc_dt, c, p, t)
        cIn = inlet_pulse(t)
        ConvDispOperatorDG.residualImpl!(
            dc_dt, c, idx_ax, 1, nNodes, nPoints, nNodes, nCells,
            deltaZ, polyDeg, invWeights_ax, polyDerM_ax, invMM_ax,
            v_axial, D, cIn, c_star_ax, h_star_ax, Dc_ax, _h_ax, mul1_ax,
            ConvDispOperatorDG.exact_integration()
        )
        return nothing
    end

    print("Solving axial flow... ")
    c_initial_ax = zeros(nPoints)
    tspan = (0.0, t_final)
    prob_ax = ODEProblem(ode_rhs_axial!, c_initial_ax, tspan)
    sol_ax = solve(prob_ax, QNDF(autodiff=false), abstol=abstol, reltol=reltol,
                   save_everystep=false, save_start=false)
    println("✓")

    c_axial = sol_ax.u[end]

    # ==================== SOLVE RADIAL FLOW ====================
    deltarho = L_actual / nCells

    nodes_rad, invWeights_rad = DGElements.lglnodes(polyDeg)
    polyDerM_rad = DGElements.derivativeMatrix(polyDeg, nodes_rad)
    invMM_rad = DGElements.invMMatrix(nodes_rad, polyDeg, 0, 0)
    MM00_rad = DGElements.MMatrix(nodes_rad, polyDeg, 0, 0)
    MM01_rad = DGElements.MMatrix(nodes_rad, polyDeg, 0, 1)

    rho_i = range(r_in, r_out, length=nCells+1) |> collect
    rho_nodes = zeros(nPoints)
    for cell in 1:nCells
        for node in 1:nNodes
            idx = (cell - 1) * nNodes + node
            rho_nodes[idx] = rho_i[cell] + (deltarho / 2) * (1 + nodes_rad[node])
        end
    end

    rMM, invrMM = DGElements.weightedMMatrix(nodes_rad, polyDeg, rho_i, deltarho)
    S_g = DGElements.dispMMatrix(nodes_rad, polyDeg, rho_i, deltarho, D, polyDerM_rad, rMM)

    idx_rad = [1]
    c_star_rad = zeros(nCells + 1)
    h_star_rad = zeros(nCells + 1)
    Dg_rad = zeros(nPoints)
    g_phys_rad = zeros(nPoints)
    mul1_rad = zeros(nNodes)

    function ode_rhs_radial!(dc_dt, c, p, t)
        cIn = inlet_pulse(t)
        RadialConvDispOperatorDG.radialresidualImpl!(
            dc_dt, c, idx_rad, 1, nNodes, nNodes, nCells,
            deltarho, polyDeg, polyDerM_rad, invMM_rad, MM01_rad, MM00_rad,
            rMM, invrMM, S_g, nodes_rad, invWeights_rad,
            v_in, D, rho_i, cIn, c_star_rad, h_star_rad, Dg_rad, g_phys_rad, mul1_rad
        )
        return nothing
    end

    print("Solving radial flow... ")
    c_initial_rad = zeros(nPoints)
    prob_rad = ODEProblem(ode_rhs_radial!, c_initial_rad, tspan)
    sol_rad = solve(prob_rad, QNDF(autodiff=false), abstol=abstol, reltol=reltol,
                    save_everystep=false, save_start=false)
    println("✓")
    println()

    c_radial = sol_rad.u[end]

    # ==================== ANALYSIS ====================
    diff = c_radial .- c_axial
    L2_error = sqrt(sum(diff.^2) / length(diff))
    Linf_error = maximum(abs.(diff))
    outlet_error = abs(c_radial[end] - c_axial[end])

    c_max = maximum(abs.([c_axial; c_radial]))
    L2_rel = c_max > 1e-15 ? 100 * L2_error / c_max : 0.0
    outlet_rel = c_max > 1e-15 ? 100 * outlet_error / c_max : 0.0

    # Velocity analysis
    v_inlet_rad = v_in
    v_outlet_rad = v_in * r_in / r_out
    v_variation = 100 * (1 - v_outlet_rad / v_inlet_rad)

    println("RESULTS:")
    @printf("  L2 error:     %.6e (%.2f%%)\n", L2_error, L2_rel)
    @printf("  L∞ error:     %.6e\n", Linf_error)
    @printf("  Outlet error: %.6e (%.2f%%)\n", outlet_error, outlet_rel)
    println()
    @printf("  Velocity variation (radial): %.2f%%\n", v_variation)
    @printf("  v_inlet = %.6f m/s, v_outlet = %.6f m/s\n", v_inlet_rad, v_outlet_rad)
    println()

    # Store results
    push!(all_results, (
        R = R,
        r_in = r_in,
        r_out = r_out,
        L = L_actual,
        tau_rad = tau_rad,
        v_axial = v_axial,
        v_inlet_rad = v_inlet_rad,
        v_outlet_rad = v_outlet_rad,
        v_variation = v_variation,
        L2_error = L2_error,
        L2_rel = L2_rel,
        Linf_error = Linf_error,
        outlet_error = outlet_error,
        outlet_rel = outlet_rel,
        c_axial_max = maximum(abs.(c_axial)),
        c_radial_max = maximum(abs.(c_radial))
    ))
end

# ==================== SAVE RESULTS ====================
println("=" ^ 80)
println("SAVING RESULTS")
println("=" ^ 80)
println()

output_file = "test12_radial_vs_axial_detailed.csv"
open(output_file, "w") do io
    println(io, "# Detailed Radial vs Axial Comparison")
    println(io, "# L=$L m, v_in=$v_in m/s, D=$D m²/s, p=$polyDeg, nCells=$nCells")
    println(io, "R,r_in,r_out,L,tau_rad,v_axial,v_inlet_rad,v_outlet_rad,v_variation,L2_error,L2_rel,Linf_error,outlet_error,outlet_rel,c_axial_max,c_radial_max")

    for res in all_results
        @printf(io, "%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.4f,%.6e,%.4f,%.6e,%.6e,%.4f,%.6e,%.6e\n",
                res.R, res.r_in, res.r_out, res.L, res.tau_rad, res.v_axial,
                res.v_inlet_rad, res.v_outlet_rad, res.v_variation,
                res.L2_error, res.L2_rel, res.Linf_error,
                res.outlet_error, res.outlet_rel, res.c_axial_max, res.c_radial_max)
    end
end

println("Saved to: $output_file")
println()

# ==================== SUMMARY TABLE ====================
println("=" ^ 80)
println("SUMMARY")
println("=" ^ 80)
println()

@printf("%-8s %-12s %-15s %-15s %-15s\n", "R", "τ_rad [s]", "v_var [%]", "L2 Error [%]", "Outlet [%]")
println("-" ^ 75)
for res in all_results
    @printf("%-8.2f %-12.4f %-15.2f %-15.4f %-15.4f\n",
            res.R, res.tau_rad, res.v_variation, res.L2_rel, res.outlet_rel)
end

println()

# Highlight R=2.0
idx_R2 = findfirst(res -> abs(res.R - 2.0) < 0.01, all_results)
if idx_R2 !== nothing
    println("✓ At R = 2.0 (literature prediction):")
    @printf("  L2 error: %.4f%%, Outlet error: %.4f%%\n",
            all_results[idx_R2].L2_rel, all_results[idx_R2].outlet_rel)
    if all_results[idx_R2].L2_rel < 5.0
        println("  EXCELLENT agreement!")
    end
end

println()
println("=" ^ 80)
println("ANALYSIS COMPLETE - Run the Python script to visualize!")
println("=" ^ 80)
