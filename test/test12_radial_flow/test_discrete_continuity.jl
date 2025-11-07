#!/usr/bin/env julia
# Test: Discrete Geometric Conservation Law
#
# For radial flow with continuity ∂(ρv)/∂ρ = 0, the discrete operator
# should satisfy: D^T * M * (v .* 1) ≈ 0
#
# This is a necessary condition for free-stream preservation.

using CADETJulia
using LinearAlgebra
using Printf

println("=" ^ 80)
println("DISCRETE CONTINUITY TEST - Geometric Conservation Law")
println("=" ^ 80)
println()

# Test parameters
polyDeg = 4
nCells = 8
rin = 0.001
rout = 0.01
v_in = 1.0e-4

# Get DG matrices from CADETJulia
nodes, weights = CADETJulia.DGElements.lglnodes(polyDeg)
D = CADETJulia.DGElements.PolynomialDerivativeMatrix(polyDeg, nodes)
M = CADETJulia.DGElements.MassMatrix(polyDeg, nodes, weights)

Δρ = (rout - rin) / nCells

println("DG Discretization:")
println("  Polynomial degree: $polyDeg")
println("  Number of cells: $nCells")
println("  Element size Δρ: $Δρ m")
println("  Number of nodes per element: $(polyDeg + 1)")
println()

println("Checking discrete continuity for each cell...")
println("-" ^ 80)

max_error = 0.0
max_error_cell = 0

for cell in 1:nCells
    # Physical coordinates of nodes in this cell
    rho_left = rin + (cell - 1) * Δρ
    rho_nodes = [rho_left + (Δρ / 2) * (1 + ξ) for ξ in nodes]

    # Velocity at nodes: v(ρ) = v_in * rin / ρ (from continuity)
    v_nodes = [v_in * rin / ρ for ρ in rho_nodes]

    # Constant concentration c = 1
    c_const = ones(polyDeg + 1)

    # Compute D^T * M * (v .* c)
    # This should be zero if discrete continuity holds
    vc = v_nodes .* c_const
    residual = D' * (M * vc)

    # Check magnitude
    max_res_cell = maximum(abs.(residual))

    if max_res_cell > max_error
        max_error = max_res_cell
        max_error_cell = cell
    end

    @printf("Cell %2d: max|D^T*M*(v*1)| = %.6e", cell, max_res_cell)
    if max_res_cell > 1e-12
        println("  ❌ FAIL")
    else
        println("  ✓ PASS")
    end
end

println("-" ^ 80)
println()

println("RESULT:")
println("  Maximum error across all cells: $max_error")
println("  Occurred in cell: $max_error_cell")
println()

if max_error > 1e-12
    println("❌ DISCRETE CONTINUITY TEST FAILED!")
    println()
    println("The discrete operator D^T * M does NOT satisfy geometric conservation.")
    println("This explains why free-stream preservation fails!")
    println()
    println("DIAGNOSIS:")
    println("  The standard DG stiffness matrix S = D^T * M is derived for")
    println("  Cartesian coordinates and does NOT account for the metric terms")
    println("  in cylindrical coordinates.")
    println()
    println("POSSIBLE FIXES:")
    println("  1. Use D^T * M_ρ instead of D^T * M for convection")
    println("     (where M_ρ includes radial weighting)")
    println("  2. Implement a flux-form discretization that enforces")
    println("     discrete conservation: ∂(ρv)/∂ρ = 0")
    println("  3. Add a correction term to enforce geometric conservation")
else
    println("✓ DISCRETE CONTINUITY TEST PASSED!")
    println()
    println("The discrete operator satisfies geometric conservation.")
    println("The bug must be elsewhere in the implementation.")
end

println("=" ^ 80)
