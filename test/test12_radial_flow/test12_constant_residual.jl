#!/usr/bin/env julia
# Constant Concentration Residual Test
# Test that when c = 1 everywhere, the residual should be zero
# This verifies the radial convection-dispersion operator implementation

using CADETJulia
using Printf
using LinearAlgebra

println("=" ^ 80)
println("CONSTANT CONCENTRATION RESIDUAL TEST")
println("=" ^ 80)
println()
println("This test verifies that when c = 1 everywhere (constant),")
println("the residual from the radial operator should be zero.")
println()

# ==================== TEST PARAMETERS ====================
ncomp = 1

# Geometry
rin = 0.01
rout = 0.1

# Physical properties
D_rad = 1.0e-4
epsilon_c = 1.0

# Discretization
polyDeg = 3
nCells = 10

# Flow
Q = 1.0e-2

println("Test Parameters:")
println("  Inner radius: $rin m")
println("  Outer radius: $rout m")
println("  Radial dispersion: $D_rad m²/s")
println("  Polynomial degree: $polyDeg")
println("  Number of cells: $nCells")
println("  Flow rate Q: $Q m³/s")
println()

# ==================== BUILD MODEL ====================

using CADETJulia.DGElements

nodes_ref, _ = DGElements.cglnodes(polyDeg)

deltarho = (rout - rin) / nCells
nNodes = polyDeg + 1
nPoints = nCells * nNodes

println("Discretization:")
println("  Number of nodes per cell: $nNodes")
println("  Total spatial points: $nPoints")
println("  Cell size Δρ: $deltarho m")
println()

# Compute nodal positions
rho_nodes = zeros(nPoints)
for cell in 1:nCells
    for node in 1:nNodes
        idx = (cell - 1) * nNodes + node
        rho_nodes[idx] = rin + (cell - 1) * deltarho + (deltarho / 2) * (1 + nodes_ref[node])
    end
end

# Create column model
col = CADETJulia.rLRM(
    nComp = ncomp,
    col_inner_radius = rin,
    col_outer_radius = rout,
    d_rad = fill(D_rad, ncomp),
    eps_c = epsilon_c,
    c0 = fill(1.0, ncomp),  # Set initial concentration to 1.0 everywhere
    q0 = fill(0.0, ncomp),
    polyDeg = polyDeg,
    nCells = nCells,
)

u_in = Q / (col.cross_section_area * col.eps_c)

println("Flow:")
println("  Cross-sectional area: $(col.cross_section_area) m²")
println("  Interstitial velocity u_in: $u_in m/s")
println()

# ==================== SETUP TEST STATE ====================

# Create state vector with c = 1.0 everywhere
# State vector has structure: [c_1, c_2, ..., c_n, q_1, q_2, ..., q_n]
x = zeros(Float64, col.unitStride)

# Set concentration part to c = 1.0
for i in 1:nPoints
    x[i] = 1.0
end

# Adsorbed phase remains zero (indices nPoints+1 to end)

# Set inlet concentration to 1.0 (matching the constant state)
col.cIn[1] = 1.0

println("Test State:")
println("  State vector: c = 1.0 for mobile phase")
println("  Adsorbed phase: q = 0 everywhere")
println("  Inlet concentration: c_in = 1.0")
println("  Total DOFs: $(col.unitStride)")
println("  Concentration DOFs: $nPoints")
println("  Adsorbed phase DOFs: $(col.unitStride - nPoints)")
println()

# ==================== COMPUTE RESIDUAL ====================

println("=" ^ 80)
println("COMPUTING RESIDUAL")
println("=" ^ 80)
println()

# Allocate residual vector
Dc = zeros(Float64, col.unitStride)

# Create dummy RHS_q (not used for this test)
RHS_q = [zeros(Float64, ncomp) for _ in 1:1]

# Setup parameters needed for residual computation
cpp = col.cpp
qq = col.qq
i = 1
t = 0.0  # Time doesn't matter for this test

# Create inlet and switches for completeness
inlet = CADETJulia.CreateInlet(nComp = ncomp, nSections = 1)
CADETJulia.modify_inlet!(
    inlet = inlet,
    nComp = ncomp,
    section = 1,
    cIn_c = fill(1.0, ncomp),
    cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp),
    cIn_cube = zeros(ncomp),
)

outlet = CADETJulia.CreateOutlet(nComp = ncomp)

idx_units = [0]
switches = CADETJulia.Switches(
    nSections = 1,
    section_times = [0.0, 1.0],
    nSwitches = 1,
    nColumns = 1,
    nComp = ncomp,
    idx_units = idx_units,
)

CADETJulia.Connection(switches, 1, 1, inlet, 1, u_in, 0.0, 0.0, 0.0, 0.0, false)
CADETJulia.Connection(switches, 1, 1, 1, outlet, 0.0, 0.0, 0.0, 0.0, 0.0, false)

# Call the residual function
p = ((col,), RHS_q, cpp, qq, i, 1, idx_units, switches, nothing)

CADETJulia.problem!(Dc, x, p, t)

println("Residual computed!")
println()

# ==================== ANALYZE RESIDUAL ====================

println("=" ^ 80)
println("RESIDUAL ANALYSIS")
println("=" ^ 80)
println()

# Compute norms
residual_norm_inf = norm(Dc, Inf)
residual_norm_2 = norm(Dc, 2)
residual_norm_1 = norm(Dc, 1)

println("Residual Norms:")
@printf("  ||Dc||_∞ (max absolute value): %.6e\n", residual_norm_inf)
@printf("  ||Dc||_2 (Euclidean norm):      %.6e\n", residual_norm_2)
@printf("  ||Dc||_1 (sum of abs values):   %.6e\n", residual_norm_1)
println()

# Find location of maximum residual
idx_max = argmax(abs.(Dc))
@printf("Maximum residual at index %d: %.6e\n", idx_max, Dc[idx_max])
if idx_max <= nPoints
    @printf("  (spatial location ρ = %.6f m)\n", rho_nodes[idx_max])
end
println()

# Show statistics
println("Residual Statistics:")
@printf("  Min: %.6e\n", minimum(Dc))
@printf("  Max: %.6e\n", maximum(Dc))
@printf("  Mean: %.6e\n", sum(Dc) / length(Dc))
@printf("  Std Dev: %.6e\n", sqrt(sum((Dc .- sum(Dc)/length(Dc)).^2) / length(Dc)))
println()

# Sample residual values at different locations
println("Sample Residual Values:")
sample_indices = [1, div(nPoints, 4), div(nPoints, 2), div(3*nPoints, 4), nPoints]
for idx in sample_indices
    if idx <= length(Dc)
        @printf("  Index %4d (ρ = %.6f m): %.6e\n", idx, rho_nodes[idx], Dc[idx])
    end
end
println()

# ==================== TEST RESULT ====================

println("=" ^ 80)
println("TEST RESULT")
println("=" ^ 80)
println()

tolerance = 1.0e-12
passed = residual_norm_inf < tolerance

if passed
    println("✓ TEST PASSED!")
    @printf("  The residual is numerically zero (||Dc||_∞ = %.6e < %.6e)\n", residual_norm_inf, tolerance)
    println("  The radial operator correctly handles constant concentration.")
else
    println("✗ TEST FAILED!")
    @printf("  The residual is NOT zero (||Dc||_∞ = %.6e >= %.6e)\n", residual_norm_inf, tolerance)
    println("  There may be an issue with the radial operator implementation.")
    println()
    println("  Largest residual components (top 10):")
    residual_sorted_idx = sortperm(abs.(Dc), rev=true)
    for j in 1:min(10, length(Dc))
        idx = residual_sorted_idx[j]
        if idx <= nPoints
            @printf("    Index %4d (ρ = %.6f m): %.6e\n", idx, rho_nodes[idx], Dc[idx])
        else
            @printf("    Index %4d: %.6e\n", idx, Dc[idx])
        end
    end
end

println()
println("=" ^ 80)
