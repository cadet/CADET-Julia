#!/usr/bin/env julia
# Linear Concentration Residual Test
# Test that when c = ρ everywhere, what is the residual
# This verifies the radial convection-dispersion operator for linear profile

using CADETJulia
using Printf
using LinearAlgebra

println("=" ^ 80)
println("LINEAR CONCENTRATION (c = ρ) RESIDUAL TEST")
println("=" ^ 80)
println()
println("This test verifies the residual when c = ρ (concentration equals radius).")
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

println("Radial positions:")
@printf("  ρ_min (inlet):  %.6f m\n", rho_nodes[1])
@printf("  ρ_max (outlet): %.6f m\n", rho_nodes[end])
println()

# Create column model with c0 = ρ
col = CADETJulia.rLRM(
    nComp = ncomp,
    col_inner_radius = rin,
    col_outer_radius = rout,
    d_rad = fill(D_rad, ncomp),
    eps_c = epsilon_c,
    c0 = rho_nodes,  # Set initial concentration to ρ at each point
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

# Create state vector with c = ρ everywhere
# State vector has structure: [c_1, c_2, ..., c_n, q_1, q_2, ..., q_n]
x = zeros(Float64, col.unitStride)

# Set concentration part to c = ρ
for i in 1:nPoints
    x[i] = rho_nodes[i]
end

# Adsorbed phase remains zero (indices nPoints+1 to end)

# Set inlet concentration to match c = ρ at inlet
col.cIn[1] = rin

println("Test State:")
println("  State vector: c(ρ) = ρ for mobile phase")
println("  Adsorbed phase: q = 0 everywhere")
@printf("  Inlet concentration: c_in = %.6f (matching ρ_inlet)\n", rin)
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
    cIn_c = fill(rin, ncomp),  # Match c = ρ at inlet
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
    @printf("  (spatial location ρ = %.6f m, c = %.6f)\n", rho_nodes[idx_max], x[idx_max])
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
        @printf("  Index %4d (ρ = %.6f m, c = %.6f): %.6e\n",
                idx, rho_nodes[idx], x[idx], Dc[idx])
    end
end
println()

# ==================== PHYSICAL INTERPRETATION ====================

println("=" ^ 80)
println("PHYSICAL INTERPRETATION")
println("=" ^ 80)
println()

println("For the linear profile c = ρ:")
println()
println("The governing equation in radial coordinates is:")
println("  ∂c/∂t = -v ∂c/∂ρ + (1/ρ) ∂/∂ρ(ρ D ∂c/∂ρ)")
println()
println("For c = ρ (steady state, ∂c/∂t = 0):")
println("  ∂c/∂ρ = 1")
println("  ∂²c/∂ρ² = 0")
println()
println("Therefore:")
println("  -v ∂c/∂ρ + (1/ρ) ∂/∂ρ(ρ D ∂c/∂ρ)")
println("  = -v(1) + (1/ρ) ∂/∂ρ(ρ D × 1)")
println("  = -v + (1/ρ) ∂/∂ρ(ρ D)")
println("  = -v + (1/ρ)(D + ρ × 0)")
println("  = -v + D/ρ")
println()
@printf("Expected residual: -v + D/ρ = -%.6e + %.6e/ρ\n", u_in, D_rad)
println()
println("At different radial positions:")
for idx in sample_indices[1:min(5, length(sample_indices))]
    if idx <= nPoints
        expected_res = -u_in + D_rad / rho_nodes[idx]
        @printf("  ρ = %.6f m: Expected = %.6e, Computed = %.6e, Diff = %.6e\n",
                rho_nodes[idx], expected_res, Dc[idx], abs(Dc[idx] - expected_res))
    end
end

println()
println("=" ^ 80)
println("CONCLUSION")
println("=" ^ 80)
println()

println("For c = ρ, the residual is NOT zero (as expected).")
println("The residual should equal: -v + D/ρ at each point.")
@printf("This gives a residual of order %.6e to %.6e\n",
        -u_in + D_rad/rho_nodes[end], -u_in + D_rad/rho_nodes[1])
println()
println("=" ^ 80)
