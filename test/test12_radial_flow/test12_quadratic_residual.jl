#!/usr/bin/env julia
# Quadratic Concentration Residual Test
# Test that when c = ρ² everywhere, what is the residual
# This verifies the radial convection-dispersion operator for quadratic profile

using CADETJulia
using Printf
using LinearAlgebra

println("=" ^ 80)
println("QUADRATIC CONCENTRATION (c = ρ²) RESIDUAL TEST")
println("=" ^ 80)
println()
println("This test verifies the residual when c = ρ² (quadratic profile).")
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

# Create column model with c0 = ρ²
col = CADETJulia.rLRM(
    nComp = ncomp,
    col_inner_radius = rin,
    col_outer_radius = rout,
    d_rad = fill(D_rad, ncomp),
    eps_c = epsilon_c,
    c0 = rho_nodes.^2,  # Set initial concentration to ρ² at each point
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

# Create state vector with c = ρ² everywhere
# State vector has structure: [c_1, c_2, ..., c_n, q_1, q_2, ..., q_n]
x = zeros(Float64, col.unitStride)

# Set concentration part to c = ρ²
for i in 1:nPoints
    x[i] = rho_nodes[i]^2
end

# Adsorbed phase remains zero (indices nPoints+1 to end)

# Set inlet concentration to match c = ρ² at inlet
col.cIn[1] = rin^2

println("Test State:")
println("  State vector: c(ρ) = ρ²")
println("  Adsorbed phase: q = 0 everywhere")
@printf("  Inlet concentration: c_in = %.6e (matching ρ_inlet²)\n", rin^2)
@printf("  Outlet concentration: c_out = %.6e (ρ_outlet²)\n", rho_nodes[end]^2)
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
    cIn_c = fill(rin^2, ncomp),  # Match c = ρ² at inlet
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
    @printf("  (spatial location ρ = %.6f m, c = %.6e)\n", rho_nodes[idx_max], x[idx_max])
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
    if idx <= length(Dc) && idx <= nPoints
        @printf("  Index %4d (ρ = %.6f m, c = %.6e): %.6e\n",
                idx, rho_nodes[idx], x[idx], Dc[idx])
    end
end
println()

# ==================== PHYSICAL INTERPRETATION ====================

println("=" ^ 80)
println("PHYSICAL INTERPRETATION")
println("=" ^ 80)
println()

println("For the quadratic profile c = ρ²:")
println()
println("The governing equation in radial coordinates is:")
println("  ∂c/∂t = -v(ρ) ∂c/∂ρ + (1/ρ) ∂/∂ρ(ρ D ∂c/∂ρ)")
println()
println("where v(ρ) = v_in × ρ_in / ρ (due to area conservation)")
println()
println("For c = ρ²:")
println("  ∂c/∂ρ = 2ρ")
println("  ∂²c/∂ρ² = 2")
println()
println("Convection term:")
println("  -v(ρ) ∂c/∂ρ = -(v_in × ρ_in / ρ) × 2ρ = -2 v_in × ρ_in")
println()
println("Diffusion term:")
println("  (1/ρ) ∂/∂ρ(ρ D ∂c/∂ρ) = (1/ρ) ∂/∂ρ(ρ D × 2ρ)")
println("                         = (1/ρ) ∂/∂ρ(2D ρ²)")
println("                         = (1/ρ) × 4D ρ")
println("                         = 4D")
println()
println("Total residual (steady state):")
println("  Residual = -2 v_in × ρ_in + 4D")
println()

v_rin = u_in * rin
expected_residual = -2.0 * v_rin + 4.0 * D_rad

@printf("  = -2 × %.6e × %.6f + 4 × %.6e\n", u_in, rin, D_rad)
@printf("  = %.6e\n", expected_residual)
println()

println("Since the residual is constant everywhere (independent of ρ),")
println("we expect the same residual value at all spatial points.")
println()

println("Comparison with computed residuals:")
for idx in sample_indices[1:min(5, length(sample_indices))]
    if idx <= nPoints
        @printf("  ρ = %.6f m: Expected = %.6e, Computed = %.6e, Diff = %.6e\n",
                rho_nodes[idx], expected_residual, Dc[idx], abs(Dc[idx] - expected_residual))
    end
end

println()
println("=" ^ 80)
println("CONCLUSION")
println("=" ^ 80)
println()

if abs(residual_norm_inf - abs(expected_residual)) < 1.0e-10
    println("✓ The computed residual matches the theoretical prediction!")
    @printf("  Expected: %.6e\n", expected_residual)
    @printf("  Computed: %.6e\n", residual_norm_inf)
else
    println("The computed residual differs from the simple theoretical prediction.")
    println("This may indicate:")
    println("  - Additional terms in the discretization")
    println("  - Numerical artifacts from the DG method")
    println("  - Boundary condition effects")
end

println()
println("=" ^ 80)
