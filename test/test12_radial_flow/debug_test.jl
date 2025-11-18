#!/usr/bin/env julia
# Debug test to understand what the radial operator computes

using CADETJulia
using Printf
using CADETJulia: RadialConvDispOperatorDG

# Simple test case
ncomp = 1
rin = 0.01
rout = 0.1
D_rad = 1.0e-8
u_in = 1.0e-4
epsilon_c = 1.0

# Small test: 4 cells, polynomial degree 2
nCells = 4
polyDeg = 2
nNodes = polyDeg + 1

println("=" ^ 70)
println("DEBUG: Radial Operator Test")
println("=" ^ 70)
println("nCells = $nCells, polyDeg = $polyDeg")
println()

# Create column
col = CADETJulia.rLRM(
    nComp = ncomp,
    col_inner_radius = rin,
    col_outer_radius = rout,
    d_rad = fill(D_rad, ncomp),
    eps_c = epsilon_c,
    c0 = fill(0.0, ncomp),
    q0 = fill(0.0, ncomp),
    polyDeg = polyDeg,
    nCells = nCells,
    cross_section_area = 1.0,
)

op = col.ConvDispOpInstance
nPoints = nCells * nNodes

# Compute nodal positions
r_nodes = zeros(nPoints)
for cell in 1:nCells
    for node in 1:nNodes
        idx = (cell - 1) * nNodes + node
        r_nodes[idx] = rin + (cell - 1) * op.deltarho + (op.deltarho / 2) * (1 + op.nodes[node])
    end
end

println("Nodal positions (r):")
for i in 1:nPoints
    @printf("  Node %2d: r = %.6f\n", i, r_nodes[i])
end
println()

# Test with manufactured solution c = r²
c_exact(r) = r^2
x_state = c_exact.(r_nodes)

println("State values (c = r²):")
for i in 1:nPoints
    @printf("  Node %2d: c = %.6f\n", i, x_state[i])
end
println()

# Analytical residual: R(r) = -v_in·r_in/r + 4D
analytical_residual(r) = -u_in * rin / r + 4.0 * D_rad
r_analytical = analytical_residual.(r_nodes)

println("Analytical residual:")
for i in 1:nPoints
    @printf("  Node %2d: R_analytical = %.6e\n", i, r_analytical[i])
end
println()

# Apply operator
fill!(op.Dc, 0.0)
fill!(op.Dg, 0.0)
fill!(op.h, 0.0)
fill!(op.c_star, 0.0)
fill!(op.g_star, 0.0)

idx = 1:nPoints
RadialConvDispOperatorDG.radialresidualImpl!(
    op.Dc, x_state, idx, op.strideNode, op.strideCell,
    op.nNodes, nCells, op.deltarho, polyDeg,
    op.polyDerM, op.invMM, op.S, op.rMM, op.invrMM,
    op.nodes, op.invWeights, u_in, D_rad, op.rho_i,
    c_exact(rin), op.c_star, op.g_star, op.Dg, op.h, op.mul1
)

println("Numerical residual (op.Dc):")
for i in 1:nPoints
    @printf("  Node %2d: R_numerical = %.6e\n", i, op.Dc[i])
end
println()

println("Error (R_numerical - R_analytical):")
error_vec = op.Dc .- r_analytical
for i in 1:nPoints
    @printf("  Node %2d: error = %.6e\n", i, error_vec[i])
end
println()

println("Max error: ", maximum(abs.(error_vec)))
println()

# Also print intermediate values
println("Intermediate values:")
println("  c_star (interface fluxes):")
for i in 1:nCells+1
    @printf("    Interface %d: %.6e\n", i, op.c_star[i])
end
println()

println("  g_star (gradient fluxes):")
for i in 1:nCells+1
    @printf("    Interface %d: %.6e\n", i, op.g_star[i])
end
println()

println("  h (auxiliary variable, scaled gradient):")
for i in 1:nPoints
    @printf("    Node %2d: h = %.6e\n", i, op.h[i])
end
println()

println("=" ^ 70)
