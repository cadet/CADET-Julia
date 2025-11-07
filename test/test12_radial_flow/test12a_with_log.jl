#!/usr/bin/env julia
# Free-Stream Preservation Test for Radial Flow WITH DEBUG LOGGING
#
# Tests that a constant state (c=1.0 everywhere) remains constant.
#
# Setup:
#   - Initial condition: c(r, t=0) = 1.0 everywhere
#   - Inlet condition:   cIn(t) = 1.0 for all time
#   - Expected result:   c(r,t) = 1.0 for all time (constant state preserved)
#
# Theory:
#   For constant state c=1.0:
#     - dc/dr = 0 (no spatial gradient)
#     - Convection term: d(v*c)/dr = d(constant)/dr = 0
#     - Dispersion term: d(rho*D*dc/dr)/dr = 0
#     - Therefore: dc/dt = 0, so c should remain constant!
#
# Any non-zero residual indicates a bug in the discretization.

using CADETJulia
using Printf

# Enable debug logging
CADETJulia.RadialConvDispOperatorDG.set_debug_log("test12a_debug.log")

println("=" ^ 70)
println("FREE-STREAM PRESERVATION TEST - RADIAL FLOW (WITH DEBUG LOG)")
println("=" ^ 70)
println()
println("Debug output will be written to: test12a_debug.log")
println()

# ==================== TEST PARAMETERS ====================
ncomp = 1           # Single component
c_const = 1.0       # Constant concentration everywhere
tend = 0.01         # Test for short time (only initial residual matters)

# Geometry
rin = 0.01         # Inner radius [m]
rout = 0.1         # Outer radius [m]

# Physical properties
D_rad = 1.0e-8      # Radial dispersion coefficient [m^2/s]
epsilon_c = 1.0     # Porosity (1.0 = no particles, pure mobile phase)

# Discretization
polyDeg = 4         # Polynomial degree
nCells = 8          # Number of cells

# Flow
u_in = 1.0e-4       # Inlet velocity [m/s]

# Solver tolerances
abstol = 1.0e-10
reltol = 1.0e-8
dtout = 0.01        # Output at end only

println("TEST CONFIGURATION:")
println("  Constant state test:")
println("    Initial c(r,t=0) = $c_const (everywhere)")
println("    Inlet cIn(t)     = $c_const (constant for all time)")
println("    Expected:        c(r,t) = $c_const (should remain constant)")
println()
println("  Geometry:")
println("    Inner radius = $rin m")
println("    Outer radius = $rout m")
println()
println("  Physics:")
println("    Dispersion D = $D_rad m^2/s")
println("    Porosity eps_c = $epsilon_c")
println("    Velocity u = $u_in m/s")
println()
println("  Discretization:")
println("    Polynomial degree = $polyDeg")
println("    Number of cells = $nCells")
println("    Total DOFs = $(nCells * (polyDeg + 1))")
println()
println("  Solver:")
println("    abstol = $abstol")
println("    reltol = $reltol")
println("    Simulation time = $tend s")
println("=" ^ 70)
println()

# ==================== BUILD MODEL ====================

# Column: Initialize with constant concentration c=1.0 everywhere
col = CADETJulia.rLRM(
    nComp = ncomp,
    col_inner_radius = rin,
    col_outer_radius = rout,
    d_rad = fill(D_rad, ncomp),
    eps_c = epsilon_c,
    c0 = fill(c_const, ncomp),  # CONSTANT initial condition
    q0 = fill(0.0, ncomp),
    polyDeg = polyDeg,
    nCells = nCells,
    cross_section_area = 1.0,
)

# Inlet: Single section with constant concentration
inlet = CADETJulia.CreateInlet(nComp = ncomp, nSections = 1)

# Section 1: Constant inlet (t = 0 -> tend)
CADETJulia.modify_inlet!(
    inlet = inlet,
    nComp = ncomp,
    section = 1,
    cIn_c = fill(c_const, ncomp),  # CONSTANT inlet c=1.0
    cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp),
    cIn_cube = zeros(ncomp),
)

# Outlet
outlet = CADETJulia.CreateOutlet(nComp = ncomp)

# ==================== SWITCHES & CONNECTIONS ====================
idx_units = [0]

switches = CADETJulia.Switches(
    nSections = 1,
    section_times = [0.0, tend],
    nSwitches = 1,
    nColumns = 1,
    nComp = ncomp,
    idx_units = idx_units,
)

# Connect INLET -> COLUMN
CADETJulia.Connection(switches, 1, 1, inlet, 1, u_in, 0.0, 0.0, 0.0, 0.0, false)

# Connect COLUMN -> OUTLET
CADETJulia.Connection(switches, 1, 1, 1, outlet, 0.0, 0.0, 0.0, 0.0, 0.0, false)

# ==================== SOLVER OPTIONS ====================
solverOptions = CADETJulia.SolverCache(
    columns = (col,), switches = switches,
    outlets = (outlet,), abstol = abstol,
    reltol = reltol, solution_times = collect(0.0:dtout:tend),
    prototypeJacobian = true, analyticalJacobian = false,
)

# ==================== INITIAL RESIDUAL CHECK ====================
println("=" ^ 70)
println("CHECKING INITIAL RESIDUAL (WITH DEBUG OUTPUT)")
println("=" ^ 70)
println("Computing residual at t=0 with constant state c=$c_const...")
println()

x0 = solverOptions.x0
dx0 = similar(x0)
fill!(dx0, 0.0)

# Set up parameter tuple (same as in solve_model)
RHS_q = col.RHS_q
cpp = col.cpp
qq = col.qq
i = 1  # section 1
p = (columns=(col,), RHS_q, cpp, qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, nothing)

# IMPORTANT: For StaticInlets, we must manually set cIn
for j = 1:ncomp
    col.cIn[j] = switches.ConnectionInstance.cIn_c[1][1][j][1]
end
println("=> Set col.cIn = $(col.cIn)")
println()

# Check initial state
println("Initial state check:")
@printf("  x0[1:5] = [%.6f, %.6f, %.6f, %.6f, %.6f]\n", x0[1], x0[2], x0[3], x0[4], x0[5])
@printf("  Expected: all values = %.6f\n", c_const)
println()

# Compute residual at t=0 - THIS WILL GENERATE DEBUG OUTPUT
println("Computing residual (debug output going to log file)...")
CADETJulia.problem!(dx0, x0, p, 0.0)

# Check residual
max_init_res = maximum(abs.(dx0))

println()
println("INITIAL RESIDUAL ANALYSIS:")
println("-" ^ 70)
@printf("  Maximum |residual| = %.6e\n", max_init_res)
@printf("  Mean |residual|    = %.6e\n", sum(abs.(dx0))/length(dx0))
@printf("  State vector size  = %d DOFs\n", length(x0))
println()

# Pass/fail criterion
TOLERANCE = 1e-10
if max_init_res > TOLERANCE
    println("[FAIL] INITIAL RESIDUAL CHECK FAILED!")
    println("   Maximum |residual| = $max_init_res exceeds tolerance $TOLERANCE")
    println()

    # Find where the largest residuals are
    idx_max = argmax(abs.(dx0))
    println("   Location of max residual: DOF index $idx_max")
    @printf("   Residual value: %.6e\n", dx0[idx_max])
    @printf("   State value:    %.6e\n", x0[idx_max])
    println()

    # Compute which cell and node
    nNodes = polyDeg + 1
    cell_num = div(idx_max - 1, nNodes) + 1
    node_num = mod(idx_max - 1, nNodes) + 1
    println("   This corresponds to: Cell $cell_num, Node $node_num")
    println()

    # Show all residuals
    println("   All residual values:")
    for j in 1:length(dx0)
        cell = div(j - 1, nNodes) + 1
        node = mod(j - 1, nNodes) + 1
        @printf("     DOF %2d (Cell %2d, Node %2d):  x = %.6e,  dx = %.6e\n",
                j, cell, node, x0[j], dx0[j])
    end
    println()

    println("   => BUG DETECTED in discretization!")
    println("   => Constant state should have zero residual.")
    println("   => Check test12a_debug.log for detailed diagnostic output")
else
    println("[PASS] INITIAL RESIDUAL CHECK PASSED!")
    println("   Maximum |residual| = $max_init_res < $TOLERANCE")
    println("   The constant state has zero time derivative (as expected).")
    println("   => Discretization is consistent!")
end
println("=" ^ 70)
println()

# Close debug log
CADETJulia.RadialConvDispOperatorDG.close_debug_log()
println("Debug log closed. Check test12a_debug.log for detailed output.")
println()
println("=" ^ 70)
println("TEST COMPLETE")
println("=" ^ 70)
