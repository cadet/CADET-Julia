# Pulse Injection Test for Radial Flow
#
# Simple test: Inject c=1.0 for 1 second into an empty column
#
# Timeline:
#   t = 0→1s: Inject c=1.0 (pulse)
#   t = 1→2s: Inject c=0.0 (wash)

using CADETJulia
using Printf

println("=" ^ 70)
println("PULSE INJECTION TEST - RADIAL FLOW")
println("=" ^ 70)
println()

# ==================== TEST PARAMETERS ====================
ncomp  = 1          # Single component
c_pulse = 1.0       # Pulse concentration
tinj = 1.0          # Injection duration: 0 → 1 second
tend = 1.1          # Total time: 0 → 2 seconds

# Geometry
rin = 0.01         # Inner radius [m]
rout = 0.2         # Outer radius [m]

# Physical properties
# Define radial dispersion coefficient as a function of radius r [m]
# Example: D_rad(r) = D0 * (1 + α*r) - increases linearly with radius
D0 = 1.0e-8         # Base dispersion coefficient [m²/s]
α = 2.0             # Scaling factor [1/m]
D_rad_func(r) = D0 * (1.0 + α * r)  # Variable dispersion coefficient
# Alternatively, use constant: D_rad = 1.0e-8
εc = 1.0            # Porosity (1.0 = no particles, pure mobile phase)

# Display dispersion coefficient information
println("Variable Dispersion Coefficient:")
println("  D_rad(r) = D0 * (1 + α*r)")
println("  D0 = $D0 m²/s")
println("  α = $α 1/m")
println()
@printf("  D_rad(r_inner=%.4f m) = %.6e m²/s\n", rin, D_rad_func(rin))
@printf("  D_rad(r_mid=%.4f m) = %.6e m²/s\n", (rin+rout)/2, D_rad_func((rin+rout)/2))
@printf("  D_rad(r_outer=%.4f m) = %.6e m²/s\n", rout, D_rad_func(rout))
println()

# Discretization
polyDeg = 4         # Polynomial degree
nCells = 16          # Number of cells

# Flow
u_in = 1.0e-5       # Inlet velocity [m/s]

# Solver tolerances
abstol = 1.0e-10
reltol = 1.0e-8
Δtout = 0.1         # Output every 0.1 seconds

# ==================== BUILD MODEL ====================

# Column: Start with empty column (c=0 everywhere)
# Use variable dispersion coefficient function
col = CADETJulia.rLRM(
    nComp = ncomp,
    col_inner_radius = rin,
    col_outer_radius = rout,
    d_rad = D_rad_func,  # Pass function for variable dispersion
    eps_c = εc,
    c0 = fill(0.0, ncomp),  # Empty initial condition
    q0 = fill(0.0, ncomp),
    polyDeg = polyDeg,
    nCells = nCells,
    cross_section_area = 1.0,
)

# Inlet: 2 sections (pulse, then wash)
inlet = CADETJulia.CreateInlet(nComp = ncomp, nSections = 2)

# Section 1: Pulse (t = 0 → tinj)
CADETJulia.modify_inlet!(
    inlet = inlet,
    nComp = ncomp,
    section = 1,
    cIn_c = fill(c_pulse, ncomp),  # Inject c=1.0
    cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp),
    cIn_cube = zeros(ncomp),
)

# Section 2: Wash (t = tinj → tend)
CADETJulia.modify_inlet!(
    inlet = inlet,
    nComp = ncomp,
    section = 2,
    cIn_c = fill(0.0, ncomp),  # Inject c=0.0
    cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp),
    cIn_cube = zeros(ncomp),
)

# Outlet
outlet = CADETJulia.CreateOutlet(nComp = ncomp)

# ==================== SWITCHES & CONNECTIONS ====================
idx_units = [0]

switches = CADETJulia.Switches(
    nSections = 2,
    section_times = [0.0, tinj, tend],
    nSwitches = 1,
    nColumns = 1,
    nComp = ncomp,
    idx_units = idx_units,
)

# Connect INLET → COLUMN
CADETJulia.Connection(
    switches, 1, 1, inlet, 1,
    u_in, 0.0, 0.0, 0.0, 0.0, false
)

# Connect COLUMN → OUTLET
CADETJulia.Connection(
    switches, 1, 1, 1, outlet,
    0.0, 0.0, 0.0, 0.0, 0.0, false
)

# ==================== SOLVER OPTIONS ====================
solverOptions = CADETJulia.SolverCache(
    columns = (col,),
    switches = switches,
    outlets = (outlet,),
    abstol = abstol,
    reltol = reltol,
    solution_times = collect(0.0:Δtout:tend),
    prototypeJacobian = true,
    analyticalJacobian = false,
)

# ==================== SOLVE ====================
println("=" ^ 70)
println("RUNNING SIMULATION")
println("=" ^ 70)
println("Solving from t=0 to t=$tend seconds...")
println()

sol = CADETJulia.solve_model(
    columns = (col,),
    switches = switches,
    outlets = (outlet,),
    solverOptions = solverOptions
)

println("✓ Simulation completed!")

println("=" ^ 70)
println()

# ==================== RESULTS ====================
println("=" ^ 70)
println("OUTLET CONCENTRATION")
println("=" ^ 70)

if !isempty(col.solution_times)
    println("\nOutlet concentration vs time:")
    println("-" ^ 70)
    @printf("%10s  %20s\n", "Time [s]", "c_out")
    println("-" ^ 70)

    for k in 1:length(col.solution_times)
        t = col.solution_times[k]
        c = col.solution_outlet[k, 1]
        @printf("%10.3f  %20.10e\n", t, c)
    end

    println("-" ^ 70)

    # Summary statistics
    c_max = maximum(col.solution_outlet[:, 1])
    c_min = minimum(col.solution_outlet[:, 1])
    c_mean = sum(col.solution_outlet[:, 1]) / length(col.solution_outlet[:, 1])

    println("\nSummary:")
    @printf("  Maximum outlet concentration: %.6e\n", c_max)
    @printf("  Minimum outlet concentration: %.6e\n", c_min)
    @printf("  Mean outlet concentration:    %.6e\n", c_mean)
    @printf("  Expected pulse height:        %.6e\n", c_pulse)

    if abs(c_max - c_pulse) > 0.1 * c_pulse
        println("\n⚠️  WARNING: Outlet concentration differs from expected!")
        println("   This may indicate a problem with the discretization.")
    end
else
    println("❌ ERROR: No solution data available!")
end

println()
println("=" ^ 70)
println("TEST COMPLETE")
println("=" ^ 70)
