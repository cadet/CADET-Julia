using Test, CADETJulia

# ==================== TEST DESCRIPTION ====================
# Test file for radial LRMP (rLRMP) model with step pulse injection
# Uses variable dispersion D(ρ) and film diffusion k_f(ρ) coefficients
# Parameters: polyDeg=5, nCells=128

# ==================== PHYSICAL PARAMETERS ====================
nComp = 1
rin = 0.1           # Inner radius [m]
rout = 1.1          # Outer radius [m]
col_height = 1.0    # Column height [m]
v = 1.0e-2          # Interstitial velocity [m/s]
eps_c = 0.6         # Column porosity [-]
eps_p = 0.5         # Particle porosity [-]
Rp = 1.0e-4         # Particle radius [m]

# Variable dispersion coefficient D(ρ) - linear variation
# D(ρ) = D0 + D_slope * (ρ - rin)
D0 = 1.0e-4         # Base dispersion coefficient [m²/s]
D_slope = 5.0e-5    # Linear slope [m/s]

# Variable film diffusion coefficient k_f(ρ) - linear variation
# k_f(ρ) = kf0 + kf_slope * (ρ - rin)
kf0 = 1.0e-3        # Base film diffusion coefficient [m/s]
kf_slope = 5.0e-4   # Linear slope [1/s]

# ==================== SIMULATION PARAMETERS ====================
t_inject = 60.0     # Injection duration [s]
t_final = 200.0     # Final simulation time [s]

# ==================== DISCRETIZATION ====================
polyDeg = 5
nCells = 128

# ==================== BUILD MODEL DICTIONARY ====================
model = OrderedDict(
    "root" => OrderedDict(
        "input" => OrderedDict(
            "model" => OrderedDict()
        )
    )
)

# Unit 000: Inlet
model["root"]["input"]["model"]["unit_000"] = OrderedDict()
model["root"]["input"]["model"]["unit_000"]["unit_type"] = "INLET"
model["root"]["input"]["model"]["unit_000"]["ncomp"] = nComp
model["root"]["input"]["model"]["unit_000"]["inlet_type"] = "PIECEWISE_CUBIC_POLY"

# Section 0: Injection (c = 1)
model["root"]["input"]["model"]["unit_000"]["sec_000"] = OrderedDict()
model["root"]["input"]["model"]["unit_000"]["sec_000"]["const_coeff"] = [1.0]

# Section 1: Wash (c = 0)
model["root"]["input"]["model"]["unit_000"]["sec_001"] = OrderedDict()
model["root"]["input"]["model"]["unit_000"]["sec_001"]["const_coeff"] = [0.0]

# Unit 001: Radial LRMP Column
model["root"]["input"]["model"]["unit_001"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["unit_type"] = "RADIAL_LUMPED_RATE_MODEL_WITH_PORES"
model["root"]["input"]["model"]["unit_001"]["ncomp"] = nComp
model["root"]["input"]["model"]["unit_001"]["col_inner_radius"] = rin
model["root"]["input"]["model"]["unit_001"]["col_outer_radius"] = rout
model["root"]["input"]["model"]["unit_001"]["col_height"] = col_height
model["root"]["input"]["model"]["unit_001"]["col_porosity"] = eps_c
model["root"]["input"]["model"]["unit_001"]["par_porosity"] = eps_p
model["root"]["input"]["model"]["unit_001"]["par_radius"] = Rp
# Variable coefficients using special dictionary format
model["root"]["input"]["model"]["unit_001"]["col_dispersion"] = OrderedDict(
    "type" => "linear",
    "base" => D0,
    "slope" => D_slope
)
model["root"]["input"]["model"]["unit_001"]["film_diffusion"] = OrderedDict(
    "type" => "linear",
    "base" => kf0,
    "slope" => kf_slope
)
model["root"]["input"]["model"]["unit_001"]["velocity"] = v
model["root"]["input"]["model"]["unit_001"]["adsorption_model"] = "LINEAR"

model["root"]["input"]["model"]["unit_001"]["adsorption"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["adsorption"]["is_kinetic"] = true
model["root"]["input"]["model"]["unit_001"]["adsorption"]["LIN_KA"] = [0.0]
model["root"]["input"]["model"]["unit_001"]["adsorption"]["LIN_KD"] = [1.0]

model["root"]["input"]["model"]["unit_001"]["init_c"] = [0.0]
model["root"]["input"]["model"]["unit_001"]["init_cp"] = [0.0]
model["root"]["input"]["model"]["unit_001"]["init_q"] = [0.0]

model["root"]["input"]["model"]["unit_001"]["discretization"] = OrderedDict()
model["root"]["input"]["model"]["unit_001"]["discretization"]["polyDeg"] = polyDeg
model["root"]["input"]["model"]["unit_001"]["discretization"]["ncol"] = nCells
model["root"]["input"]["model"]["unit_001"]["discretization"]["nbound"] = ones(Bool, nComp)

# Unit 002: Outlet
model["root"]["input"]["model"]["unit_002"] = OrderedDict()
model["root"]["input"]["model"]["unit_002"]["unit_type"] = "OUTLET"
model["root"]["input"]["model"]["unit_002"]["ncomp"] = nComp

# Solver settings - 2 sections: injection + wash
model["root"]["input"]["solver"] = OrderedDict("sections" => OrderedDict())
model["root"]["input"]["solver"]["sections"]["nsec"] = 2
model["root"]["input"]["solver"]["sections"]["section_times"] = [0.0, t_inject, t_final]
model["root"]["input"]["solver"]["sections"]["section_continuity"] = [0]

# Connections
model["root"]["input"]["model"]["connections"] = OrderedDict()
model["root"]["input"]["model"]["connections"]["nswitches"] = 1
model["root"]["input"]["model"]["connections"]["switch_000"] = OrderedDict()
model["root"]["input"]["model"]["connections"]["switch_000"]["section"] = 0
model["root"]["input"]["model"]["connections"]["switch_000"]["connections"] = [0, 1, -1, -1, v,
                                                                               1, 2, -1, -1, v]

# User solution times
model["root"]["input"]["solver"]["user_solution_times"] = LinRange(0, t_final, Int(t_final) + 1)

# Time integrator settings
model["root"]["input"]["solver"]["time_integrator"] = OrderedDict()
model["root"]["input"]["solver"]["time_integrator"]["abstol"] = 1e-10
model["root"]["input"]["solver"]["time_integrator"]["algtol"] = 1e-8
model["root"]["input"]["solver"]["time_integrator"]["reltol"] = 1e-8

# ==================== CREATE UNITS ====================
# Variable coefficients are automatically parsed from the dictionary format
inlets, outlets, columns, switches, solverOptions = create_units(model)

# ==================== SOLVE MODEL ====================
solve_model(
    columns = columns,
    switches = switches,
    solverOptions = solverOptions,
    outlets = outlets,
    alg = QNDF(autodiff=AutoFiniteDiff()),
)

# ==================== OUTPUT RESULTS ====================
println("rLRMP Test with Variable D(ρ) and k_f(ρ)")
println("=========================================")
println("Parameters:")
println("  Inner radius: $rin m")
println("  Outer radius: $rout m")
println("  D(ρ) = $D0 + $D_slope * (ρ - $rin)")
println("  k_f(ρ) = $kf0 + $kf_slope * (ρ - $rin)")
println("  polyDeg: $polyDeg, nCells: $nCells")
println()

# Get solution
c_outlet = columns[1].solution_outlet[:, 1]
t_sol = columns[1].solution_times

println("Simulation completed successfully!")
println("  Final time: $(t_sol[end]) s")
println("  Max outlet concentration: $(maximum(c_outlet))")
println("  Time of max concentration: $(t_sol[argmax(c_outlet)]) s")
println("  Outlet concentration at t_final: $(c_outlet[end])")

# ==================== TESTS ====================
@test length(c_outlet) > 0
@test maximum(c_outlet) > 0  # Should have some concentration at outlet
@test c_outlet[1] ≈ 0.0 atol=1e-10  # Initial concentration should be zero

# Physical radial coordinates are now available from the model
rho_nodes = columns[1].ConvDispOpInstance.rho_nodes
@test length(rho_nodes) == columns[1].ConvDispOpInstance.nPoints
@test rho_nodes[1] ≈ rin atol=1e-10  # First node should be near inner radius
@test rho_nodes[end] ≈ rout atol=1e-10  # Last node should be near outer radius

println()
println("All tests passed!")
