#!/usr/bin/env julia
# EOC (Experimental Order of Convergence) Test for Radial Flow using rLRMP Model
# Tests spatial convergence using a smooth Gaussian pulse injection
# with VARIABLE film diffusion coefficient k_f(ρ)

using Printf
using LinearAlgebra
using OrdinaryDiffEqBDF

# Add the CADETJulia module path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))

# Load CADETJulia module
using CADETJulia

# Also import the specific modules needed
using CADETJulia.DGElements
using CADETJulia.RadialConvDispOperatorDG

# Import the model constructor and support functions
using CADETJulia: rLRMP, CreateInlet, modify_inlet!, CreateOutlet, Switches, Connection, compute_transport!

println("=" ^ 80)
println("EOC TEST - GAUSSIAN PULSE INJECTION (rLRMP Model)")
println("Variable Film Diffusion Coefficient k_f(ρ)")
println("=" ^ 80)
println()

# ==================== PHYSICAL PARAMETERS ====================
rin = 0.1
rout = 1.1
v = 2/60       # Interstitial velocity [m/s]
D = 1.0e-4     # Dispersion coefficient [m²/s]
eps_c = 0.6    # Column porosity [-]
eps_p = 0.5    # Particle porosity [-]
Rp = 1.0e-4    # Particle radius [m]
col_height = 1.0  # Column height [m]

# Variable film diffusion coefficient k_f(ρ) - linear variation
kf_base = 1.0e-3       # Base film diffusion coefficient [m/s]
kf_slope = 5.0e-4      # Linear slope [1/s]

# k_f(ρ) = kf_base + kf_slope * (ρ - rin)
# This gives k_f(rin) = kf_base and k_f(rout) = kf_base + kf_slope * (rout - rin)
function kf_variable(rho)
    return kf_base + kf_slope * (rho - rin)
end

println("Physical Parameters:")
println("  Inner radius: $rin m")
println("  Outer radius: $rout m")
println("  Column height: $col_height m")
println("  Velocity: $v m/s")
println("  Dispersion: $D m²/s")
println("  Column porosity: $eps_c")
println("  Particle porosity: $eps_p")
println("  Particle radius: $Rp m")
println()
println("Variable Film Diffusion:")
println("  k_f(ρ) = $kf_base + $kf_slope * (ρ - $rin)")
println("  k_f(rin=$rin) = $(kf_variable(rin))")
println("  k_f(rout=$rout) = $(kf_variable(rout))")
println()

# ==================== GAUSSIAN PULSE PARAMETERS ====================
t_center = 10.0      # Center time of Gaussian pulse [s]
t_width = 3.0        # Width (standard deviation) of Gaussian [s]
c_max = 1.0          # Maximum concentration

function inlet_gaussian(t)
    """Smooth Gaussian pulse inlet boundary condition."""
    return c_max * exp(-((t - t_center) / t_width)^2)
end

println("Inlet Boundary Condition:")
println("  Gaussian pulse: c(t) = $c_max * exp(-((t - $t_center) / $t_width)²)")
println()

# ==================== SIMULATION PARAMETERS =================
t_final = 25.0
abstol = 1e-16
reltol = 1e-14

println("Simulation Parameters:")
println("  Final time: $t_final s")
println("  ODE tolerances: abstol=$abstol, reltol=$reltol")
println()

# ==================== EOC STUDY CONFIGURATION ====================
polyDeg_list = [1, 2, 3, 4, 5]
nCells_list = [1, 2, 4, 8, 16, 32, 64]

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
    solutions_c = []   # Mobile phase
    solutions_cp = []  # Pore phase
    rho_nodes_list = []
    nPoints_list = []
    compute_times = []

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

        println("  Total DOFs (mobile): $nPoints")
        println("  Total DOFs (pore): $nPoints")
        println("  Total DOFs: $(2*nPoints)")

        # Create rLRMP model with variable k_f
        model = rLRMP(
            nComp = 1,
            col_Rho_c = rin,
            col_Rho = rout,
            col_height = col_height,
            d_rad = D,
            eps_c = eps_c,
            eps_p = eps_p,
            kf = kf_variable,  # Variable film diffusion coefficient
            Rp = Rp,
            polyDeg = polyDeg,
            nCells = nCells,
            c0 = 0.0,
            cp0 = 0.0,
            q0 = 0.0
        )

        # Get radial node positions (now available directly from model)
        rho_nodes = model.ConvDispOpInstance.rho_nodes

        # Setup switches infrastructure for compute_transport!
        inlet_obj = CreateInlet(nComp = 1, nSections = 1)
        modify_inlet!(
            inlet = inlet_obj,
            nComp = 1,
            section = 1,
            cIn_c = [0.0],
            cIn_l = [0.0],
            cIn_q = [0.0],
            cIn_cube = [0.0],
        )

        outlet_obj = CreateOutlet(nComp = 1)

        idx_units = [0]
        switches = Switches(
            nSections = 1,
            section_times = [0.0, t_final],
            nSwitches = 1,
            nColumns = 1,
            nComp = 1,
            idx_units = idx_units,
        )

        Connection(switches, 1, 1, inlet_obj, 1, v, 0.0, 0.0, 0.0, 0.0, false)
        Connection(switches, 1, 1, 1, outlet_obj, 0.0, 0.0, 0.0, 0.0, 0.0, false)

        # ODE right-hand side using the model interface
        # For rLRMP: state = [c (mobile), cp (pore)]
        function ode_rhs!(dy, y, p, t)
            # Update inlet concentration (time-dependent Gaussian pulse)
            cIn = inlet_gaussian(t)
            model.cIn[1] = cIn

            # Allocate RHS with proper size
            RHS = zeros(Float64, model.unitStride)
            RHS_q = zeros(Float64, nPoints)

            # State vector for compute_transport
            # For rLRMP: first nPoints = mobile, next nPoints = pore
            x = zeros(Float64, model.unitStride)
            x[1:nPoints] .= y[1:nPoints]           # Mobile phase
            x[nPoints+1:2*nPoints] .= y[nPoints+1:2*nPoints]  # Pore phase

            # Call compute_transport!
            cpp = model.cpp
            compute_transport!(RHS, RHS_q, cpp, x, model, t, 1, 1, switches, idx_units)

            # Extract residuals
            dy[1:nPoints] .= RHS[1:nPoints]           # Mobile phase
            dy[nPoints+1:2*nPoints] .= RHS[nPoints+1:2*nPoints]  # Pore phase

            return nothing
        end

        # Initial condition: zero everywhere (mobile + pore)
        y_initial = zeros(2*nPoints)

        # Solve ODE
        println("  Solving...")
        tspan = (0.0, t_final)
        prob = ODEProblem(ode_rhs!, y_initial, tspan)

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

        # Store solutions (mobile and pore phases separately)
        push!(solutions_c, sol.u[end][1:nPoints])
        push!(solutions_cp, sol.u[end][nPoints+1:2*nPoints])
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
    c_ref = solutions_c[end]
    cp_ref = solutions_cp[end]
    rho_ref = rho_nodes_list[end]

    # Compute errors relative to finest solution
    errors_L2_c = Float64[]
    errors_Linf_c = Float64[]
    errors_L2_cp = Float64[]
    errors_Linf_cp = Float64[]
    h_list = Float64[]

    for level in 1:(n_levels-1)
        nCells = nCells_list[level]
        h = (rout - rin) / nCells

        c_coarse = solutions_c[level]
        cp_coarse = solutions_cp[level]
        rho_coarse = rho_nodes_list[level]

        # Interpolate reference solution to coarse grid points
        c_ref_interp = zeros(length(rho_coarse))
        cp_ref_interp = zeros(length(rho_coarse))

        for i in 1:length(rho_coarse)
            rho_target = rho_coarse[i]
            idx_right = findfirst(r -> r >= rho_target, rho_ref)

            if idx_right === nothing
                c_ref_interp[i] = c_ref[end]
                cp_ref_interp[i] = cp_ref[end]
            elseif idx_right == 1
                c_ref_interp[i] = c_ref[1]
                cp_ref_interp[i] = cp_ref[1]
            else
                idx_left = idx_right - 1
                alpha = (rho_target - rho_ref[idx_left]) / (rho_ref[idx_right] - rho_ref[idx_left])
                c_ref_interp[i] = (1 - alpha) * c_ref[idx_left] + alpha * c_ref[idx_right]
                cp_ref_interp[i] = (1 - alpha) * cp_ref[idx_left] + alpha * cp_ref[idx_right]
            end
        end

        # Compute errors for mobile phase
        diff_c = c_coarse .- c_ref_interp
        L2_error_c = sqrt(sum(diff_c.^2) / length(diff_c))
        Linf_error_c = maximum(abs.(diff_c))

        # Compute errors for pore phase
        diff_cp = cp_coarse .- cp_ref_interp
        L2_error_cp = sqrt(sum(diff_cp.^2) / length(diff_cp))
        Linf_error_cp = maximum(abs.(diff_cp))

        push!(errors_L2_c, L2_error_c)
        push!(errors_Linf_c, Linf_error_c)
        push!(errors_L2_cp, L2_error_cp)
        push!(errors_Linf_cp, Linf_error_cp)
        push!(h_list, h)

        println("Level $level (nCells = $nCells, h = $(round(h, digits=6))):")
        @printf("  Mobile phase - L2: %.6e, L∞: %.6e\n", L2_error_c, Linf_error_c)
        @printf("  Pore phase   - L2: %.6e, L∞: %.6e\n", L2_error_cp, Linf_error_cp)
        println()
    end

    # ==================== COMPUTE CONVERGENCE RATES ====================
    println("-" ^ 80)
    println("CONVERGENCE RATES - MOBILE PHASE")
    println("-" ^ 80)
    println()

    @printf("%-8s %-12s %-14s %-14s %-10s %-10s\n",
            "Level", "nCells", "L2 error", "L∞ error", "EOC(L2)", "EOC(L∞)")
    println("-" ^ 80)

    for i in 1:(n_levels-1)
        nCells = nCells_list[i]

        if i == 1
            @printf("%-8d %-12d %.6e  %.6e  %-10s %-10s\n",
                    i, nCells, errors_L2_c[i], errors_Linf_c[i], "--", "--")
        else
            eoc_L2 = log(errors_L2_c[i-1] / errors_L2_c[i]) / log(h_list[i-1] / h_list[i])
            eoc_Linf = log(errors_Linf_c[i-1] / errors_Linf_c[i]) / log(h_list[i-1] / h_list[i])

            @printf("%-8d %-12d %.6e  %.6e  %-10.2f %-10.2f\n",
                    i, nCells, errors_L2_c[i], errors_Linf_c[i], eoc_L2, eoc_Linf)
        end
    end

    println()
    println("-" ^ 80)
    println("CONVERGENCE RATES - PORE PHASE")
    println("-" ^ 80)
    println()

    @printf("%-8s %-12s %-14s %-14s %-10s %-10s\n",
            "Level", "nCells", "L2 error", "L∞ error", "EOC(L2)", "EOC(L∞)")
    println("-" ^ 80)

    for i in 1:(n_levels-1)
        nCells = nCells_list[i]

        if i == 1
            @printf("%-8d %-12d %.6e  %.6e  %-10s %-10s\n",
                    i, nCells, errors_L2_cp[i], errors_Linf_cp[i], "--", "--")
        else
            eoc_L2 = log(errors_L2_cp[i-1] / errors_L2_cp[i]) / log(h_list[i-1] / h_list[i])
            eoc_Linf = log(errors_Linf_cp[i-1] / errors_Linf_cp[i]) / log(h_list[i-1] / h_list[i])

            @printf("%-8d %-12d %.6e  %.6e  %-10.2f %-10.2f\n",
                    i, nCells, errors_L2_cp[i], errors_Linf_cp[i], eoc_L2, eoc_Linf)
        end
    end

    println()
    println("Expected convergence rate: p+1 = $(polyDeg+1) for smooth solutions")
    println()

end  # End of polyDeg loop

println()
println("=" ^ 80)
println("EOC TEST COMPLETE - ALL POLYNOMIAL DEGREES (rLRMP Model with Variable k_f)")
println("=" ^ 80)
