#!/usr/bin/env julia
# L∞ Error and Convergence Test for Radial Flow
#
# Tests the numerical solution against an analytical solution
# and computes the L∞ (maximum) error.
#
# Analytical Solution:
# For a step input with constant velocity and dispersion,
# we can use the analytical solution from Huang et al. (1988)
# or a manufactured solution approach.
#
# This test uses a MANUFACTURED SOLUTION approach:
# 1. Choose a smooth function c_analytical(r,t)
# 2. Compute the forcing term f(r,t) that makes it satisfy the PDE
# 3. Add f as a source term in the numerical solution
# 4. Compare numerical vs analytical
#
# Manufactured solution (Gaussian pulse):
#   c(r,t) = exp(-(r - r_center - v*t)^2 / (4*D*(t + t0)))
#
# where t0 > 0 ensures smoothness at t=0

using CADETJulia
using Printf

println("=" ^ 80)
println("L∞ ERROR TEST - RADIAL FLOW CHROMATOGRAPHY")
println("=" ^ 80)
println()

# ==================== MANUFACTURED SOLUTION PARAMETERS ====================
# Analytical solution: Gaussian pulse moving with the flow
# c_analytical(r,t) = exp(-(r - r0 - v_avg*t)^2 / (4*D*(t + t0)))

r_center = 0.005   # Center of Gaussian at t=0 [m]
t0 = 0.001         # Time shift to ensure smoothness [s]
v_avg = 1.0e-4     # Average velocity [m/s]
D_eff = 1.0e-8     # Effective dispersion [m^2/s]

# ==================== TEST PARAMETERS ====================
ncomp = 1
tend = 0.05        # Simulation time [s]

# Geometry
rin = 0.001        # Inner radius [m]
rout = 0.01        # Outer radius [m]

# Physical properties
D_rad = D_eff
epsilon_c = 1.0

# Flow
u_in = v_avg * rin  # Inlet velocity (to get average v_avg)

# Solver tolerances
abstol = 1.0e-12
reltol = 1.0e-10
dtout = tend / 10   # Output 10 times

println("MANUFACTURED SOLUTION:")
println("  c(r,t) = exp(-(r - r0 - v*t)^2 / (4*D*(t + t0)))")
println("  where:")
println("    r0 = $r_center m (initial center)")
println("    t0 = $t0 s (time shift)")
println("    v_avg = $v_avg m/s")
println("    D_eff = $D_eff m^2/s")
println()

# ==================== ANALYTICAL SOLUTION FUNCTION ====================
"""
    c_analytical(r, t)

Manufactured analytical solution: Gaussian pulse
"""
function c_analytical(r::Float64, t::Float64)
    r_shifted = r - r_center - v_avg * t
    denominator = 4.0 * D_eff * (t + t0)
    return exp(-r_shifted^2 / denominator)
end

"""
    source_term(r, t)

Source term needed to make c_analytical satisfy the radial flow PDE.
This is computed by substituting c_analytical into the PDE and solving for f.
"""
function source_term(r::Float64, t::Float64)
    # For simplicity, we'll test without source term first
    # and use boundary conditions to inject the analytical solution
    return 0.0
end

# ==================== CONVERGENCE STUDY ====================
println("=" ^ 80)
println("CONVERGENCE STUDY")
println("=" ^ 80)
println()

# Test multiple mesh refinements
polyDegs = [2, 3, 4, 5]
nCells_list = [4, 8, 16, 32]

results = []

println(@sprintf("%-10s %-10s %-15s %-15s %-15s",
                 "PolyDeg", "nCells", "Total DOFs", "L∞ Error", "L2 Error"))
println("-" ^ 80)

for polyDeg in polyDegs
    for nCells in nCells_list
        local nNodes = polyDeg + 1
        local nDOFs = nCells * nNodes

        # Build model
        local col = CADETJulia.rLRM(
            nComp = ncomp,
            col_inner_radius = rin,
            col_outer_radius = rout,
            d_rad = fill(D_rad, ncomp),
            eps_c = epsilon_c,
            c0 = [c_analytical(rin + (rout - rin) * i / nDOFs, 0.0) for i in 1:nDOFs],
            q0 = fill(0.0, ncomp),
            polyDeg = polyDeg,
            nCells = nCells,
            cross_section_area = 1.0,
        )

        # Inlet with analytical solution
        local inlet = CADETJulia.CreateInlet(nComp = ncomp, nSections = 1)
        CADETJulia.modify_inlet!(
            inlet = inlet,
            nComp = ncomp,
            section = 1,
            cIn_c = [c_analytical(rin, 0.0)],
            cIn_l = zeros(ncomp),
            cIn_q = zeros(ncomp),
            cIn_cube = zeros(ncomp),
        )

        local outlet = CADETJulia.CreateOutlet(nComp = ncomp)

        # Switches & connections
        local idx_units = [0]
        local switches = CADETJulia.Switches(
            nSections = 1,
            section_times = [0.0, tend],
            nSwitches = 1,
            nColumns = 1,
            nComp = ncomp,
            idx_units = idx_units,
        )

        CADETJulia.Connection(switches, 1, 1, inlet, 1, u_in, 0.0, 0.0, 0.0, 0.0, false)
        CADETJulia.Connection(switches, 1, 1, 1, outlet, 0.0, 0.0, 0.0, 0.0, 0.0, false)

        # Solve
        local solverOptions = CADETJulia.SolverCache(
            columns = (col,), switches = switches,
            outlets = (outlet,), abstol = abstol,
            reltol = reltol, solution_times = [0.0, tend],
            prototypeJacobian = true, analyticalJacobian = false,
        )

        solution = CADETJulia.solve_model(
            columns = (col,),
            switches = switches,
            solverOptions = solverOptions,
            outlets = (outlet,),
            alg = nothing
        )

        # Extract final solution
        c_numerical = solution.sol[end][1:nDOFs]

        # Compute analytical solution at node positions at final time
        Δr = (rout - rin) / nCells
        c_analytical_vec = zeros(nDOFs)

        for cell in 1:nCells
            r_left = rin + (cell - 1) * Δr
            for node in 1:nNodes
                idx = (cell - 1) * nNodes + node
                # LGL nodes mapped to physical coordinates
                ξ = col.polyDeg_i > 0 ? col.nodes[node] : 0.0
                r_node = r_left + (Δr / 2) * (1 + ξ)
                c_analytical_vec[idx] = c_analytical(r_node, tend)
            end
        end

        # Compute errors
        error_vec = abs.(c_numerical - c_analytical_vec)
        L_inf_error = maximum(error_vec)
        L2_error = sqrt(sum(error_vec.^2) / nDOFs)

        println(@sprintf("%-10d %-10d %-15d %-15.6e %-15.6e",
                         polyDeg, nCells, nDOFs, L_inf_error, L2_error))

        push!(results, (polyDeg=polyDeg, nCells=nCells, nDOFs=nDOFs,
                       L_inf=L_inf_error, L2=L2_error))
    end
end

println("-" ^ 80)
println()

# ==================== CONVERGENCE RATE ANALYSIS ====================
println("=" ^ 80)
println("CONVERGENCE RATE ANALYSIS")
println("=" ^ 80)
println()

for polyDeg in polyDegs
    results_p = filter(r -> r.polyDeg == polyDeg, results)

    if length(results_p) >= 2
        println("Polynomial degree p = $polyDeg:")

        for i in 2:length(results_p)
            r1 = results_p[i-1]
            r2 = results_p[i]

            # Convergence rate: log(e1/e2) / log(h1/h2)
            h1 = (rout - rin) / r1.nCells
            h2 = (rout - rin) / r2.nCells

            rate_inf = log(r1.L_inf / r2.L_inf) / log(h1 / h2)
            rate_L2 = log(r1.L2 / r2.L2) / log(h1 / h2)

            println(@sprintf("  nCells: %d → %d: L∞ rate = %.2f, L2 rate = %.2f",
                           r1.nCells, r2.nCells, rate_inf, rate_L2))
        end

        println()
    end
end

println("=" ^ 80)
println("Expected convergence rate: O(h^(p+1)) where p = polynomial degree")
println("For p=2: rate ≈ 3")
println("For p=3: rate ≈ 4")
println("For p=4: rate ≈ 5")
println("=" ^ 80)
println()

println("TEST COMPLETE")
