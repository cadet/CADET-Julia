#!/usr/bin/env julia
# Pulse Injection Test - Direct ODE Solver Approach
# This version uses direct ODE solver access for full solution data

using CADETJulia
using Printf
using Statistics

# ==================== TEST PARAMETERS ====================
ncomp = 1
tend = 130.0

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

# Solver tolerances
abstol = 1.0e-10
reltol = 1.0e-8
dtout = 0.1

println("=" ^ 80)
println("PULSE INJECTION TEST (Direct ODE Solver)")
println("=" ^ 80)
println()

# ==================== BUILD MODEL ====================

using CADETJulia.DGElements
nodes_ref, _ = DGElements.cglnodes(polyDeg)

deltarho = (rout - rin) / nCells
nNodes = polyDeg + 1
nPoints = nCells * nNodes

# Compute nodal positions
rho_nodes = zeros(nPoints)
for cell in 1:nCells
    for node in 1:nNodes
        idx = (cell - 1) * nNodes + node
        rho_nodes[idx] = rin + (cell - 1) * deltarho + (deltarho / 2) * (1 + nodes_ref[node])
    end
end

col = CADETJulia.rLRM(
    nComp = ncomp,
    col_inner_radius = rin,
    col_outer_radius = rout,
    d_rad = fill(D_rad, ncomp),
    eps_c = epsilon_c,
    c0 = fill(0.0, ncomp),  # Clean column
    q0 = fill(0.0, ncomp),
    polyDeg = polyDeg,
    nCells = nCells,
)

u_in = Q / (col.cross_section_area * col.eps_c)

println("Model setup:")
println("  u_in = $u_in m/s")
println("  nPoints = $nPoints")
println()

# Inlet: Three sections
inlet = CADETJulia.CreateInlet(nComp = ncomp, nSections = 3)

CADETJulia.modify_inlet!(inlet = inlet, nComp = ncomp, section = 1,
    cIn_c = fill(0.0, ncomp), cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp), cIn_cube = zeros(ncomp))

CADETJulia.modify_inlet!(inlet = inlet, nComp = ncomp, section = 2,
    cIn_c = fill(1.0, ncomp), cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp), cIn_cube = zeros(ncomp))

CADETJulia.modify_inlet!(inlet = inlet, nComp = ncomp, section = 3,
    cIn_c = fill(0.0, ncomp), cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp), cIn_cube = zeros(ncomp))

outlet = CADETJulia.CreateOutlet(nComp = ncomp)

# ==================== RUN SIMULATION SECTION BY SECTION ====================

solution_times_all = Float64[]
c_outlet_all = Float64[]

section_times = [0.0, 1.0, 60.0, 130.0]
section_c_in = [0.0, 1.0, 0.0]  # Inlet concentration for each section

# Initialize state
x0 = zeros(Float64, col.unitStride)

for isection in 1:3
    t_start = section_times[isection]
    t_end = section_times[isection + 1]
    c_in = section_c_in[isection]

    println("Running section $isection: t = $t_start to $t_end s, c_in = $c_in")

    # Create a simple single-section inlet for this section
    section_inlet = CADETJulia.CreateInlet(nComp = ncomp, nSections = 1)
    CADETJulia.modify_inlet!(inlet = section_inlet, nComp = ncomp, section = 1, cIn_c = fill(c_in, ncomp), cIn_l = zeros(ncomp), cIn_q = zeros(ncomp), cIn_cube = zeros(ncomp))

    # Setup switches for this section
    idx_units = [0]
    switches = CADETJulia.Switches(nSections = 1, section_times = [t_start, t_end], nSwitches = 1, nColumns = 1, nComp = ncomp, idx_units = idx_units,)

    CADETJulia.Connection(switches, 1, 1, section_inlet, 1, u_in, 0.0, 0.0, 0.0, 0.0, false)
    CADETJulia.Connection(switches, 1, 1, 1, outlet, 0.0, 0.0, 0.0, 0.0, 0.0, false)

    section_solution_times = collect(t_start:dtout:t_end)

    # Direct ODE solve
    RHS_q = col.RHS_q
    cpp = col.cpp
    qq = col.qq
    i = 1
    p = ((col,), RHS_q, cpp, qq, i, 1, idx_units, switches, nothing)

    # Set inlet concentration
    col.cIn[1] = c_in

    fun = SciMLBase.ODEFunction(CADETJulia.problem!)
    prob = SciMLBase.ODEProblem(fun, x0, (t_start, t_end), p)
    sol = SciMLBase.solve(prob, QNDF(autodiff=AutoFiniteDiff()), saveat=section_solution_times,
                          abstol=abstol, reltol=reltol)

    println("  Section complete, $(length(sol.t)) time points")

    # Extract outlet concentration (last spatial point)
    start_idx = (isection == 1) ? 1 : 2  # Include first point only for first section
    for j in start_idx:length(sol.t)
        push!(solution_times_all, sol.t[j])
        # Outlet concentration is at the last spatial node (nPoints)
        c_out = sol.u[j][nPoints]
        push!(c_outlet_all, c_out)
    end

    # Update initial condition for next section
    x0 .= sol.u[end]
end

println()
println("Simulation complete!")
println("  Total time points: $(length(solution_times_all))")
println()

# ==================== SAVE RESULTS ====================

output_dir = @__DIR__
csv_file = joinpath(output_dir, "test12_pulse_direct_results.csv")

open(csv_file, "w") do io
    println(io, "time,c_outlet,c_inlet")

    for i in 1:length(solution_times_all)
        t = solution_times_all[i]
        c_out = c_outlet_all[i]

        # Determine inlet
        if t < 1.0
            c_in = 0.0
        elseif t < 60.0
            c_in = 1.0
        else
            c_in = 0.0
        end

        @printf(io, "%.3f,%.6e,%.6e\n", t, c_out, c_in)
    end
end

println("Saved: $csv_file")
println()

# ==================== ANALYSIS ====================

println("=" ^ 80)
println("RESULTS ANALYSIS")
println("=" ^ 80)
println()

# Find peak
if length(c_outlet_all) > 0
    idx_max = argmax(c_outlet_all)
    c_max = c_outlet_all[idx_max]
    t_max = solution_times_all[idx_max]

    @printf("Peak concentration: %.6f at t = %.2f s\n", c_max, t_max)
    @printf("Initial concentration: %.6f\n", c_outlet_all[1])
    @printf("Final concentration: %.6f\n", c_outlet_all[end])
    println()

    # Sample output
    println("Sample data points:")
    for i in [1, div(length(c_outlet_all), 4), div(length(c_outlet_all), 2),
              div(3*length(c_outlet_all), 4), length(c_outlet_all)]
        @printf("  t = %.1f s: c = %.6f\n", solution_times_all[i], c_outlet_all[i])
    end
end

println()
println("=" ^ 80)
println("Use plot_pulse_results.py to visualize (after renaming CSV file)")
println("=" ^ 80)
