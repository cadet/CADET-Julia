#!/usr/bin/env julia
# Extract concentration vs time data from test12d (free stream) for plotting
# Outputs CSV files that can be plotted with external tools

using CADETJulia
using Printf
using Statistics

# ==================== TEST PARAMETERS ====================
ncomp = 1
tend = 0.5

# Geometry
rin = 0.01
rout = 0.1

# Physical properties
D_rad = 1.0e-4
epsilon_c = 1.0 

# Discretization
polyDeg = 4
nCells = 100

# Flow
u_in = 1.0e-2

# Solver tolerances
abstol = 1.0e-10
reltol = 1.0e-8
dtout = 0.001

println("Running free stream simulation and extracting data for plotting...")

# ==================== BUILD MODEL ====================

using CADETJulia.DGElements
nodes_ref, _ = DGElements.cglnodes(polyDeg)

deltarho = (rout - rin) / nCells
nNodes = polyDeg + 1
nPoints = nCells * nNodes

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
    c0 = rho_nodes,
    q0 = fill(0.0, ncomp),
    polyDeg = polyDeg,
    nCells = nCells,
)

inlet = CADETJulia.CreateInlet(nComp = ncomp, nSections = 1)
CADETJulia.modify_inlet!(
    inlet = inlet,
    nComp = ncomp,
    section = 1,
    cIn_c = fill(rin, ncomp),  # Match initial condition at inlet: c = ρ_inlet
    cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp),
    cIn_cube = zeros(ncomp),
)

outlet = CADETJulia.CreateOutlet(nComp = ncomp)

idx_units = [0]
switches = CADETJulia.Switches(
    nSections = 1,
    section_times = [0.0, tend],
    nSwitches = 1,
    nColumns = 1,
    nComp = ncomp,
    idx_units = idx_units,
)

CADETJulia.Connection(switches, 1, 1, inlet, 1, u_in, 0.0, 0.0, 0.0, 0.0, false)
CADETJulia.Connection(switches, 1, 1, 1, outlet, 0.0, 0.0, 0.0, 0.0, 0.0, false)

solution_times = collect(0.0:dtout:tend)
solverOptions = CADETJulia.SolverCache(
    columns = (col,), switches = switches,
    outlets = (outlet,), abstol = abstol,
    reltol = reltol, solution_times = solution_times,
    prototypeJacobian = true, analyticalJacobian = false,
)

# ==================== RUN SIMULATION ====================

RHS_q = col.RHS_q
cpp = col.cpp
qq = col.qq
i = 1
p = (columns=(col,), RHS_q, cpp, qq, i, solverOptions.nColumns, solverOptions.idx_units, switches, nothing)

for j = 1:ncomp
    col.cIn[j] = switches.ConnectionInstance.cIn_c[1][1][j][1]
end

fun = SciMLBase.ODEFunction(CADETJulia.problem!)
prob = SciMLBase.ODEProblem(fun, solverOptions.x0, (0.0, tend), p)
sol = SciMLBase.solve(prob, QNDF(autodiff=AutoFiniteDiff()), saveat=solution_times, abstol=abstol, reltol=reltol)

println("Simulation complete! Extracting data...")

# ==================== EXTRACT DATA ====================

times = sol.t
n_times = length(times)

# Select specific spatial locations
idx_inlet = 1
idx_quarter = div(nPoints, 4)
idx_mid = div(nPoints, 2)
idx_3quarter = div(3*nPoints, 4)
idx_outlet = nPoints

# ==================== WRITE CSV FILES ====================

# File 1: Concentration vs Time at specific locations (chromatogram data)
output_dir = @__DIR__
csv_file_1 = joinpath(output_dir, "test12a_c_vs_t.csv")

open(csv_file_1, "w") do io
    # Header
    println(io, "time,c_inlet,c_quarter,c_mid,c_3quarter,c_outlet")
    println(io, "# rho_inlet=$(rho_nodes[idx_inlet]), rho_quarter=$(rho_nodes[idx_quarter]), rho_mid=$(rho_nodes[idx_mid]), rho_3quarter=$(rho_nodes[idx_3quarter]), rho_outlet=$(rho_nodes[idx_outlet])")
    println(io, "# Free stream test: D=$D_rad, u_in=$u_in m^3/s")

    # Data
    for i in 1:n_times
        @printf(io, "%.3f,%.6e,%.6e,%.6e,%.6e,%.6e\n",
            times[i],
            sol.u[i][idx_inlet],
            sol.u[i][idx_quarter],
            sol.u[i][idx_mid],
            sol.u[i][idx_3quarter],
            sol.u[i][idx_outlet])
    end
end

println("Saved: $csv_file_1")
# ==================== SIMPLE TEXT-BASED VISUALIZATION ====================

println("\n" * "="^80)
println("SIMPLE TEXT VISUALIZATION")
println("="^80)
println()

# Plot concentration at inlet, mid, and outlet vs time
println("Concentration Evolution (Text Plot)")
println("-"^80)
println("Legend: I=Inlet(ρ=$(round(rho_nodes[idx_inlet],digits=4))), M=Mid(ρ=$(round(rho_nodes[idx_mid],digits=4))), O=Outlet(ρ=$(round(rho_nodes[idx_outlet],digits=4)))")
println()

# Create a simple ASCII plot
plot_width = 70
plot_height = 20

c_inlet_data = [sol.u[i][idx_inlet] for i in 1:n_times]
c_mid_data = [sol.u[i][idx_mid] for i in 1:n_times]
c_outlet_data = [sol.u[i][idx_outlet] for i in 1:n_times]

all_c = vcat(c_inlet_data, c_mid_data, c_outlet_data)
c_min_all = minimum(all_c)
c_max_all = maximum(all_c)

println(@sprintf("Concentration range: [%.6f, %.6f]", c_min_all, c_max_all))
println()
println(@sprintf("Time [s] | %s", "Concentration"))
println("-"^80)

for i in 1:n_times
    t = times[i]
    c_i = c_inlet_data[i]
    c_m = c_mid_data[i]
    c_o = c_outlet_data[i]

    # Normalize to plot width
    pos_i = round(Int, (c_i - c_min_all) / (c_max_all - c_min_all) * plot_width)
    pos_m = round(Int, (c_m - c_min_all) / (c_max_all - c_min_all) * plot_width)
    pos_o = round(Int, (c_o - c_min_all) / (c_max_all - c_min_all) * plot_width)

    line = fill(' ', plot_width + 1)
    if pos_i >= 1 && pos_i <= plot_width
        line[pos_i] = 'I'
    end
    if pos_m >= 1 && pos_m <= plot_width
        line[pos_m] = 'M'
    end
    if pos_o >= 1 && pos_o <= plot_width
        line[pos_o] = 'O'
    end

    @printf("%7.2f  | %s\n", t, String(line))
end

println("-"^80)
println()

println("Summary:")
println("  - Inlet concentration (I) changes from $(round(c_inlet_data[1], digits=6)) to $(round(c_inlet_data[end], digits=6))")
println("  - Mid concentration (M) changes from $(round(c_mid_data[1], digits=6)) to $(round(c_mid_data[end], digits=6))")
println("  - Outlet concentration (O) changes from $(round(c_outlet_data[1], digits=6)) to $(round(c_outlet_data[end], digits=6))")
println()
println("CSV files contain full data for plotting with Python, MATLAB, Excel, etc.")
println("="^80)