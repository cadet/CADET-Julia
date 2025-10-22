

using CADETJulia
using Printf

# ------------------ User knobs ------------------
ncomp  = 1
c0     = 1.0          # inlet concentration during pulse
tinj   = 60.0         # pulse duration [s]
tend   = 130.0        # final time [s]
rin, rout = 0.025, 0.5
D_rad  = 1.0e-8
εc     = 0.40
polyDeg, nCells = 4, 40
u_in   = 1.0e-4       # superficial velocity at inner radius [m/s]

abstol = 1.0e-10
reltol = 1.0e-8
Δtout  = 0.5

# ------------------ Build units manually ------------------
# Column (radial LRM)
col = CADETJulia.rLRM(nComp = ncomp, col_inner_radius = rin, col_outer_radius = rout, d_rad = fill(D_rad, ncomp), eps_c = εc, c0 = fill(0.0, ncomp), q0 = fill(0.0, ncomp), polyDeg = polyDeg, nCells = nCells, cross_section_area = 1.0,)

# Inlet with 2 sections (pulse then wash)
inlet = CADETJulia.CreateInlet(nComp = ncomp, nSections = 2)
CADETJulia.modify_inlet!(
    inlet = inlet,
    nComp = ncomp,
    section = 1,                 # section 1: t ∈ [0, tinj)
    cIn_c = fill(c0, ncomp),
    cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp),
    cIn_cube = zeros(ncomp),
)
CADETJulia.modify_inlet!(
    inlet = inlet,
    nComp = ncomp,
    section = 2,                 # section 2: t ∈ [tinj, tend]
    cIn_c = fill(0.0, ncomp),
    cIn_l = zeros(ncomp),
    cIn_q = zeros(ncomp),
    cIn_cube = zeros(ncomp),
)

# Outlet
outlet = CADETJulia.CreateOutlet(nComp = ncomp)

# ------------------ Switching & connections ------------------
# idx_units offsets: for a single column, it starts at 0
idx_units = [0]

switches = CADETJulia.Switches(
    nSections = 2,
    section_times = [0.0, tinj, tend],
    nSwitches = 1,
    nColumns = 1,
    nComp = ncomp,
    idx_units = idx_units,
)

# Register connections
# 1) INLET (unit_000) -> COLUMN (column index 1). For columns, the 5th arg is velocity u
CADETJulia.Connection(
    switches,           # switches object to mutate
    1,                  # switch index (1-based)
    1,                  # section index (1: applies to both unless you add more switches)
    inlet,              # source (inlet instance)
    1,                  # sink (column number in [1..nColumns])
    u_in,               # superficial velocity at inner face; radial profile handled internally
    0.0, 0.0, 0.0, 0.0, # dynamic flow coeffs (unused)
    false,              # dynamic_flow_check
)

# 2) COLUMN -> OUTLET. For outlets, the value is ignored
CADETJulia.Connection(
    switches,
    1,
    1,
    1,                  # source (column number)
    outlet,             # sink (outlet instance)
    0.0,
    0.0, 0.0, 0.0, 0.0,
    false,
)

# ------------------ Solver options ------------------
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

# ------------------ Solve ------------------
sol = CADETJulia.solve_model(columns=(col,), switches=switches, outlets=(outlet,), solverOptions=solverOptions)

# ------------------ Report outlet ------------------
# The column struct stores outlet trace & times; print quick summary
@printf("\nPulse injection run finished.\n")
@printf("  tinj = %.3f s, u_in = %.3e m/s, D_rad = %.3e m^2/s\n", tinj, u_in, D_rad)



# If available, dump last few samples
if !isempty(col.solution_times)
    nt = length(col.solution_times)
    nshow = min(nt, 400)
    @printf("\nLast %d samples at outlet (time, c_out):\n", nshow)
    for k in (nt - nshow + 1):nt
        t = col.solution_times[k]
        c = col.solution_outlet[k, 1]  # first component
        @printf("  %10.3f  %14.6e\n", t, c)
    end
end


#############################################
# Plot outlet concentration for [0,60,130] #
#############################################
try
    @eval begin
        using Plots
    end
    tmin, tmid, tmax = 0.0, 60.0, 130.0
    mask = (col.solution_times .>= tmin) .& (col.solution_times .<= tmax)
    plt = plot(
        col.solution_times[mask],
        col.solution_outlet[mask, 1],
        xlabel = "Time [s]",
        ylabel = "Outlet concentration",
        title = "Radial column pulse elution (sections: [0,60,130])",
        legend = false,
        lw = 2,
    )
    vline!([tmin, tmid, tmax], l = (:dash, 1))
    display(plt)
catch err
    @warn "Plotting failed (Plots.jl not available?)." exception=(err, catch_backtrace())
end