#!/usr/bin/env julia
# EOC (Experimental Order of Convergence) Test for Radial Flow rLRMP Model
# Tests spatial convergence using a smooth Gaussian pulse injection
# with BOTH variable D(ρ) and variable k_f(ρ)

using Printf
using LinearAlgebra
using OrdinaryDiffEqBDF
using SciMLBase

# Add the CADETJulia module path
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "src"))

# Load CADETJulia module
using CADETJulia
using CADETJulia.DGElements
using CADETJulia.RadialConvDispOperatorDG
using CADETJulia: rLRMP, CreateInlet, modify_inlet!, CreateOutlet, Switches, Connection, compute_transport!

# ==================== PHYSICAL PARAMETERS ====================
rin = 0.1
rout = 1.1
v = 1.0e-2       # Interstitial velocity [m/s]
eps_c = 0.6      # Column porosity [-]
eps_p = 0.5      # Particle porosity [-]
Rp = 1.0e-4      # Particle radius [m]

# Variable dispersion coefficient D(ρ) - linear variation
D0 = 1.0e-4      # Base dispersion coefficient [m²/s]
D_slope = 5.0e-5 # Linear slope [m/s]

function D_variable(rho)
    return D0 + D_slope * (rho - rin)
end

# Variable film diffusion coefficient k_f(ρ) - linear variation
kf0 = 1.0e-3       # Base film diffusion coefficient [m/s]
kf_slope = 5.0e-4  # Linear slope [1/s]

function kf_variable(rho)
    return kf0 + kf_slope * (rho - rin)
end

# ==================== GAUSSIAN PULSE PARAMETERS ====================
t_center = 50.0      # Center time of Gaussian pulse [s]
t_width = 10.0       # Width (standard deviation) of Gaussian [s]
c_max = 1.0          # Maximum concentration

function inlet_gaussian(t)
    return c_max * exp(-((t - t_center) / t_width)^2)
end

# ==================== SIMULATION PARAMETERS ====================
t_final = 60.0
abstol = 1e-14
reltol = 1e-12

# ==================== EOC STUDY CONFIGURATION ====================
polyDeg_list = [1, 2, 3, 4, 5]
nCells_list = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]

# ==================== PRINT HEADER ====================
println("EOC Test - Gaussian Pulse (rLRMP) with Variable D(ρ) and k_f(ρ)")
println("D(ρ) = $D0 + $D_slope * (ρ - $rin)")
println("D(rin=$rin) = $(D_variable(rin)), D(rout=$rout) = $(D_variable(rout))")
println("k_f(ρ) = $kf0 + $kf_slope * (ρ - $rin)")
println("k_f(rin=$rin) = $(kf_variable(rin)), k_f(rout=$rout) = $(kf_variable(rout))")
println("Parameters: rin=$rin, rout=$rout, v=$v, eps_c=$eps_c, eps_p=$eps_p, Rp=$Rp")
println("t_final=$t_final")
println()

# Print headers for mobile phase
println("MOBILE PHASE (c):")
@printf("%-8s %-8s %-14s %-14s %-10s %-10s\n", "polyDeg", "nCells", "L2_err", "Linf_err", "EOC_L2", "EOC_Linf")
println("-" ^ 70)

# Store all results for printing
all_results_c = Dict{Int, Vector{Tuple{Int, Float64, Float64, String, String}}}()
all_results_cp = Dict{Int, Vector{Tuple{Int, Float64, Float64, String, String}}}()

# ==================== LOOP OVER POLYNOMIAL DEGREES ====================
for polyDeg in polyDeg_list
    n_levels = length(nCells_list)
    solutions_c = []   # Mobile phase
    solutions_cp = []  # Pore phase
    rho_nodes_list = []

    # Solve at each refinement level
    for (level, nCells) in enumerate(nCells_list)
        nNodes = polyDeg + 1
        nPoints = nCells * nNodes

        # Create rLRMP model with variable D and k_f
        model = rLRMP(
            nComp = 1,
            col_Rho_c = rin,
            col_Rho = rout,
            col_height = 1.0,
            d_rad = D_variable,
            eps_c = eps_c,
            eps_p = eps_p,
            kf = kf_variable,
            Rp = Rp,
            polyDeg = polyDeg,
            nCells = nCells,
            c0 = 0.0,
            cp0 = 0.0,
            q0 = 0.0
        )

        # Get radial node positions (now available directly from model)
        rho_nodes = model.ConvDispOpInstance.rho_nodes

        # Setup switches infrastructure
        inlet_obj = CreateInlet(nComp = 1, nSections = 1)
        modify_inlet!(inlet = inlet_obj, nComp = 1, section = 1,
                      cIn_c = [0.0], cIn_l = [0.0], cIn_q = [0.0], cIn_cube = [0.0])
        outlet_obj = CreateOutlet(nComp = 1)

        idx_units = [0]
        switches = Switches(nSections = 1, section_times = [0.0, t_final],
                           nSwitches = 1, nColumns = 1, nComp = 1, idx_units = idx_units)

        Connection(switches, 1, 1, inlet_obj, 1, v, 0.0, 0.0, 0.0, 0.0, false)
        Connection(switches, 1, 1, 1, outlet_obj, 0.0, 0.0, 0.0, 0.0, 0.0, false)

        # ODE right-hand side
        function ode_rhs!(dy, y, p, t)
            cIn = inlet_gaussian(t)
            model.cIn[1] = cIn

            RHS = zeros(Float64, model.unitStride)
            RHS_q = zeros(Float64, nPoints)

            x = zeros(Float64, model.unitStride)
            x[1:nPoints] .= y[1:nPoints]
            x[nPoints+1:2*nPoints] .= y[nPoints+1:2*nPoints]

            cpp = model.cpp
            compute_transport!(RHS, RHS_q, cpp, x, model, t, 1, 1, switches, idx_units)

            dy[1:nPoints] .= RHS[1:nPoints]
            dy[nPoints+1:2*nPoints] .= RHS[nPoints+1:2*nPoints]

            return nothing
        end

        # Solve ODE
        y_initial = zeros(2*nPoints)
        tspan = (0.0, t_final)
        prob = ODEProblem(ode_rhs!, y_initial, tspan)
        sol = solve(prob, QNDF(autodiff=false), abstol=abstol, reltol=reltol,
                    save_everystep=false, save_start=false)

        push!(solutions_c, sol.u[end][1:nPoints])
        push!(solutions_cp, sol.u[end][nPoints+1:2*nPoints])
        push!(rho_nodes_list, rho_nodes)
    end

    # Compute EOC using finest mesh as reference
    c_ref = solutions_c[end]
    cp_ref = solutions_cp[end]
    rho_ref = rho_nodes_list[end]
    nCells_ref = nCells_list[end]
    nNodes_ref = polyDeg + 1
    deltarho_ref = (rout - rin) / nCells_ref

    # High-order Lagrange interpolation
    function lagrange_interp(rho_target, c_ref, rho_ref_nodes, cell_start, nNodes)
        result = 0.0
        for j in 1:nNodes
            Lj = 1.0
            for k in 1:nNodes
                if k != j
                    Lj *= (rho_target - rho_ref_nodes[cell_start + k - 1]) /
                          (rho_ref_nodes[cell_start + j - 1] - rho_ref_nodes[cell_start + k - 1])
                end
            end
            result += c_ref[cell_start + j - 1] * Lj
        end
        return result
    end

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

        c_ref_interp = zeros(length(rho_coarse))
        cp_ref_interp = zeros(length(rho_coarse))

        for i in 1:length(rho_coarse)
            rho_target = rho_coarse[i]
            cell_ref = min(nCells_ref, max(1, Int(floor((rho_target - rin) / deltarho_ref)) + 1))
            if rho_target <= rin
                cell_ref = 1
            elseif rho_target >= rout
                cell_ref = nCells_ref
            end
            cell_start = (cell_ref - 1) * nNodes_ref + 1
            c_ref_interp[i] = lagrange_interp(rho_target, c_ref, rho_ref, cell_start, nNodes_ref)
            cp_ref_interp[i] = lagrange_interp(rho_target, cp_ref, rho_ref, cell_start, nNodes_ref)
        end

        diff_c = c_coarse .- c_ref_interp
        L2_error_c = sqrt(sum(diff_c.^2) / length(diff_c))
        Linf_error_c = maximum(abs.(diff_c))

        diff_cp = cp_coarse .- cp_ref_interp
        L2_error_cp = sqrt(sum(diff_cp.^2) / length(diff_cp))
        Linf_error_cp = maximum(abs.(diff_cp))

        push!(errors_L2_c, L2_error_c)
        push!(errors_Linf_c, Linf_error_c)
        push!(errors_L2_cp, L2_error_cp)
        push!(errors_Linf_cp, Linf_error_cp)
        push!(h_list, h)
    end

    # Store results for this polynomial degree
    results_c = Tuple{Int, Float64, Float64, String, String}[]
    results_cp = Tuple{Int, Float64, Float64, String, String}[]

    for i in 1:(n_levels-1)
        nCells = nCells_list[i]
        if i == 1
            push!(results_c, (nCells, errors_L2_c[i], errors_Linf_c[i], "--", "--"))
            push!(results_cp, (nCells, errors_L2_cp[i], errors_Linf_cp[i], "--", "--"))
        else
            eoc_L2_c = log(errors_L2_c[i-1] / errors_L2_c[i]) / log(h_list[i-1] / h_list[i])
            eoc_Linf_c = log(errors_Linf_c[i-1] / errors_Linf_c[i]) / log(h_list[i-1] / h_list[i])
            eoc_L2_cp = log(errors_L2_cp[i-1] / errors_L2_cp[i]) / log(h_list[i-1] / h_list[i])
            eoc_Linf_cp = log(errors_Linf_cp[i-1] / errors_Linf_cp[i]) / log(h_list[i-1] / h_list[i])
            push!(results_c, (nCells, errors_L2_c[i], errors_Linf_c[i], @sprintf("%.2f", eoc_L2_c), @sprintf("%.2f", eoc_Linf_c)))
            push!(results_cp, (nCells, errors_L2_cp[i], errors_Linf_cp[i], @sprintf("%.2f", eoc_L2_cp), @sprintf("%.2f", eoc_Linf_cp)))
        end
    end

    all_results_c[polyDeg] = results_c
    all_results_cp[polyDeg] = results_cp

    # Print mobile phase results
    for (nCells, L2_err, Linf_err, eoc_L2, eoc_Linf) in results_c
        @printf("%-8d %-8d %.6e  %.6e  %-10s %-10s\n", polyDeg, nCells, L2_err, Linf_err, eoc_L2, eoc_Linf)
    end
end

println("-" ^ 70)
println("Expected: EOC ~ p+1 for smooth solutions")
println()

# Print pore phase results
println("PORE PHASE (cp):")
@printf("%-8s %-8s %-14s %-14s %-10s %-10s\n", "polyDeg", "nCells", "L2_err", "Linf_err", "EOC_L2", "EOC_Linf")
println("-" ^ 70)

for polyDeg in polyDeg_list
    for (nCells, L2_err, Linf_err, eoc_L2, eoc_Linf) in all_results_cp[polyDeg]
        @printf("%-8d %-8d %.6e  %.6e  %-10s %-10s\n", polyDeg, nCells, L2_err, Linf_err, eoc_L2, eoc_Linf)
    end
end

println("-" ^ 70)
println("Expected: EOC ~ p+1 for smooth solutions with variable D(ρ) and k_f(ρ)")
