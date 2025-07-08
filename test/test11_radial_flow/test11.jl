#-------------------------------------------------------------------
# Test 11: Radial Discontinuous Galerkin (DG) + LRM mass-conservation
# This test verifies that the DG+LRM discretization for radial flow
# conserves total mass over time by comparing initial and final integrals.
#-------------------------------------------------------------------
using Test, Plots
include(joinpath(@__DIR__,"..","..","include.jl"))

"""
N, Ne = 3, 20
rhomin, rhomax = 0.01, 0.1
D_rad, tau = 1e-5, 1.0
t_end = 10.0
vel_fn(r) = 1.0
u0_fn(r) = 1.0      # non-zero initial condition
kf = 0.0    # no linear sink

# --- Binding model (no binding for pure mass conservation test) ---
# If  want to test with linear binding, uncomment below and update your Model and f! accordingly
# using .binding_base: Linear
# ka = [0.0]; kd = [0.0]; is_kinetic = true; kkin = [1.0]; nBound = [true]
# bindStride = N * Ne
# binding = Linear(; ka=ka, kd=kd, is_kinetic=is_kinetic, kkin=kkin, bindStride=bindStride, nBound=nBound)
# nComp = 1
# m = Model(N, Ne, rhomin, rhomax, D_rad, tau, vel_fn, kf, u0_fn, binding, nComp, bindStride)

# --- Standard no-binding model ---

m = Model(N, Ne, rhomin, rhomax, D_rad, tau, vel_fn, kf, u0_fn)
prob = ODEProblem((du,u,p,t)->f!(du,u,m,t), m.u0, (0.0, t_end))
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

xi, w = lglnodes(N)
deltarho = (rhomax - rhomin)/Ne

function dg_mass(u, N, Ne, rhomin, deltarho, w, xi)
    mass = 0.0; offset = 0
    for e in 1:Ne
        rho_left = rhomin + (e-1)*deltarho
        rho_nodes = ((xi .+ 1)*(deltarho/2)) .+ rho_left
        u_nodes = u[offset+1 : offset+N]
        mass += sum(u_nodes .* rho_nodes .* w) * (deltarho/2)
        offset += N
    end
    return mass
end
mass0 = dg_mass(m.u0, N, Ne, rhomin, deltarho, w, xi)
mass1 = dg_mass(sol.u[end], N, Ne, rhomin, deltarho, w, xi)

@test isapprox(mass1, mass0, rtol=1e-8)
"""
# Convergence Test
# Problem setup
D_rad = 1e-3
t_end = 0.1
rhomin, rhomax = 0.01, 0.1

analytic_u(r, t) = sin(pi * (r - rhomin)/(rhomax - rhomin)) * exp(-D_rad * (pi/(rhomax - rhomin))^2 * t)

# Choose polynomial order for DG (e.g. 3)
N = 3
err_list = Float64[]
h_list = Float64[]

for Ne in [5, 10, 20, 40, 80]
    deltarho = (rhomax - rhomin)/Ne
    vel_fn(r) = 0.0  # Pure diffusion
    u0_fn(r) = sin(pi * (r - rhomin)/(rhomax - rhomin))

    m = Model(N, Ne, rhomin, rhomax, D_rad, 1.0, vel_fn, 0.0, u0_fn)
    prob = ODEProblem((du,u,p,t)->f!(du,u,m,t), m.u0, (0.0, t_end))
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    xi, w = lglnodes(N)
    # Compute L2 error norm
    offset = 0
    err2 = 0.0
    for e in 1:Ne
        r_left = rhomin + (e-1)*deltarho
        r_nodes = ((xi .+ 1)*(deltarho/2)) .+ r_left
        u_num = sol.u[end][offset+1:offset+N]
        u_ex = analytic_u.(r_nodes, t_end)
        # L2 error on this element, weighted by quadrature and r
        err2 += sum(((u_num - u_ex).^2) .* r_nodes .* w) * (deltarho/2)
        offset += N
    end
    push!(err_list, sqrt(err2))
    push!(h_list, deltarho)
end

# Plot error vs. h
plot(h_list, err_list; xaxis=:log, yaxis=:log, marker=:o, xlabel="Element size h", ylabel="L2 error", label="DG error")
plot!(h_list, h_list .^ (N), l=:dash, label="O(h^$N)")

# Print observed order
rates = log.(err_list[1:end-1]) .- log.(err_list[2:end]) ./ (log.(h_list[1:end-1]) .- log.(h_list[2:end]))
println("Observed order: ", rates)