using Test, Plots
include(joinpath(@__DIR__, "..", "..", "include.jl"))

# --- Problem setup ---
D_rad = 1e-5
t_end = 1
rhomin, rhomax = 0.01, 0.1

###
analytic_u(r, t) = sin(pi * (r - rhomin)/(rhomax - rhomin)) * exp(-D_rad * (pi/(rhomax - rhomin))^2 * t)
###

# Discretization settings
N = 3      # Polynomial order
Ne = 20    # Number of elements

vel_fn(r) = 0.0
u0_fn(r) = sin(pi * (r - rhomin)/(rhomax - rhomin))
kf = 0.0

m = Model(N, Ne, rhomin, rhomax, D_rad, 1.0, vel_fn, kf, u0_fn)
prob = ODEProblem((du,u,p,t)->f!(du,u,m,t), m.u0, (0.0, t_end))
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# --- Reshape solution vector for easy indexing: ---
u_mat = reshape(sol.u[end], N, Ne)

xi, _ = lglnodes(N)
deltarho = (rhomax - rhomin) / Ne

r_nodes_all = Float64[]
u_num_all = Float64[]
u_analytic_all = Float64[]

for e in 1:Ne
    rho_left = rhomin + (e-1)*deltarho
    rho_nodes = ((xi .+ 1)*(deltarho/2)) .+ rho_left
    for i in 1:N
        push!(r_nodes_all, rho_nodes[i])
        push!(u_num_all, u_mat[i, e])
        push!(u_analytic_all, analytic_u(rho_nodes[i], t_end))
    end
end

# --- Sort for smooth plotting ---
inds = sortperm(r_nodes_all)
r_sorted = r_nodes_all[inds]
u_num_sorted = u_num_all[inds]
u_analytic_sorted = u_analytic_all[inds]

# --- Plot DG vs Analytic ---
plot(
    r_sorted, u_num_sorted,
    marker=:circle, label="DG Solution",
    xlabel="r", ylabel="u",
    title="Analytical Solution Match (t = $t_end)", legend=:best
)
plot!(
    r_sorted, u_analytic_sorted,
    lw=2, ls=:dash, label="Analytical Solution"
)
savefig("analytical_match.png")

# --- (Optional) Print L2 error norm ---
l2err = sqrt(sum((u_num_sorted .- u_analytic_sorted).^2) / length(r_sorted))
println("L2 error norm: ", l2err)