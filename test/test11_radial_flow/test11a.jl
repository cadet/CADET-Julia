using Test
include(joinpath(@__DIR__, "..", "..", "include.jl"))

@testset "DG Diffusion: Analytical Solution Match" begin
    # Parameters
    N = 3
    Ne = 40
    rhomin, rhomax = 0.01, 0.1
    D_rad = 1e-8
    t_end = 0.1

    # Analytic solution function
    analytic_u(r, t) = sin(pi * (r - rhomin) / (rhomax - rhomin)) *
        exp(-D_rad * (pi / (rhomax - rhomin))^2 * t)

    # Problem setup
    vel_fn(r) = 0.0
    u0_fn(r) = sin(pi * (r - rhomin) / (rhomax - rhomin))
    kf = 0.0

    m = Model(N, Ne, rhomin, rhomax, D_rad, 1.0, vel_fn, kf, u0_fn)
    prob = ODEProblem((du, u, p, t) -> f!(du, u, m, t), m.u0, (0.0, t_end))
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    # Reconstruct solution at all DG nodes
    xi, _ = lglnodes(N)
    deltarho = (rhomax - rhomin) / Ne
    u_mat = reshape(sol.u[end], N, Ne)

    err_max = 0.0
    for e in 1:Ne
        r_left = rhomin + (e - 1) * deltarho
        r_nodes = ((xi .+ 1) * (deltarho / 2)) .+ r_left
        for i in 1:N
            u_num = u_mat[i, e]
            u_ex = analytic_u(r_nodes[i], t_end)
            err_max = max(err_max, abs(u_num - u_ex))
        end
    end

    @test err_max < 1e-4

    @info "Maximum error vs analytic: $err_max"
end