#-------------------------------------------------------------------
# Test 11: Radial Discontinuous Galerkin (DG) + LRM mass-conservation
# This test verifies that the DG+LRM discretization for radial flow
# conserves total mass over time by comparing initial and final integrals.
#-------------------------------------------------------------------
using Test
include(joinpath(@__DIR__,"..","..","include.jl"))

N, Ne = 3, 20
rhomin, rhomax = 0.01, 0.1
D_a, tau = 1e-5, 1.0
t_end = 10.0
vel_fn(r) = 1.0
u0_fn(r) = 1.0   # <-- nonzero initial condition
kf = 0.0

m = Model(N, Ne, rhomin, rhomax, D_a, tau, vel_fn, kf, u0_fn)
prob = ODEProblem((du,u,p,t)->f!(du,u,m,t), m.u0, (0.0, t_end))
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

xi, w = lglnodes(N)
deltarho = (rhomax - rhomin)/Ne

function dg_mass(u, N, Ne, rhomin, deltarho, w)
    mass = 0.0; offset = 0
    for e in 1:Ne
        r_left = rhomin + (e-1)*deltarho
        r_nodes = ((xi .+ 1)*(deltarho/2)) .+ r_left
        u_nodes = u[offset+1 : offset+N]
        mass += sum(u_nodes .* r_nodes .* w) * (deltarho/2)
        offset += N
    end
    return mass
end

mass0 = dg_mass(m.u0, N, Ne, rhomin, deltarho, w)
mass1 = dg_mass(sol.u[end], N, Ne, rhomin, deltarho, w)

@test isapprox(mass1, mass0, rtol=1e-8)