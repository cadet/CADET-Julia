#-------------------------------------------------------------------
# Test 11: Radial Discontinuous Galerkin (DG) + LRM mass-conservation
# This test verifies that the DG+LRM discretization for radial flow
# conserves total mass over time by comparing initial and final integrals.
#-------------------------------------------------------------------
using Test
include(joinpath(@__DIR__,"..","..","include.jl"))

const N      =  3    # polynomial order: number of basis functions per element
const Ne     = 20    # number of radial elements discretizing the domain
const rhomin = 0.01  # inner radius of the annular domain
const rhomax = 0.1   # outer radius of the annular domain
const D_a    = 1e-5  # disperho_sion coefficient in radial coordinate
const tau    = 1.0   # penalty parameter for DG interface flux stabilization
const t_end  = 10.0   # final time for the simulation
vel_fn(r)   = 1.0    # velocity profile (constant radial velocity)
u0_fn(r)    = 0.0    # initial concentration profile as a function of r

@testset "Radial DG + LRM mass-conservation" begin
  # parameter
  m = Model(N, Ne, rhomin, rhomax, D_a, tau, vel_fn, t_end, u0_fn)
  prob = ODEProblem((du,u,p,t)->f!(du,u,m,t), m.u0, (0.0, t_end))
  sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

  # To account for cylindrical symmetry, total mass is integrated with weight r:
  #   mass ≈ ∫_{rmin}^{rmax} u(r) * r dr
  # We approximate this integral using midpoint rule over Ne elements:
  deltarho = (rhomax - rhomin)/Ne
  rho_s = rhomin .+ ((collect(1:Ne) .- 0.5) .* deltarho)
  mass0 = sum(u0_fn.(rho_s) .* rho_s) * deltarho
  mass1 = sum(sol.u[end] .* rho_s) * deltarho

  @test isapprox(mass1, mass0, atol=1e-6)
end