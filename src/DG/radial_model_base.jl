
# Defines a Model wrapping
#   1) a Radial DG operator (convective/dispersive)
#   2) a linear rate mass sink (–u/tau)
# and provides f!(du,u,model,t) to plug straight into ODEProblem.
# -------------------------------------------------------------------
# the Model: operator + initial condition
# -------------------------------------------------------------------
# so that Main.RadialConvDispOperatorDG is guaranteed to exist
using .RadialConvDispOperatorDG: RadialConvDispOperator, radial_convdisp_op!
mutable struct Model{T}
    op::RadialConvDispOperator{T}
    u0::Vector{T}
    kf::T
end

"""
    Model(N, Ne, rhomin, rhomax, D_a, tau, vel_fn, u0_fn)

Builds a radial DG + LRM model with:
- `N`        : polynomial degree per element
- `Ne`       : number of elements
- `rhomin,rhomax`: radial domain
- `D_a`       : axial dispersion coefficient
- `tau`        : LRM time constant
- `vel_fn(rho)`: velocity profile
- `u0_fn(rho)` : initial condition profile

Returns `Model{Float64}` with precomputed operator and u₀ vector.
"""
function Model(N::Integer, Ne::Integer,
               rhomin::Real, rhomax::Real,
               D_a::Real, tau::Real,
               vel_fn::Function, kf::Real,
               u0_fn::Function)
    # build the DG operator
    op = RadialConvDispOperator(N, Ne, rhomin, rhomax, D_a, tau, vel_fn)

    # generate cell-centers for initial condition
    deltarho = (rhomax - rhomin)/Ne
    # cell‐centered nodes: midpoint of each element
    rhoc = collect(rhomin .+ deltarho/2 .+ deltarho*(0:(Ne-1)))

    # build u0
    u0 = [u0_fn(rho) for rho in rhoc]

    Model{eltype(u0)}(op, u0, kf)
end

# -------------------------------------------------------------------
# RHS function: du = convdisp_operator(u) – u/tau
# -------------------------------------------------------------------

"""
    f!(du, u, m::Model, t)

Computes the right-hand side of the ODE system:
    ∂u/∂t = L_op(u) - kf u
where L_op is the radial convective-dispersive operator (DG discretization), and -kf u is a linear rate mass sink (first-order decay).

Mathematically:
    du = m.op(u, t)   # radial convective-dispersive DG operator
    du .-= m.kf * u   # linear rate mass sink
"""
function f!(du::Vector{T}, u::Vector{T}, m::Model{T}, t::Real) where T
    # Zero the derivative vector: corresponds to initializing ∂u/∂t to zero
    fill!(du, zero(T))
    # Apply DG convective-dispersive operator: computes radial flux divergence term ∇·(ρ u v) + diffusion with axial dispersion
    radial_convdisp_op!(du, u, m.op, t)
    # Add linear rate mass sink: -kf * u term representing first-order decay
    @. du -= m.kf * u

    return nothing
end