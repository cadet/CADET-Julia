
# Defines a Model wrapping
#   1) a Radial DG operator (convective/dispersive)
#   2) a linear rate mass sink (-kf * u )
# and provides f!(du,u,model,t) to plug straight into ODEProblem.
# -------------------------------------------------------------------
# the Model: operator + initial condition
# -------------------------------------------------------------------
using .RadialConvDispOperatorDG: RadialConvDispOperator, radial_convdisp_op!
using RadialDGElements: lglnodes
#using .binding_base
mutable struct Model{T}
    op::RadialConvDispOperator{T}
    u0::Vector{T}
    kf::T
    #binding
    #nComp::Int
    #bindStride::Int
end

"""
    Model(N, Ne, rhomin, rhomax, D_rad, tau, vel_fn, u0_fn)

Builds a radial DG + LRM model with:
- `N`        : polynomial degree per element
- `Ne`       : number of elements
- `rhomin,rhomax`: radial domain
- `D_rad`       : axial dispersion coefficient
- `tau`        : LRM time constant
- `vel_fn(rho)`: velocity profile
- `u0_fn(rho)` : initial condition profile

Returns `Model{Float64}` with precomputed operator and u₀ vector.
"""
function Model(N::Integer, Ne::Integer,
               rhomin::Real, rhomax::Real,
               D_rad::Real, tau::Real,
               vel_fn::Function, kf::Real,
               u0_fn::Function)
    # build the DG operator
    op = RadialConvDispOperator(N, Ne, rhomin, rhomax, D_rad, tau, vel_fn)

    # Generate initial condition at DG nodes in each element (mapped from reference element to physical radius)
    deltarho = (rhomax - rhomin)/Ne
    xi, _ = lglnodes(N)  # Nodal points on reference element

    u0 = zeros(Float64, N, Ne)
    for e in 1:Ne
        r_left = rhomin + (e-1)*deltarho
        r_nodes = ((xi .+ 1)*(deltarho/2)) .+ r_left
        u0[:,e] = u0_fn.(r_nodes)
    end
    u0 = vec(u0)  # flatten for ODE

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
    # Apply DG convective-dispersive operator: computes radial flux divergence term (in cylindrical coordinates) and radial diffusion
    radial_convdisp_op!(du, u, m.op, t)
    # Add linear rate mass sink: -kf * u term representing first-order decay
    @. du -= m.kf * u

    # Add binding source term
    #RHS_q = similar(u)
    # For Linear, qq can be zero if not using actual solid-phase dynamics
    #qq = zeros(size(u))  # Or pass the true stationary phase vector if you have a two-phase model
    #compute_binding!(RHS_q, u, qq, m.binding, m.nComp, m.bindStride, t)
    #du .-= RHS_q   # *subtract* binding rate: mobile phase loses mass as it binds

    return nothing
end