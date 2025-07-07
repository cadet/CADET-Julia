"""
module RadialConvDispOperatorDG
using SpecialFunctions, LinearAlgebra
using Main.RadialDGElements: lglnodes, gauss_legendre, Vandermonde1D, invVandermonde1D, Dmatrix1D, MassMatrix1D, ExtractE1D, LiftMatrix1D

struct RadialConvDispOperator{T}
    N::Int
    Ne::Int
    rho_min::T
    rho_max::T
    D::T
    tau::T
    vel_fn::Function

    # reference element data
    nodes::Vector{T}       # Gauss-Legendre nodes in [-1,1]
    weights::Vector{T}     # Gauss-Legendre weights
    Vand::Matrix{T}
    invV::Matrix{T}
    Dr::Matrix{T}
    M::Matrix{T}
    E_int::Matrix{T}
    LIFT::Matrix{T}
end

function RadialConvDispOperator(N::Int, Ne::Int, rho_min::T, rho_max::T, D::T, tau::T, vel_fn::Function) where T
    # reference nodes & weights on [-1,1]
    xi, w = lglnodes(N)

    # Vandermonde and operators
    V    = Vandermonde1D(N, xi)
    invV = invVandermonde1D(xi)
    Dr   = Dmatrix1D(N, xi)
    M    = MassMatrix1D(w)
    E    = ExtractE1D(N)
    L    = LiftMatrix1D(w)

    return RadialConvDispOperator{T}(N, Ne, rho_min, rho_max, D, tau, vel_fn, xi, w, V, invV, Dr, M, E, L)
end

function radial_convdisp_op!(du::AbstractVector{T}, u::AbstractVector{T}, op::RadialConvDispOperator{T}, t) where T
    N, Ne = op.N, op.Ne
    rho_min, rho_max = op.rho_min, op.rho_max
    D, vel = op.D, op.vel_fn
    xi, w = op.nodes, op.weights
    Dr, M = op.Dr, op.M
    E, L  = op.E_int, op.LIFT

    # physical cell width
    # Physical cell width Δρ = (ρ_max - ρ_min) / Ne
    deltarho = (rho_max - rho_min)/Ne

    # loop over elements
    for e in 1:Ne
        # extract local solution
        u_e = @view u[(e-1)*(N+1)+1 : e*(N+1)]
        du_e = @view du[(e-1)*(N+1)+1 : e*(N+1)]

        # map xi-> physical r
        r_left = rho_min + (e-1)*deltarho
        r_e = ((xi .+ 1) .* (deltarho/2)) .+ r_left

        # convective flux f_c = v(r)*u
        f_c = vel.(r_e) .* u_e

        # diffusive flux f_d = -D * du/dr
        du_dr = Dr * u_e .* (2/deltarho)
        f_d = -D .* du_dr

        # total flux at nodes
        f_tot = f_c .+ f_d

        # volume integral term: elementwise scalar multiplication moved outside broadcast
        vol = Dr' * (M * f_tot) * (2/deltarho)

        # numerical flux at interfaces (simple upwind for conv, central for diff)
        # left interface between e-1 and e
        ul = u_e[1];  ur = u_e[end]
        # convective upwind
        vL = vel(r_e[1]);  vR = vel(r_e[end])
        FcL = vL >= 0 ? vL*ul : vL*ur
        FcR = vR >= 0 ? vR*ul : vR*ur
        # diffusive central
        FdL = -D*((u_e[2]-u_e[1])*(2/deltarho))
        FdR = -D*((u_e[end]-u_e[end-1])*(2/deltarho))

        # total interface flux
        FL = FcL + FdL
        FR = FcR + FdR

        # lift flux to element
        lift_term = L * vcat(FL, FR)

        # assemble du
        du_e .= -vol .+ lift_term
        du[(e-1)*(N+1)+1 : e*(N+1)] = du_e
    end

    return du
end

end # module
"""

module RadialConvDispOperatorDG
using SpecialFunctions, LinearAlgebra
using Main.RadialDGElements: lglnodes, gauss_legendre, Vandermonde1D, invVandermonde1D, Dmatrix1D, MassMatrix1D, ExtractE1D, LiftMatrix1D

struct RadialConvDispOperator{T}
    N::Int; Ne::Int; rho_min::T; rho_max::T; D::T; tau::T; vel_fn::Function
    nodes::Vector{T}; weights::Vector{T}; Vand::Matrix{T}; invV::Matrix{T}; Dr::Matrix{T}; M::Matrix{T}; E_int::Matrix{T}; LIFT::Matrix{T}
end

function RadialConvDispOperator(N::Int, Ne::Int, rho_min::T, rho_max::T, D::T, tau::T, vel_fn::Function) where T
    xi, w = lglnodes(N)
    V    = Vandermonde1D(N, xi)
    invV = invVandermonde1D(xi)
    Dr   = Dmatrix1D(N, xi)
    M    = MassMatrix1D(w)
    E    = ExtractE1D(N)
    L    = LiftMatrix1D(w)
    RadialConvDispOperator{T}(N, Ne, rho_min, rho_max, D, tau, vel_fn, xi, w, V, invV, Dr, M, E, L)
end

function radial_convdisp_op!(du::AbstractVector{T}, u::AbstractVector{T}, op::RadialConvDispOperator{T}, t) where T
    N, Ne = op.N, op.Ne
    rho_min, rho_max = op.rho_min, op.rho_max
    D, vel = op.D, op.vel_fn
    xi, w = op.nodes, op.weights
    Dr, M = op.Dr, op.M; E, L = op.E_int, op.LIFT
    deltarho = (rho_max - rho_min)/Ne

    ndof = N+1
    total_dof = length(u)
    nvar = total_dof ÷ (Ne*ndof)

    for e in 1:Ne
        # local slice for element and variables
        idx = ((e-1)*ndof*nvar + 1) : (e*ndof*nvar)
        u_loc = reshape(view(u, idx), ndof, nvar)
        du_loc = zeros(T, ndof, nvar)

        # Map reference nodes ξ ∈ [-1,1] to physical radii r = (ξ+1)*(Δρ/2) + ρ_left
        r_left = rho_min + (e-1)*deltarho
        r_e = ((xi .+ 1)*(deltarho/2)) .+ r_left

        # convective and diffusive per var
        for j in 1:nvar
            u_e = u_loc[:, j]
            # Convective flux f_c = v(r) * u
            vel_vec = vel.(r_e)
            f_c = vel_vec .* u_e
            # Compute spatial derivative du/dr ≈ D_r * u * (2/Δρ)
            # where Dr is derivative matrix on reference element and scaling factor accounts for mapping.
            du_dr = Dr * u_e * (2/deltarho)
            # Diffusive flux f_d = -D * (du/dr)
            f_d = -D * du_dr
            # total flux
            f_tot = f_c .+ f_d
            # Volume integral term: ∫ Dr^T M f_tot dξ scaled by (2/Δρ)
            vol = Dr' * (M * f_tot) * (2/deltarho)
            # Numerical interface fluxes:
            # - Convective: upwind based on sign of velocity.
            # - Diffusive: central flux approximation.
            ul, ur = u_e[1], u_e[end]
            vL, vR = vel(r_e[1]), vel(r_e[end])
            FcL = vL >= 0 ? vL*ul : vL*ur;  FcR = vR >= 0 ? vR*ul : vR*ur
            FdL = -D*((u_e[2]-u_e[1])*(2/deltarho)); FdR = -D*((u_e[end]-u_e[end-1])*(2/deltarho))
            FL, FR = FcL+FdL, FcR+FdR
            # Lift interface fluxes into element via LIFT matrix.
            lift_term = L * vcat(FL, FR)
            # Combine volume and surface terms to update time derivative du_loc.
            du_loc[:, j] = -vol .+ lift_term
        end

        # write back
        du[idx] .= reshape(du_loc, ndof*nvar)
    end
    du
end

end # module
