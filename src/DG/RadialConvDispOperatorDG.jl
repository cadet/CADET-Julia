module RadialConvDispOperatorDG
using SpecialFunctions, LinearAlgebra
using Main.RadialDGElements: lglnodes, gauss_legendre, Vandermonde1D, invVandermonde1D, Dmatrix1D, MassMatrix1D, ExtractE1D, LiftMatrix1D

struct RadialConvDispOperator{T} # Holds all DG matrices, geometry, and parameters for the radial conv-disp operator.
    N::Int; Ne::Int; rho_min::Float64; rho_max::Float64; D::T; tau::T; vel_fn::Function
    nodes::Vector{T}; weights::Vector{T}; Vand::Matrix{T}; invV::Matrix{T}; D_rad::Matrix{T}; M::Matrix{T}; E_int::Matrix{T}; LIFT::Matrix{T}
end

function RadialConvDispOperator(N::Int, Ne::Int, rho_min::T, rho_max::T, D::T, tau::T, vel_fn::Function) where T
    xi, w = lglnodes(N) # “Legendre-Gauss-Lobatto nodes in [-1,1]”
    V    = Vandermonde1D(N, xi)
    invV = invVandermonde1D(xi)
    D_rad   = Dmatrix1D(N, xi)
    M    = MassMatrix1D(w)
    E    = ExtractE1D(N)
    L    = LiftMatrix1D(w)
    RadialConvDispOperator{T}(N, Ne, rho_min, rho_max, D, tau, vel_fn, xi, w, V, invV, D_rad, M, E, L)
end

function radial_convdisp_op!(du::AbstractVector{T}, u::AbstractVector{T}, op::RadialConvDispOperator{T}, t) where T
    N, Ne = op.N, op.Ne
    rho_min, rho_max = op.rho_min, op.rho_max
    D, vel = op.D, op.vel_fn
    xi, w = op.nodes, op.weights
    D_rad, M = op.D_rad, op.M; E, L = op.E_int, op.LIFT
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
        rho_left = rho_min + (e-1)*deltarho
        rho_e = ((xi .+ 1)*(deltarho/2)) .+ rho_left # physical quadrature nodes
        
        # convective and diffusive per var
        for j in 1:nvar
            u_e = u_loc[:, j]
            # Convective flux f_c = v(r) * u
            vel_vec = vel.(rho_e)
            f_c = vel_vec .* u_e
            # Compute spatial derivative du/dr ≈ D_rad * u * (2/Δρ)
            # where D_rad is derivative matrix on reference element and scaling factor accounts for mapping.
            du_dr = D_rad * u_e * (2/deltarho)
            # Diffusive flux f_d = -D * (du/dr)
            f_d = -D * du_dr
            # total flux
            f_tot = f_c .+ f_d
            # Volume integral term: ∫ D_rad^T M f_tot dξ scaled by (2/Δρ)
            """Note: For full cylindrical accuracy, M should include a weighting by the physical radius r at each node."""
            # vol = D_rad' * (M * f_tot) * (2/deltarho)
            w_weighted = w .* r_e  # Quadrature weights times physical radius
            M_weighted = MassMatrix1D(w_weighted)  # rho_e are physical node positions in this element
            vol = D_rad' * (M_weighted * f_tot) * (2/deltarho)
            # Numerical interface fluxes:
            # - Convective: upwind based on sign of velocity.
            # - Diffusive: central flux approximation.
            ul, ur = u_e[1], u_e[end]
            vL, vR = vel(rho_e[1]), vel(rho_e[end])
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
