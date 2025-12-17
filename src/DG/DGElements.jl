#Functions used for DG Elements from Jan
module DGElements
using SpecialFunctions, LinearAlgebra

"""
    lglnodes(N)

Computes the Legendre-Gauss-Lobatto (LGL) nodes and weights for polynomial degree `N`.

# Arguments
- `N`: Polynomial degree.

# Returns
A tuple `(x, w)` where `x` are the LGL nodes and `w` are the corresponding weights.
"""
function lglnodes(N)
    # Truncation + 1
    N1 = N + 1

    # Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x = cos.(pi .* (0:N) ./ N)

    # The Legendre Vandermonde Matrix
    P = zeros(N1, N1)

    # Compute P_(N) using the recursion relation
    # Compute its first and second derivatives and
    # update x using the Newton-Raphson method.
    xold = ones(length(x)) .* 2.0

    while maximum(abs.(x .- xold)) > eps()
        xold[:] = x

        P[:, 1] .= 1
        P[:, 2] = x

        for k = 2:N
            P[:, k + 1] = ((2 * k - 1) .* x .* P[:, k] - (k - 1) .* P[:, k - 1]) ./ k
        end

        x[:] = xold - (x .* P[:, N1] - P[:, N]) ./ (N1 .* P[:, N1])
    end

    w = 2 ./(N .* N1 .* P[:, N1].^2)
    
    return reverse(x), 1 ./ w
end


"""
    barycentricWeights(baryWeights, _polyDeg, _nodes)

Computes barycentric weights.

# Arguments
- `baryWeights`: Vector to store barycentric weights (should be initialized with ones).
- `_polyDeg`: Polynomial degree.
- `_nodes`: Vector of interpolation nodes.

# Returns
The vector of barycentric weights.
"""
function barycentricWeights(baryWeights,_polyDeg,_nodes)
    for j in 1:_polyDeg
        for k in 0:j-1
            baryWeights[k+1] *= (_nodes[k+1] - _nodes[j+1]) * 1.0
            baryWeights[j+1] *= (_nodes[j+1] - _nodes[k+1]) * 1.0
        end
    end
    for j in 1:_polyDeg+1
        baryWeights[j] = 1 / baryWeights[j]
    end
    return baryWeights
end

"""
    derivativeMatrix(_polyDeg, _nodes)

Computes the nodal (Lagrange) polynomial derivative matrix for a given set of nodes.

# Arguments
- `_polyDeg`: Polynomial degree.
- `_nodes`: Vector of interpolation nodes.

# Returns
A matrix of derivatives of the Lagrange basis polynomials at the nodes.
"""
function derivativeMatrix(_polyDeg,_nodes)
    baryWeights = ones(_polyDeg + 1)
    baryWeights = barycentricWeights(baryWeights,_polyDeg,_nodes)
    _polyDerM = zeros(Float64,_polyDeg+1,_polyDeg+1)
    for i in 1:_polyDeg+1
        for j in 1:_polyDeg+1
            if i != j
                _polyDerM[i, j] = baryWeights[j] / (baryWeights[i] * (_nodes[i] - _nodes[j]))
                _polyDerM[i, i] += -_polyDerM[i, j]
            end
        end
    end
    return _polyDerM
end


"""
    orthonFactor(polyDeg, a=0.0, b=0.0)

Computes the normalization factor for Jacobi or Legendre polynomials.

# Arguments
- `polyDeg`: Polynomial degree.
- `a`, `b`: Jacobi polynomial parameters (default 0.0 for Legendre).

# Returns
The normalization factor as a Float64.
"""
function orthonFactor(polyDeg, a = 0.0, b = 0.0)
    # a = alpha, b = beta
    n = polyDeg

    return sqrt(((2.0 * n + a + b + 1.0) * gamma(n + 1.0) * gamma(n + a + b + 1.0)) / (2.0^(a + b + 1.0) * gamma(n + a + 1.0) * gamma(n + b + 1.0)))
end

"""
    jacVandermondeMatrix(_nodes, _polyDeg, alpha=0.0, beta=0.0)

Calculates the Vandermonde matrix of normalized Jacobi polynomials at the given nodes.

# Arguments
- `_nodes`: Vector of interpolation nodes.
- `_polyDeg`: Polynomial degree.
- `alpha`, `beta`: Jacobi polynomial parameters (optional).

# Returns
The Vandermonde matrix as a 2D array.
"""
function jacVandermondeMatrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    
    V = zeros(Float64, length(_nodes), length(_nodes))
    
    # degree 0
    V[:, 1] .= orthonFactor(0, alpha, beta)
    
    # degree 1
    V[:, 2] .= ((_nodes .- 1) ./ 2 .* (alpha .+ beta .+ 2) .+ (alpha .+  1.0)) .* orthonFactor(1, alpha, beta)


    for deg in 2:_polyDeg
        for node in 1:size(_nodes)[1] #length(_nodes)
            orthn_1 = orthonFactor(deg, alpha, beta) / orthonFactor(deg - 1, alpha, beta)
            orthn_2 = orthonFactor(deg, alpha, beta) / orthonFactor(deg - 2, alpha, beta)

            V[node, deg + 1] = orthn_1 * ((2.0 * deg + alpha + beta - 1.0) * ((2.0 * deg + alpha + beta) * (2.0 * deg + alpha + beta - 2.0) * _nodes[node] + alpha * alpha - beta * beta) * V[node, deg])
            V[node, deg + 1] -= orthn_2 * (2.0 * (deg + alpha - 1.0) * (deg + beta - 1.0) * (2.0 * deg + alpha + beta) * V[node, deg - 1]) 
            V[node, deg + 1] /= (2.0 * deg * (deg + alpha + beta) * (2.0 * deg + alpha + beta - 2.0))
 
        end
    end
    
    return V
end

"""
    getVandermonde_LEGENDRE(_nodes, _polyDeg, alpha=0.0, beta=0.0)

Calculates the Vandermonde matrix of normalized Legendre polynomials at the given nodes.

# Arguments
- `_nodes`: Vector of interpolation nodes.
- `_polyDeg`: Polynomial degree.
- `alpha`, `beta`: Parameters (optional, default to Legendre).

# Returns
The Vandermonde matrix as a 2D array.
"""
function getVandermonde_LEGENDRE(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    
    V = zeros(Float64, length(_nodes), length(_nodes))
    
    # degree 0
    V[:, 1] .= orthonFactor(0, alpha, beta)
    
    # degree 1
    V[:, 2] .= _nodes .* orthonFactor(1, alpha, beta)


    for deg in 2:_polyDeg
        for node in 1:size(_nodes)[1] #length(_nodes)
            orthn_1 = orthonFactor(deg, alpha, beta) / orthonFactor(deg - 1, alpha, beta)
            orthn_2 = orthonFactor(deg, alpha, beta) / orthonFactor(deg - 2, alpha, beta)

            fac_1 = ((2.0 * deg - 1.0) * 2.0 * deg * (2.0 * deg - 2.0) * _nodes[node]) / (2.0 * deg * deg * (2.0 * deg - 2.0))
            fac_2 = (2.0 * (deg - 1.0) * (deg - 1.0) * 2.0 * deg) / (2.0 * deg * deg * (2.0 * deg - 2.0))

            V[node, deg + 1] = orthn_1 * fac_1 * V[node, deg] - orthn_2 * fac_2 * V[node, deg - 1]
        end
    end
    
    return V
end

"""
    invMMatrix(_nodes, _polyDeg, alpha=0.0, beta=0.0)

Computes the inverse mass matrix for the given nodes and polynomial degree.

# Arguments
- `_nodes`: Vector of interpolation nodes.
- `_polyDeg`: Polynomial degree.
- `alpha`, `beta`: Jacobi polynomial parameters (optional).

# Returns
The inverse mass matrix as a 2D array.
"""
function invMMatrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    _invMM = jacVandermondeMatrix(_nodes,_polyDeg, alpha, beta) * transpose(jacVandermondeMatrix(_nodes,_polyDeg, alpha, beta))
    return _invMM
end

"""
    MMatrix(_nodes, _polyDeg, alpha=0.0, beta=0.0)

Computes the mass matrix for the given nodes and polynomial degree.

# Arguments
- `_nodes`: Vector of interpolation nodes.
- `_polyDeg`: Polynomial degree.
- `alpha`, `beta`: Jacobi polynomial parameters (optional).

# Returns
The mass matrix as a 2D array.
"""
function MMatrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    return inv(invMMatrix(_nodes,_polyDeg, alpha, beta))
end

# #Stiffness matrix
# function steifMatrix(_nodes,_polyDeg)
#     #Mass matrix * Derivative matrix
#     _steifMatrix = (getVandermonde_LEGENDRE(_nodes,_polyDeg) * transpose(getVandermonde_LEGENDRE(_nodes,_polyDeg)))^-1 * derivativeMatrix(_polyDeg,_nodes) 
#     return _steifMatrix
# end

"""
    second_order_stiff_matrix(_nodes, _polyDeg, alpha=0.0, beta=0.0)

Computes the second-order stiffness matrix for the given nodes and polynomial degree.

# Arguments
- `_nodes`: Vector of interpolation nodes.
- `_polyDeg`: Polynomial degree.
- `alpha`, `beta`: Jacobi polynomial parameters (optional).

# Returns
The second-order stiffness matrix as a 2D array.
"""
function second_order_stiff_matrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    #Mass matrix * Derivative matrix
    return Transpose(derivativeMatrix(_polyDeg,_nodes)) * MMatrix(_nodes,_polyDeg, alpha, beta) * derivativeMatrix(_polyDeg,_nodes)
end

# Legendre-Gauss (LG) nodes and weights
function lgnodes(N)
    N1 = N + 1

    x = zeros(N1)
    for i in 1:N1
        x[i] = -cos((2*i - 1) * pi / (2*N1))
    end

    # Legendre polynomial values
    P = zeros(N1, N1 + 1)

    # Newton-Raphson iteration to find roots of P_{N+1}(x)
    xold = ones(N1) .* 2.0

    while maximum(abs.(x .- xold)) > eps()
        xold[:] = x

        P[:, 1] .= 1.0
        P[:, 2] = x

        # Three-term recurrence for Legendre polynomials
        for k = 2:N1
            P[:, k + 1] = ((2*k - 1) .* x .* P[:, k] - (k - 1) .* P[:, k - 1]) ./ k
        end

        # Newton-Raphson update: x_{n+1} = x_n - P_{N+1}(x_n) / P'_{N+1}(x_n)
        # where P'_{N+1}(x) = (N+1) / (x^2 - 1) * [x * P_{N+1}(x) - P_N(x)]
        for i in 1:N1
            dP = (N1) / (x[i]^2 - 1.0) * (x[i] * P[i, N1 + 1] - P[i, N1])
            x[i] = xold[i] - P[i, N1 + 1] / dP
        end
    end

    # Compute weights: w_i = 2 / [(1 - x_i^2) * (P'_{N+1}(x_i))^2]
    w = zeros(N1)
    for i in 1:N1
        dP = (N1) / (x[i]^2 - 1.0) * (x[i] * P[i, N1 + 1] - P[i, N1])
        w[i] = 2.0 / ((1.0 - x[i]^2) * dP^2)
    end

    return x, 1 ./ w
end

# Chebyshev-Gauss-Lobatto nodes and weights
function cglnodes(N)
    N1 = N + 1

    # Chebyshev-Gauss-Lobatto nodes: x_j = -cos(π*j/N)
    x = zeros(N1)
    for j in 0:N
        x[j+1] = -cos(pi * j / N)
    end

    # Weights for Chebyshev-Gauss-Lobatto quadrature
    w = zeros(N1)
    w[1] = pi / (2 * N)
    w[N1] = pi / (2 * N)
    for j in 1:(N-1)
        w[j+1] = pi / N
    end

    return x, 1 ./ w
end

# Chebyshev-Gauss nodes and weights
function cgnodes(N)
    # N+1 interior nodes
    N1 = N + 1

    # Chebyshev-Gauss nodes: x_j = -cos(π*(2j+1)/(2N+2))
    x = zeros(N1)
    for j in 0:N
        x[j+1] = -cos(pi * (2*j + 1) / (2 * N1))
    end

    # Weights for Chebyshev-Gauss quadrature (all equal)
    w = fill(pi / N1, N1)

    return x, 1 ./ w
end

# hatrho Weighted mass matrix and its inverse
function weightedMMatrix(_nodes,_polyDeg, rho_i::Vector{Float64}, _deltarho::Float64)
    # Precompute base mass matrices once
    M00 = MMatrix(_nodes, _polyDeg, 0, 0)
    M01 = MMatrix(_nodes, _polyDeg, 0, 1)

    # Compute weighted mass matrix and its inverse for each cell
    nCells = length(rho_i) - 1
    rMM = Vector{Matrix{Float64}}(undef, nCells)
    invrMM = Vector{Matrix{Float64}}(undef, nCells)

    for Cell in 1:nCells
        # M_ρ = ρ_i * M(0,0) + (Δρ/2) * M(0,1)
        rMM[Cell] = rho_i[Cell] .* M00 .+ (_deltarho/2) .* M01
        # Precompute inverse
        invrMM[Cell] = inv(rMM[Cell])
    end

    return rMM, invrMM
end

function dispMMatrix(_nodes, _polyDeg, rho_i::Vector{Float64}, _deltarho::Float64, d_rad::Union{Float64, Function}, _polyDerM::Matrix{Float64}, rMM::Vector{Matrix{Float64}})
    nCells = length(rho_i) - 1
    nNodes = _polyDeg + 1
    S_g = Vector{Matrix{Float64}}(undef, nCells)

    # Check if d_rad is a constant or a function
    if isa(d_rad, Float64)
        # Use analytical formula for constant d_rad: S_g = d_rad * D^T * M_ρ
        # This computes: S_g[j,k] = ∫ (dL_j/dρ) * ρ * D * L_k dρ = D * D^T * M_ρ
        for Cell in 1:nCells
            S_g[Cell] = d_rad * transpose(_polyDerM) * rMM[Cell]
        end
    else
        # Use quadrature integration for spatially varying d_rad
        quad_nodes, quad_invWeights = lgnodes(_polyDeg)
        quad_weights = 1.0 ./ quad_invWeights
        nQuad = length(quad_nodes)

        lagrange_at_quad = zeros(Float64, nNodes, nQuad)
        for k in 1:nNodes
            for q in 1:nQuad
                lagrange_at_quad[k, q] = 1.0
                for m in 1:nNodes
                    if m != k
                        lagrange_at_quad[k, q] *= (quad_nodes[q] - _nodes[m]) / (_nodes[k] - _nodes[m])
                    end
                end
            end
        end

        dlagrange_at_quad = zeros(Float64, nNodes, nQuad)
        for q in 1:nQuad
            for j in 1:nNodes
                sum_terms = 0.0
                for n in 1:nNodes
                    if n != j
                        term = 1.0 / (_nodes[j] - _nodes[n])
                        for p in 1:nNodes
                            if p != j && p != n
                                term *= (quad_nodes[q] - _nodes[p]) / (_nodes[j] - _nodes[p])
                            end
                        end
                        sum_terms += term
                    end
                end
                dlagrange_at_quad[j, q] = sum_terms
            end
        end

        for Cell in 1:nCells
                S_g[Cell] = zeros(Float64, nNodes, nNodes)
                jacobian = _deltarho / 2  # dρ/dξ

                for q in 1:nQuad
                    rho_q = rho_i[Cell] + jacobian * (1 + quad_nodes[q])
                    d_rad_q = d_rad(rho_q)
                    # We need: S_g[j,k] = ∫ ρ*D*(dL_j/dρ)*L_k dρ
                    # dL/dρ = (dL/dξ) / J where J = Δρ/2
                    # Transform to ref coords: ∫ ρ*D*(dL/dξ / J)*L * J dξ = ∫ ρ*D*(dL/dξ)*L dξ
                    # The Jacobians cancel! So we integrate dL/dξ in reference space with weight ρ*D
                    weight_factor = quad_weights[q] * rho_q * d_rad_q
                    S_g[Cell] .+= weight_factor .* (dlagrange_at_quad[:, q] * lagrange_at_quad[:, q]')
                end
        end
    end
    return S_g
end

end
