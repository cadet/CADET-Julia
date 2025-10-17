#Functions used for DG Elements from Jan
module DGElements
using SpecialFunctions, LinearAlgebra
using SpecialFunctions, LinearAlgebra

#LGL nodes and weights
function lglnodes(N)
    # Truncation + 1
    N1 = N + 1
#LGL nodes and weights
function lglnodes(N)
    # Truncation + 1
    N1 = N + 1

    # Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x = cos.(pi .* (0:N) ./ N)
    # Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x = cos.(pi .* (0:N) ./ N)

    # The Legendre Vandermonde Matrix
    P = zeros(N1, N1)
    # The Legendre Vandermonde Matrix
    P = zeros(N1, N1)

    # Compute P_(N) using the recursion relation
    # Compute its first and second derivatives and
    # update x using the Newton-Raphson method.
    xold = ones(length(x)) .* 2.0
    # Compute P_(N) using the recursion relation
    # Compute its first and second derivatives and
    # update x using the Newton-Raphson method.
    xold = ones(length(x)) .* 2.0

    while maximum(abs.(x .- xold)) > eps()
        xold[:] = x
    while maximum(abs.(x .- xold)) > eps()
        xold[:] = x

        P[:, 1] .= 1
        P[:, 2] = x
        P[:, 1] .= 1
        P[:, 2] = x

        for k = 2:N
            P[:, k + 1] = ((2 * k - 1) .* x .* P[:, k] - (k - 1) .* P[:, k - 1]) ./ k
        end
        for k = 2:N
            P[:, k + 1] = ((2 * k - 1) .* x .* P[:, k] - (k - 1) .* P[:, k - 1]) ./ k
        end

        x[:] = xold - (x .* P[:, N1] - P[:, N]) ./ (N1 .* P[:, N1])
    end
        x[:] = xold - (x .* P[:, N1] - P[:, N]) ./ (N1 .* P[:, N1])
    end

    w = 2 ./(N .* N1 .* P[:, N1].^2)
    
    return reverse(x), 1 ./ w
end
    w = 2 ./(N .* N1 .* P[:, N1].^2)
    
    return reverse(x), 1 ./ w
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

# Chebyshev polynomial evaluation using three-term recurrence (Algorithm 28 from Kopriva 2009)
function chebyshev_poly(x, N)
    # Evaluate Chebyshev polynomials T_0, T_1, ..., T_N at point x
    # Returns vector of length N+1
    T = zeros(N + 1)

    # Initial values
    T[1] = 1.0  # T_0(x) = 1
    if N >= 1
        T[2] = x  # T_1(x) = x
    end

    # Three-term recurrence: T_{n+1}(x) = 2*x*T_n(x) - T_{n-1}(x)
    for n in 1:(N-1)
        T[n+2] = 2.0 * x * T[n+1] - T[n]
    end

    return T
end

# Chebyshev polynomial first derivative evaluation
function chebyshev_poly_derivative(x, N)
    # Evaluate first derivatives T'_0, T'_1, ..., T'_N at point x
    # Returns vector of length N+1
    T = chebyshev_poly(x, N)
    dT = zeros(N + 1)

    # Initial values
    dT[1] = 0.0  # T'_0(x) = 0
    if N >= 1
        dT[2] = 1.0  # T'_1(x) = 1
    end

    # Recurrence relation: T'_{n+1}(x) = 2*T_n(x) + 2*x*T'_n(x) - T'_{n-1}(x)
    for n in 1:(N-1)
        dT[n+2] = 2.0 * T[n+1] + 2.0 * x * dT[n+1] - dT[n]
    end

    return dT
end

# Chebyshev polynomial second derivative evaluation
function chebyshev_poly_second_derivative(x, N)
    # Evaluate second derivatives T''_0, T''_1, ..., T''_N at point x
    # Returns vector of length N+1
    T = chebyshev_poly(x, N)
    dT = chebyshev_poly_derivative(x, N)
    d2T = zeros(N + 1)

    # Initial values
    d2T[1] = 0.0  # T''_0(x) = 0
    if N >= 1
        d2T[2] = 0.0  # T''_1(x) = 0
    end
    if N >= 2
        d2T[3] = 4.0  # T''_2(x) = 4
    end

    # Recurrence: T''_{n+1}(x) = 4*T'_n(x) + 2*x*T''_n(x) - T''_{n-1}(x)
    for n in 2:(N-1)
        d2T[n+2] = 4.0 * dT[n+1] + 2.0 * x * d2T[n+1] - d2T[n]
    end

    return d2T
end

# Chebyshev polynomial arbitrary derivative evaluation
function chebyshev_poly_nth_derivative(x, N, n_deriv)
    # Evaluate nth derivatives at point x
    # n_deriv: order of derivative
    if n_deriv == 0
        return chebyshev_poly(x, N)
    elseif n_deriv == 1
        return chebyshev_poly_derivative(x, N)
    elseif n_deriv == 2
        return chebyshev_poly_second_derivative(x, N)
    else
        # For higher derivatives, use the general recurrence
        error("Derivatives of order > 2 not yet implemented")
    end
end

# Chebyshev Vandermonde matrix (based on Kopriva 2009)
function getVandermonde_CHEBYSHEV(_nodes, _polyDeg)
    N1 = length(_nodes)
    V = zeros(Float64, N1, N1)

    # Fill Vandermonde matrix: V[i,j] = T_{j-1}(x_i)
    for i in 1:N1
        T = chebyshev_poly(_nodes[i], _polyDeg)
        V[i, :] = T
    end

    return V
end

# Normalized/Orthogonal Chebyshev Vandermonde matrix
function getVandermonde_CHEBYSHEV_NORMALIZED(_nodes, _polyDeg)
    N1 = length(_nodes)
    V = zeros(Float64, N1, N1)

    # Normalization factors for orthogonal Chebyshev polynomials
    # ∫_{-1}^{1} T_n(x) T_m(x) / √(1-x²) dx = 0 if n≠m, π/2 if n=m>0, π if n=m=0
    norm = zeros(N1)
    norm[1] = sqrt(pi)  # ||T_0||
    for n in 1:_polyDeg
        norm[n+1] = sqrt(pi / 2.0)  # ||T_n|| for n > 0
    end

    # Fill normalized Vandermonde matrix
    for i in 1:N1
        T = chebyshev_poly(_nodes[i], _polyDeg)
        for j in 1:N1
            V[i, j] = T[j] / norm[j]
        end
    end

    return V
end

# Chebyshev derivative matrix (Vandermonde approach)
function chebyshev_derivative_vandermonde(_nodes, _polyDeg)
    N1 = length(_nodes)
    DV = zeros(Float64, N1, N1)

    # Fill derivative Vandermonde matrix: DV[i,j] = T'_{j-1}(x_i)
    for i in 1:N1
        dT = chebyshev_poly_derivative(_nodes[i], _polyDeg)
        DV[i, :] = dT
    end

    return DV
end

# Chebyshev mass matrix (for Gauss or Gauss-Lobatto quadrature)
function chebyshev_mass_matrix(_nodes, _polyDeg)
    # M = V^T * V for normalized Chebyshev basis
    V = getVandermonde_CHEBYSHEV_NORMALIZED(_nodes, _polyDeg)
    return V' * V
end

# Chebyshev stiffness matrix
function chebyshev_stiffness_matrix(_nodes, _polyDeg)
    # S = V^T * D * V where D is the derivative matrix
    V = getVandermonde_CHEBYSHEV_NORMALIZED(_nodes, _polyDeg)
    D = derivativeMatrix(_polyDeg, _nodes)
    return V' * D * V
end


# computation of barycentric weights for fast polynomial evaluation
# @param [in] baryWeights vector to store barycentric weights. Must already be initialized with ones!
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
# computation of barycentric weights for fast polynomial evaluation
# @param [in] baryWeights vector to store barycentric weights. Must already be initialized with ones!
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

# @brief computation of nodal (lagrange) polynomial derivative matrix
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
# @brief computation of nodal (lagrange) polynomial derivative matrix
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


#  factor to normalize legendre polynomials
function orthonFactor(polyDeg, a = 0.0, b = 0.0)
    # a = alpha, b = beta
    n = polyDeg

    return sqrt(((2.0 * n + a + b + 1.0) * gamma(n + 1.0) * gamma(n + a + b + 1.0)) / (2.0^(a + b + 1.0) * gamma(n + a + 1.0) * gamma(n + b + 1.0)))
end
#  factor to normalize legendre polynomials
function orthonFactor(polyDeg, a = 0.0, b = 0.0)
    # a = alpha, b = beta
    n = polyDeg

    return sqrt(((2.0 * n + a + b + 1.0) * gamma(n + 1.0) * gamma(n + a + b + 1.0)) / (2.0^(a + b + 1.0) * gamma(n + a + 1.0) * gamma(n + b + 1.0)))
end

# calculates the Vandermonde matrix of the normalized jacobi polynomials
function jacVandermondeMatrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    
    V = zeros(Float64, length(_nodes), length(_nodes))
    
    # degree 0
    V[:, 1] .= orthonFactor(0, alpha, beta)
    
    # degree 1
    V[:, 2] .= ((_nodes .- 1) ./ 2 .* (alpha .+ beta .+ 2) .+ (alpha .+  1.0)) .* orthonFactor(1, alpha, beta)
# calculates the Vandermonde matrix of the normalized jacobi polynomials
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
            V[node, deg + 1] = orthn_1 * ((2.0 * deg + alpha + beta - 1.0) * ((2.0 * deg + alpha + beta) * (2.0 * deg + alpha + beta - 2.0) * _nodes[node] + alpha * alpha - beta * beta) * V[node, deg])
            V[node, deg + 1] -= orthn_2 * (2.0 * (deg + alpha - 1.0) * (deg + beta - 1.0) * (2.0 * deg + alpha + beta) * V[node, deg - 1]) 
            V[node, deg + 1] /= (2.0 * deg * (deg + alpha + beta) * (2.0 * deg + alpha + beta - 2.0))
 
        end
    end
    
    return V
end

# calculates the Vandermonde matrix of the normalized legendre polynomials
function getVandermonde_LEGENDRE(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    
    V = zeros(Float64, length(_nodes), length(_nodes))
    
    # degree 0
    V[:, 1] .= orthonFactor(0, alpha, beta)
    
    # degree 1
    V[:, 2] .= _nodes .* orthonFactor(1, alpha, beta)
# calculates the Vandermonde matrix of the normalized legendre polynomials
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
    for deg in 2:_polyDeg
        for node in 1:size(_nodes)[1] #length(_nodes)
            orthn_1 = orthonFactor(deg, alpha, beta) / orthonFactor(deg - 1, alpha, beta)
            orthn_2 = orthonFactor(deg, alpha, beta) / orthonFactor(deg - 2, alpha, beta)

            fac_1 = ((2.0 * deg - 1.0) * 2.0 * deg * (2.0 * deg - 2.0) * _nodes[node]) / (2.0 * deg * deg * (2.0 * deg - 2.0))
            fac_2 = (2.0 * (deg - 1.0) * (deg - 1.0) * 2.0 * deg) / (2.0 * deg * deg * (2.0 * deg - 2.0))
            fac_1 = ((2.0 * deg - 1.0) * 2.0 * deg * (2.0 * deg - 2.0) * _nodes[node]) / (2.0 * deg * deg * (2.0 * deg - 2.0))
            fac_2 = (2.0 * (deg - 1.0) * (deg - 1.0) * 2.0 * deg) / (2.0 * deg * deg * (2.0 * deg - 2.0))

            V[node, deg + 1] = orthn_1 * fac_1 * V[node, deg] - orthn_2 * fac_2 * V[node, deg - 1]
        end
    end
    
    return V
end
            V[node, deg + 1] = orthn_1 * fac_1 * V[node, deg] - orthn_2 * fac_2 * V[node, deg - 1]
        end
    end
    
    return V
end

#Inverse mass matrix
function invMMatrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    _invMM = jacVandermondeMatrix(_nodes,_polyDeg, alpha, beta) * transpose(jacVandermondeMatrix(_nodes,_polyDeg, alpha, beta))
    return _invMM
end
#Inverse mass matrix
function invMMatrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    _invMM = jacVandermondeMatrix(_nodes,_polyDeg, alpha, beta) * transpose(jacVandermondeMatrix(_nodes,_polyDeg, alpha, beta))
    return _invMM
end

# Mass matrix
function MMatrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    return inv(invMMatrix(_nodes,_polyDeg, alpha, beta))
end
# Mass matrix
function MMatrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    return inv(invMMatrix(_nodes,_polyDeg, alpha, beta))
end

# #Stiffness matrix
# function steifMatrix(_nodes,_polyDeg)
#     #Mass matrix * Derivative matrix
#     _steifMatrix = (getVandermonde_LEGENDRE(_nodes,_polyDeg) * transpose(getVandermonde_LEGENDRE(_nodes,_polyDeg)))^-1 * derivativeMatrix(_polyDeg,_nodes) 
#     return _steifMatrix
# end

# Second order stiffness matrix
function second_order_stiff_matrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    #Mass matrix * Derivative matrix
    return Transpose(derivativeMatrix(_polyDeg,_nodes)) * MMatrix(_nodes,_polyDeg, alpha, beta) * derivativeMatrix(_polyDeg,_nodes)
end
# Second order stiffness matrix
function second_order_stiff_matrix(_nodes,_polyDeg, alpha=0.0, beta=0.0)
    #Mass matrix * Derivative matrix
    return Transpose(derivativeMatrix(_polyDeg,_nodes)) * MMatrix(_nodes,_polyDeg, alpha, beta) * derivativeMatrix(_polyDeg,_nodes)
end

    # Weighted Mass Matrix
    # M_rho = (Δrho_i/2) * M^{(0,1)} + rho_i * M^(0,0), where hatrho(ξ) = rho_i + (Δrho_i/2)(1+ξ)
    function MrhoMatrix(_nodes, _polyDeg, _deltarho, rho_i)
        return (_deltarho/2) * MMatrix(_nodes, _polyDeg, 0.0, 1.0) + rho_i * MMatrix(_nodes, _polyDeg, 0.0, 0.0)
    end

    function invMrhoMatrix(_nodes,_polyDeg, _deltarho, rho_i)
        return inv(MrhoMatrix(_nodes, _polyDeg, _deltarho, rho_i))
    end

    # Weighted Matrix must accept a vector of physical radii ρ and return a vector of weights of same size
    function weightedQuadrature(_nodes, rho_i, _deltarho, weight)
        # Map reference nodes ξ∈[-1,1] to physical radii hatρ = ρ_i + (ξ+1)·Δρ_i/2
        hatrho = rho_i .+ (_nodes .+ 1) .* (_deltarho/2)
        # Fetch LGL quadrature weights for the same N and assume ordering equals `nodes`
        _, invw = lglnodes(length(_nodes) - 1)
        w = 1.0 ./ invw
        w_phys = weight(hatrho)
        return Diagonal(w .* w_phys)
    end

    function weighted_stiff_Matrix(_nodes, _polyDeg, rho_i, _deltarho, weight)
        return Transpose(derivativeMatrix(_polyDeg,_nodes)) * weightedQuadrature(_nodes, rho_i, _deltarho, hatrho -> hatrho .* weight.(hatrho))
    end

end
