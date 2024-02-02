#Functions used for DG Elements from Jan
module DGElements
using SpecialFunctions,LinearAlgebra


#Basic stuff

#LGL nodes and weights
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


#  factor to normalize legendre polynomials
function orthonFactor(polyDeg)
    n = polyDeg
    a = 0.0
    b = 0.0
    return sqrt(((2.0 * n + a + b + 1.0) * gamma(n + 1.0) * gamma(n + a + b + 1.0)) / (2.0^(a + b + 1.0) * gamma(n + a + 1.0) * gamma(n + b + 1.0)))
end

# calculates the Vandermonde matrix of the normalized legendre polynomials
function getVandermonde_LEGENDRE(_nodes,_polyDeg)
    
    V = zeros(Float64, length(_nodes), length(_nodes))

    alpha = 0.0
    beta = 0.0
    
    # degree 0
    V[:, 1] .= orthonFactor(0)
    
    # degree 1
    V[:, 2] .= _nodes .* orthonFactor(1)


    for deg in 2:_polyDeg
        for node in 1:size(_nodes)[1] #length(_nodes)
            orthn_1 = orthonFactor(deg) / orthonFactor(deg - 1)
            orthn_2 = orthonFactor(deg) / orthonFactor(deg - 2)

            fac_1 = ((2.0 * deg - 1.0) * 2.0 * deg * (2.0 * deg - 2.0) * _nodes[node]) / (2.0 * deg * deg * (2.0 * deg - 2.0))
            fac_2 = (2.0 * (deg - 1.0) * (deg - 1.0) * 2.0 * deg) / (2.0 * deg * deg * (2.0 * deg - 2.0))

            V[node, deg + 1] = orthn_1 * fac_1 * V[node, deg] - orthn_2 * fac_2 * V[node, deg - 1]
        end
    end
    
    return V
end

#Inverse mass matrix
function invMMatrix(_nodes,_polyDeg)
    _invMM = getVandermonde_LEGENDRE(_nodes,_polyDeg) * transpose(getVandermonde_LEGENDRE(_nodes,_polyDeg))
    return _invMM
end

# #Stiffness matrix
# function steifMatrix(_nodes,_polyDeg)
#     #Mass matrix * Derivative matrix
#     _steifMatrix = (getVandermonde_LEGENDRE(_nodes,_polyDeg) * transpose(getVandermonde_LEGENDRE(_nodes,_polyDeg)))^-1 * derivativeMatrix(_polyDeg,_nodes) 
#     return _steifMatrix
# end

# Second order stiffness matrix



end
