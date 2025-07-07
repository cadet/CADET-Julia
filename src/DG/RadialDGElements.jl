module RadialDGElements
export gauss_legendre, lglnodes, Vandermonde1D, invVandermonde1D, Dmatrix1D, MassMatrix1D, ExtractE1D, LiftMatrix1D
using LinearAlgebra

"""
Compute the Legendre polynomial P_n(x) of degree n at point x.

The Legendre polynomials are defined by the three-term recurrence:
  P₀(x) = 1,
  P₁(x) = x,
  (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}(x),
for k ≥ 1. They satisfy orthogonality on [-1,1]:
  ∫_{-1}^1 P_m(x) P_n(x) dx = 2/(2n+1) δ_{mn}.
"""
function legendreP(n::Integer, x::Real)
    # recurrence: P₀(x)=1, P₁(x)=x
    if n == 0
        return one(x)
    elseif n == 1
        return x
    else
        p_nm2 = one(x)        # P₀
        p_nm1 = x             # P₁
        p_n = zero(x)
        for k in 2:n
            p_n = ((2k-1)*x*p_nm1 - (k-1)*p_nm2) / k
            p_nm2, p_nm1 = p_nm1, p_n
        end
        return p_n
    end
end

#
# Gauss–Legendre quadrature finds nodes x_i and weights w_i such that
# ∫_{-1}^1 f(x) dx ≈ Σ_{i=1}^N w_i f(x_i),
# with exactness for polynomials up to degree 2N-1.
# We compute nodes as eigenvalues of the symmetric tridiagonal Jacobi matrix
# from the Legendre recurrence, and weights from the first component of eigenvectors.
"""
Gauss-Legendre nodes and weights on [-1,1] for N points.
"""
function gauss_legendre(N::Integer)
    # Build symmetric tridiagonal matrix from recurrence
    beta = [k / sqrt(4k^2 - 1) for k in 1:N-1]
    J = SymTridiagonal(zeros(N), beta)
    F = eigen(J)
    x = F.values
    V = F.vectors
    w = 2.0 * (V[1, :].^2)
    return x, w
end

#
# Legendre–Gauss–Lobatto (LGL) quadrature includes endpoints ±1.
# Nodes are the roots of (1 − x²) P'_{N-1}(x), given by cosines,
# and weights w_i = 2/(N(N-1)[P_{N-1}(x_i)]²), exact for polynomials ≤ 2N-3.
"""
Legendre-Gauss-Lobatto (LGL) nodes and weights on [-1,1] for N points.
"""
function lglnodes(N::Integer)
    # uniformly spaced cosines for endpoints + interior
    x = [-cos(pi*(i-1)/(N-1)) for i in 1:N]
    # weights based on (N-1)th Legendre polynomial
    w = [2/(N*(N-1)*legendreP(N-1, xi)^2) for xi in x]
    return x, w
end

#
# Construct the Vandermonde matrix V where V_{i,j} = P_{j-1}(x_i)
# maps modal coefficients (in Legendre basis) to nodal values at x_i.
"""
Vandermonde matrix (1D) using Legendre basis up to degree N-1.
"""
function Vandermonde1D(N::Integer, x::AbstractVector{<:Real})
    V = [legendreP(j-1, xi) for xi in x, j in 1:N]
    return V
end

#
# Inverse Vandermonde matrix to recover modal coefficients from nodal values.
"""
Inverse Vandermonde matrix.
"""
function invVandermonde1D(x::AbstractVector{<:Real})
    N = length(x)
    V = Vandermonde1D(N, x)
    return inv(V)
end

#
# Differentiation matrix D approximates d/dx at nodes x_i:
# D_{i,j} = P'_{N-1}(x_i) / [P'_{N-1}(x_j) (x_i − x_j)] for i ≠ j,
# with diagonal entries chosen so each row sums to zero.
"""
Differentiation matrix on the spectral nodes.
"""
function Dmatrix1D(N::Integer, x::AbstractVector{<:Real})
    D = zeros(N, N)
    for i in 1:N, j in 1:N
        if i != j
            D[i, j] = legendreP(N-1, x[i]) / (legendreP(N-1, x[j]) * (x[i] - x[j]))
        end
    end
    for i in 1:N
        D[i, i] = -sum(D[i, k] for k in 1:N if k != i)
    end
    return D
end

#
# Mass matrix for 1D DG is diagonal with quadrature weights w_i:
# M_{ii} = ∫_{-1}^1 ℓ_i(x) ℓ_i(x) dx ≈ w_i, where ℓ_i are Lagrange basis.
"""
Mass matrix (diagonal) from quadrature weights.
"""
function MassMatrix1D(w::AbstractVector{<:Real})
    return Diagonal(w)
end

#
# Face extraction operator picks left (i=1) and right (i=N) nodal values:
# E_{:,k} = ℓ at boundary k.
"""
Face extraction operator (left, right) as a Nx2 matrix.
"""
function ExtractE1D(N::Integer)
    eL = zeros(N); eL[1] = 1.0
    eR = zeros(N); eR[N] = 1.0
    return hcat(eL, eR)
end

#
# Lift matrix L = M^{-1} E spreads boundary flux contributions into volume
# ensuring the DG formulation conserves mass and enforces flux through faces.
"""
Lift operator L = M^{-1} E for flux lifting.
"""
function LiftMatrix1D(w::AbstractVector{<:Real})
    M = MassMatrix1D(w)
    E = ExtractE1D(length(w))
    return M \ E
end

end # module
