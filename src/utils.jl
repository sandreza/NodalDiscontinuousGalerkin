
using Revise
using SpecialFunctions # for gamma function reasons
using LinearAlgebra    # for Guass quadrature

"""
unimesh1D(xmin, xmax, K)

# Description

    Generates a uniform 1D mesh

# Arguments

    xmin: smallest value of array

    xmax: largest values of array

    K: number of elements in an array

# Return Values: VX, EtoV

    VX: vertex values | an Array of size K+1

    EtoV: element to node connectivity | a Matrix of size Kx2

# Example
xmin = -1
xmax =  1
K    =  4
VX, EtoV = unimesh1D(xmin, xmax, K)

"""
function unimesh1D(xmin, xmax, K)
    VX = @. collect(0:K) / K * (xmax - xmin) + xmin
    EtoV = Int.(ones(K, 2))
    for i = 1:K
        EtoV[i,1] = Int(i)
        EtoV[i,2] = Int(i+1)
    end
    return VX, EtoV
end

#code checked against the matlab code
"""
jacobi(x, α, β, n)

# Description

- Evaluates the jacobi polynomial at the point x

# Arguments

- `x`: point at which you will evaluate the jacobi polynomial
- `α`: first parameter for Jacobi polynomials
- `β`: second parameter for Jacobi polynomials
- `n` : order

# Return

-  `y`: the value of the of the Jacobi polynomial

"""
function jacobi(x, α, β, n::Int)
    γ0 = 2^(α + β + 1) / (α + β + 1) * gamma(α+1) * gamma(β+1)
    γ0 /= gamma(α+β+1)
    γ1 = (α + 1) * (β + 1) / (α + β + 3) * γ0
    #create array to return
    PL = zeros(n+1)
    PL[1] = 1 / sqrt(γ0)
    if n==0
        return PL[end]
    elseif n==1
        PL[2] = ( (α + β + 2) * x / 2 + (α - β) / 2) / sqrt(γ1)
        return PL[end]
    else
        PL[2] = ( (α + β + 2) * x / 2 + (α - β) / 2) / sqrt(γ1)
        aold = 2 / (2 + α + β) * sqrt((α+1)*(β+1)/(α + β + 3))
        for i in 1:(n-1)
            h1 = 2 * i + α + β
            anew = 2 /(h1 + 2)*sqrt((i+1)*(i+1+α+β)*(i+1+α)*(i+1+β)/(h1+1)/(h1+3))
            bnew = - (α^2 - β^2)/h1/(h1+2)
            PL[i+2] = 1 / anew * (-aold*PL[i] + (x-bnew)*PL[i+1])
            aold = anew
        end
        return PL[end]
    end
end

"""
djacobi(x, α, β, n)

# Description

- Evaluates the derivative of the jacobi polynomial at the point x

# Arguments

- `x`: point at which you will evaluate the derivative of the jacobi polynomial
- `α`: first parameter for Jacobi polynomials
- `β`: second parameter for Jacobi polynomials
- `n` : order

# Return

-  `y`: the derivative of the of the Jacobi polynomial

"""
function djacobi(x, α, β, n::Int)
    if n==0
        dp = 0.0
        return dp
    end
    dp = sqrt(n * (n + α + β + 1)) * jacobi(x, α + 1, β + 1, n-1)
    return dp
end



"""
vandermonde(x, α, β, N)

# Description

    Return vandermonde matrix of order N at the values x
    Allocates a little bit of memory

# Arguments

-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second parameter for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include

# Return Values

-   `v`: vandermonde matrix

# Example

    See LegendreTests.jl

"""
function vandermonde(x, α, β, N)
    # compute first two coefficients
    γ0 = 2^(α + β + 1) * factorial(α) * factorial(β) / ((α + β + 1) * factorial(α + β))
    γ1 = (α + 1) * (β + 1) / (α + β + 3) * γ0

    # create view to assign values
    v = zeros(length(x), N+1)
    v1 = view(v, :, 1)
    @. v1 = 1 / sqrt(γ0)

    # explicitly compute second coefficient
    if N == 0
        return v
    end

    v2 = view(v, :, 2)
    @. v2 = ( (α + β + 2) * x/2 + (α - β)/2) / sqrt(γ1)

    if N == 1
        return v
    end

    aʲ = 2 / (2 + α + β) * sqrt((α+1) * (β+1) / (α + β + 3))

    for i in 3:(N+1)
        # get views for ith, i-1th, and i-2th columns
        vi = view(v, :, i)
        vM1 = view(v, :, i-1)
        vM2 = view(v, :, i-2)

        # compute new a and b values
        h1 = 2 * (i-2) + α + β
        aⁱ = 2 / (h1 + 2) * sqrt((i-1) * (i-1 + α + β) * (i-1 + α) * (i-1 + β) / ((h1 + 1) * (h1 + 3)))
        bⁱ = - (α^2 - β^2) / (h1 * (h1 + 2))

        # compute coefficients for ith column
        @. vi = 1 / aⁱ * (-aʲ * vM2 + (x - bⁱ) * vM1)

        # save a coefficient for next iteration
        aʲ = aⁱ
    end

    return v
end

"""
dvandermonde(x, α, β, N)

# Description

    Return the gradient of the vandermonde matrix of order N at the values x
    Allocates a little bit of memory

# Arguments

-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second paramater for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include

# Return Values

-   `vr`: gradient of vandermonde matrix

# Example

    See LegendreTests.jl

"""
function dvandermonde(x, α, β, N)
    # create empty matrix (also handles first set of derivatives)
    vr = zeros(length(x), N+1)

    if N == 0
        return vr
    end

    # set values using vandermonde matrix
    v = vandermonde(x, α+1, β+1, N)
    for i in 1:N
        vi = view(v, :, i)
        vrP1 = view(vr, :, i+1)
        @. vrP1 = sqrt(i * (α + β + i+1)) * vi
    end

    return vr
end

"""
dmatrix(x, α, β, N)

# Description

    Return the differentiation matrix of order N at the values x
    Allocates too much memory

# Arguments

-   `x`: points at which to evaluate the Jacobi polynomials
-   `α`: first parameter for Jacobi polynomials
-   `β`: second paramater for Jacobi polynomials
-   `N`: maximum order of Jacobi polynomial to include

# Return Values

-   `D`: the differentiation matrix

# Example

    See LegendreTests.jl

"""
function dmatrix(x, α, β, N)
    # calculate vandermonde matrix and grad of vandermonde matrix
    vr = dvandermonde(x, α, β, N)
    v  =  vandermonde(x, α, β, N)

    # calculate values using D = vr * v^-1
    d = vr / v

    return d
end

"""
lift1D(V, y)
for computing fluxes

helps compute a surface integral of a quantity
note that the parantheses are necessary to prevent too much multiplcation
the E function takes the surface integrals are presents it
with respect to the full space inside an element
the entire operator represents how fluxes flow
into the interior of an element
"""
function lift1D(V)
    m,n = size(V)

    E = zeros(m , 2)
    E[1,1] = 1.0
    E[m,2] = 1.0

    return V * (transpose(V) * E)
end

"""
lift1D_v2(V, y)
for computing fluxes

nodal form
helps compute a surface integral of a quantity
note that the parantheses are necessary to prevent too much multiplcation
the E function takes the surface integrals are presents it
with respect to the full space inside an element
the entire operator represents how fluxes flow
into the interior of an element
"""
function lift1D_v2(V)
    m,n = size(V)

    E = zeros(m , 2)
    E[1,1] = 1.0
    E[m,2] = 1.0

    return E
end


"""
jacobiGQ(α, β, N)
# Description
    Guass Quadrature points and weights for the Jacobi Polynomial (α,β)
# Input
α, β: Jacobi polynomial descriptors
N:    order of quadrature points
# Return: x,w
x: quadrature points | array of size N+1
w: quadrature weights | array of size N+1

#Example
α = 0
β = 0
N = 4
x, w = jacobiGQ(α, β, N)
"""

function jacobiGQ(α, β, N)
    x = zeros(N+1)
    w = zeros(N+1)

    # explicit if N=0
    if N == 0
        x[1] = (α - β) / (α + β + 2)
        w[1] = 2
    end

    # form symmetric matrix from recurrence
    h1 = @. 2 * collect(0:N) + α + β;

    # construct diagonal matrix
    diag = @. - (α^2 - β^2) / (h1 + 2) / h1

    # construct super diagonal matrix
    h1view = view(h1, 1:N)
    cf = collect(1:N) # common factor that shows up a lot
    superdiag = @. 2 / (h1view + 2)
    @. superdiag *= sqrt( cf * (cf + α + β) * (cf + α) * (cf + β) )
    @. superdiag *= sqrt( 1 / (h1view + 1) / (h1view + 3) )

    # create full matrix combining the two
    J = SymTridiagonal(diag, superdiag)
    if (α + β) ≈ 0.0
        J[1,1] = -0.0
    end

    # compute quadrature by eigenvalue solve
    x,V = eigen(J)
    w = @. (V[1,:] ^ 2 ) * 2^(α + β + 1) / (α + β + 1)
    @. w *= factorial(α) * factorial(β) / factorial(α + β)
    return x,w
end

"""
jacobiGL(α, β, N)
# Description

    Guass Labatto quadrature points for the Jacobi Polynomial (α,β)
    The quadrature weights are computed as well (but not returned)

# Arguments

- `α, β`: Jacobi polynomial descriptors
- `N`:    order of quadrature


# Return: x

- `x`: quadrature points  | array of size N+1

# Examples
```julia-repl
julia> α = 0
0
julia> β = 0
0
julia> N = 4
4
julia> x = jacobiGL(α, β, N)
5-element Array{Float64,1}:
 -1.0
 -0.6546536707079759
  4.440892098500626e-16
  0.6546536707079771
  1.0
```
"""
function jacobiGL(α, β, N)
    x = zeros(N+1)
    w = zeros(N+1)

    # set end points
    x[1] = -1.0
    x[end] =  1.0

    # need at least two nodes
    if N == 0
        error("What are you doing?")
    end

    # have exactly two nodes
    if N == 1
        return x
    end

    # compute inner nodes
    xview = view(x, 2:N)
    xtmp,w = jacobiGQ(α+1, β+1, N-2)
    @. xview = xtmp
    return x
end



"""
Some nice documentation here.

# Examples
```jldoctest
julia> a = [1 2; 3 4]
2×2 Array{Int64,2}:
 1  2
 3  4
```
"""
function tt(a)
    return 1
end


"""
rel_error(u,v)

# Description

- calculate the relative error between u and v with respect to v

# Arguments

- `u` : a structure of numbers
- `v` : a structure of numbers

# return

- `relative error`:
"""
function rel_error(u,v)
    return maximum(abs.(u[:] .- v[:])) / maximum(abs.(u[:]))
end
