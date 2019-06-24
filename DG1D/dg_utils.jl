
using Revise
using SpecialFunctions #for gamma function reasons
using LinearAlgebra    #for Guass quadrature

#code checked against the matlab code
"""
jacobi(x, α, β, n)

evaluate a jacobi polynomial of degree n at the value x ∈ [-1, 1]
From Nodal Discontinuous Galerkin Methods by Hesthaven and Warburton
uses the gamma function
returns an array of size n+1, all the jacobi polynomials
up to value n
"""
function jacobi(x, α, β, n::Int)
    γ0 = 2^(α + β + 1) / (α + β + 1) * factorial(α) * factorial(β)
    γ0 /= factorial(α + β)
    γ1 = (α + 1) * (β + 1) / (α + β + 3) * γ0
    #create array to return
    PL = zeros(n+1)
    PL[1] = 1 / sqrt(γ0)
    if n==0
        return PL
    elseif n==1
        PL[2] = ( (α + β + 2) * x / 2 + (α - β) / 2) / sqrt(γ1)
        return PL
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
        return PL
    end
end

"""
jacobi!(y, x, α, β)

evaluate a jacobi polynomial of degree n at the values x ∈ [-1, 1]
From Nodal Discontinuous Galerkin Methods by Hesthaven and Warburton
uses the gamma function
returns an array of size n+1, all the jacobi polynomials
up to value n
overwrites the matrix y, an (length(x)) by (m=(n+1)) matrix
note that y is the "vandermonde matrix"

Example:
x = collect(-1:0.01:1);
jacobi!(y,x,0,0)

This gives all the Legendre polynomials up to 9 on the interval [-1,1]
i.e. y[:,3] is the Legendre polynomial of degree 2 evaluated at the
x values
plot(x,y[:,8])
Allocates a little bit of memory
"""
function jacobi!(y, x, α, β)
    γ0 = 2^(α + β + 1) / (α + β + 1) * factorial(α) * factorial(β)
    γ0 /= factorial(α + β)
    γ1 = (α + 1) * (β + 1) / (α + β + 3) * γ0
    #create view to assigne values
    m, n = size(y)
    if m != length(x)
        println("Make sure the size of the arrays is correct")
        return error
    end
    yview1 = view(y, :, 1)
    yview1 .= 1 / sqrt(γ0)
    if n == 2
        yview1 = view(y, :, 2)
        @. yview1 = ( (α + β + 2) * x / 2 + (α - β) / 2) / sqrt(γ1)
        return nothing
    elseif n > 2
        yview2 = view(y, :, 2)
        @. yview2 = ( (α + β + 2) * x / 2 + (α - β) / 2) / sqrt(γ1)
        aold = 2 / (2 + α + β) * sqrt((α+1)*(β+1)/(α + β + 3))

        for i in 1:(n-2)
            yview1 = view(y, :, i)
            yview2 = view(y, :, i+1)
            yview3 = view(y, :, i+2)
            h1 = 2 * i + α + β
            anew = 2 /(h1 + 2)*sqrt((i+1)*(i+1+α+β)*(i+1+α)*(i+1+β)/(h1+1)/(h1+3))
            bnew = - (α^2 - β^2)/h1/(h1+2)
            @. yview3 = 1 / anew * (-aold*yview1 + (x-bnew)*yview2)
            aold = anew
        end
        return nothing
    end
end

"""
djacobi!(y, x, α, β)

evaluate the derivative jacobi polynomial (α, β) of degree n at the values x ∈ [-1, 1].
From Nodal Discontinuous Galerkin Methods by Hesthaven and Warburton
returns an array of size n+1, the derivative of all the jacobi polynomials
up to of value n
overwrites the matrix y, an (length(x)) by (m=(n+1)) matrix

Example:

x = collect(-1:0.01:1);

x = collect(-1:0.5:1)

dy = ones(length(x),10);

djacobi!(dy,x,0,0)

This gives the derivative of Legendre polynomials up to 9 on the interval [-1,1]
i.e. dy[:,3] is the derivative of the Legendre polynomial of degree 2 evaluated at the
x values

plot(x,dy[:,8])

Allocates a little bit of memory
"""
function djacobi!(y, x, α, β)
    yview = view(y, :, 1)
    @. yview = 0
    m, n = size(y)
    if n>1
        yview = view(y, :, 2:n)
        jacobi!(yview, x, α+1, β+1)
        for i in 2:n
            yview2 = view(y, :, i)
            @. yview2 *= sqrt( (i-1) * (α + β + i) )
        end
    end
    return nothing
end

"""
dmatrix(x, α, β)

output is the differentiation matrix for nodal values of jacobi polynomial (α, β) of degree n at the values x ∈ [-1, 1]
From Nodal Discontinuous Galerkin Methods by Hesthaven and Warburton
returns a matrix of size length(x) by length(x)

Example:
x = collect(-1:0.5:1);
d = dmatrix(x,0,0)
y = x.^2 .+ 2 .* x .+ 8
dy  = 2 .* x .+ 2
dy .- d*y

allocates too much memory
"""
function dmatrix(x, α, β)
    n = length(x)
    vr = ones(n,n)
    v = ones(n,n)
    d = ones(n,n)
    djacobi!(vr, x, α, β)
    jacobi!(v, x, α, β)
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
#Description
    Guass Quadrature points and weights for the Jacobi Polynomial (α,β)
#Input
α, β: Jacobi polynomial descriptors
N:    order of quadrature points
#Return: x,w
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
    if N==0
        x[1] = (α - β) / (α + β + 2)
        w[1] = 2
    end
    h1 = @. 2 * collect(0:N) + α + β;
    diag = @. - (α^2 - β^2) / (h1 + 2) / h1
    h1view = view(h1, 1:N)
    cf = collect(1:N) #common factor that shows up a lot
    superdiag = @. 2 / (h1view + 2)
    @. superdiag *= sqrt( cf * (cf + α + β) * (cf + α) * (cf + β) )
    @. superdiag *= sqrt( 1 / (h1view + 1) / (h1view + 3) )
    J = SymTridiagonal(diag, superdiag)
    if (α + β) ≈ 0.0
        J[1,1] = -0.0
    end
    x , V = eigen(J)
    w = @. (V[1,:] ^ 2 ) * 2^(α + β + 1) / (α + β + 1)
    @. w *= factorial(α) * factorial(β) / factorial(α + β)
    return x, w
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
    x[1] = -1.0
    x[end] =  1.0
    if N==0
        error("What are you doing?")
    end
    if N==1
        return x
    end
    xview = view(x, 2:N)
    xtmp, w = jacobiGQ(α+1, β+1, N-2)
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
