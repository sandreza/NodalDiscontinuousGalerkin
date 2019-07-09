include("../src/utils.jl")

abstract type AbstractElement2D end
"""
element2D(k, N, M, vmap, EtoV)

# Description

    initialize 2D element struct

# Arguments

-   `k`: element number in global map
-   `EtoV`: element to vertex map

# Return Values: x

    return index and vertices

"""
struct Element2D{N, S, T, U, V} <: AbstractElement2D
    # identifying features
    index::S
    vertices::T

    # ideal coordinates
    r::U
    s::U

    # physical coordinates
    x::U
    y::U

    # matrices for computation
    Dʳ::V
    Dˢ::V
    lift::V
    nˣ::V
    nʸ::V

    # geometric factors
    J::U
    xʳ::U
    xˢ::U
    yʳ::U
    yˢ::U
    rˣ::U
    rʸ::U
    sˣ::U
    sʸ::U

    function Element2D{N}(index,vertices, r,s, x,y, Dʳ,Dˢ,lift, nˣ,nʸ) where N
        xʳ = Dʳ * x
        xˢ = Dˢ * x

        # partial derivatives of y
        yʳ = Dʳ * y
        yˢ = Dˢ * y

        # Jacobian
        J =  @. - xˢ * yʳ + xʳ * yˢ # determinant

        # partial derivatives of r
        rˣ = @.   yˢ / J
        rʸ = @. - xˢ / J

        # partial derivatives of s
        sˣ = @. - yʳ / J
        sʸ = @.   xʳ / J

        return new{N,typeof(index),typeof(vertices),typeof(r),typeof(lift)}(index,vertices, r,s, x,y, Dʳ,Dˢ,lift, nˣ,nʸ, J, xʳ,xˢ,yʳ,yˢ, rˣ,rʸ,sˣ,sʸ)
    end
end

### exampleeeee
# function nfaces(::Element2D{N}) where N
#     return N
# end

"""
∇(u, Ω)

# Description

    Compute gradient of u wrt physical grid

# Arguments

-   `u`: scalar to take gradient of
-   `Ω`: element to compute in

# Return Values

-   `uˣ`: first component of the gradient
-   `uʸ`: second component of the gradient

"""
function ∇(u, Ω)
    # compute partial derivatives on ideal grid
    uʳ = Ω.Dʳ * u
    uˢ = Ω.Dˢ * u

    # compute partial derivatives on physical grid
    uˣ = @. Ω.rˣ * uʳ + Ω.sˣ * uˢ
    uʸ = @. Ω.rʸ * uʳ + Ω.sʸ * uˢ

    return uˣ,uʸ
end

"""
∇⨀(x, y, Ω)

# Description

    Compute the divergence of u=(x,y) wrt physical grid

# Arguments

-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ω`: element to compute in

# Return Values

-   `∇⨀u`: the divergence of u

"""
function ∇⨀(x, y, Ω)
    # compute partial derivatives on ideal grid
    xʳ = Ω.Dʳ * x
    xˢ = Ω.Dˢ * x
    yʳ = Ω.Dʳ * y
    yˢ = Ω.Dˢ * y

    # compute gradient on physical grid
    ∇⨀u = @. Ω.rˣ * xʳ + Ω.sˣ * xˢ + Ω.rʸ * yʳ + Ω.sʸ * yˢ

    return ∇⨀u
end

"""
∇⨂(x, y, Ω)

# Description

    Compute the curl of u=(x,y) wrt physical grid

# Arguments

-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ω`: element to compute in

# Return Values

-   `∇⨂u`: the curl of u

"""
function ∇⨂(x, y, Ω)
    # compute partial derivatives on ideal grid
    xʳ = Ω.Dʳ * x
    xˢ = Ω.Dˢ * x
    yʳ = Ω.Dʳ * y
    yˢ = Ω.Dˢ * y

    # compute gradient on physical grid
    ∇⨂u = @. Ω.rˣ * yʳ + Ω.sˣ * yˢ - Ω.rʸ * xʳ - Ω.sʸ * xˢ

    return ∇⨂u
end

"""
phys2ideal(x, y, Ω)

# Description

    Converts from physical rectangle Ω to ideal [-1,1]⨂[-1,1] square for legendre interpolation

# Arguments

-   `x`: first physical coordinate
-   `y`: second physical coordinate
-   `Ω`: element to compute in

# Return Values

-   `r`: first ideal coordinate
-   `s`: second ideal coordinate

# Example

"""
function phys2ideal(x, y, Ω)
    r = Ω.rˣ * (x - Ω.xmin) + Ω.rʸ * (y - Ω.ymin) - 1
    s = Ω.sˣ * (x - Ω.xmin) + Ω.sʸ * (y - Ω.ymin) - 1

    return r,s
end

"""
ideal2phys(r, s, Ω)

# Description

    Converts from ideal [-1,1]⨂[-1,1] square to physical rectangle Ω

# Arguments

-   `r`: first ideal coordinate
-   `s`: second ideal coordinate
-   `Ω`: element to compute in

# Return Values

-   `x`: first physical coordinate
-   `y`: second physical coordinate

# Example

"""
function ideal2phys(r, s, Ω)
    x = Ω.xʳ * (r + 1) + Ω.xˢ * (s + 1) + Ω.xmin
    y = Ω.yʳ * (r + 1) + Ω.yˢ * (s + 1) + Ω.ymin

    return x,y
end
