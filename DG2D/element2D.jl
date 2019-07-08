include("../utils.jl")

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
mutable struct element2D <: AbstractElement2D
    # identifying features
    index::Int
    vertices::Array{Int,1}

    # ideal coordinates
    r::Array{Float64,1}
    s::Array{Float64,1}

    # physical coordinates
    x::Array{Float64,1}
    y::Array{Float64,1}

    # matrices for computation
    Dʳ::Array{Float64,2}
    Dˢ::Array{Float64,2}
    lift::Array{Float64,2}

    # geometric factors
    J::Array{Float64,1}
    xʳ::Array{Float64,1}
    xˢ::Array{Float64,1}
    yʳ::Array{Float64,1}
    yˢ::Array{Float64,1}
    rˣ::Array{Float64,1}
    rʸ::Array{Float64,1}
    sˣ::Array{Float64,1}
    sʸ::Array{Float64,1}


    # only index and vertices are well defined for an element of arbitrary shape
    function element2D(k, EtoV)
        element = new()
        element.index = k
        element.vertices = view(EtoV, k, :)

        return element
    end
end

"""
geometricfactors2D!(Ω)

# NEEDS TO BE TESTED

# Description

-   Compute metric factors for an element

# Arguments

-   `Ω`: element to compute metric factors for

"""
function geometricfactors2D!(Ω)
    # partial derivatives of x
    Ω.xʳ = Ω.Dʳ * Ω.x
    Ω.xˢ = Ω.Dˢ * Ω.x

    # partial derivatives of y
    Ω.yʳ = Ω.Dʳ * Ω.y
    Ω.yˢ = Ω.Dˢ * Ω.y

    # Jacobian
    Ω.J = - Ω.xˢ .* Ω.yʳ + Ω.xʳ .* Ω.yˢ #determinant

    # partial derivatives of r
    Ω.rˣ =   Ω.yˢ ./ Ω.J
    Ω.rʸ = - Ω.xˢ ./ Ω.J

    # partial derivatives of s
    Ω.sˣ = - Ω.yʳ ./ Ω.J
    Ω.sʸ =   Ω.xʳ ./ Ω.J
    return
end

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
