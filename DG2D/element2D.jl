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
struct Element2D{S, T, U, V, W, X} <: AbstractElement2D
    # identifying features
    index::S
    vertices::T

    # GL points and normals
    r::U # ideal coordinates
    x::U # physical coordinates
    n̂::U # normal vectors

    # matrices for computation
    D::V
    lift::U

    # geometric factors
    J::W
    xʳ::X
    rˣ::X

    function Element2D(index,vertices, r̃,x̃, D,lift,n̂)
        # partial derivatives of x
        nGL,nDim = size(x̃)
        x̃ʳ = zeros(nGL, 2, 2)
        r̃ˣ = similar(x̃ʳ)
        J = zeros(nGL)

        # compute the derivates component wise
        xʳ = D[1] * x̃[:, 1]
        xˢ = D[2] * x̃[:, 1]
        yʳ = D[1] * x̃[:, 2]
        yˢ = D[2] * x̃[:, 2]

        # save partials as jacobian matrix, inverse, and determinant
        for i in 1:nGL
            𝒥 = [ [xʳ[i] xˢ[i]]; [yʳ[i] yˢ[i]]]
            x̃ʳ[i, :, :] = 𝒥
            r̃ˣ[i, :, :] = inv(𝒥)
            J[i] = det(𝒥)
        end

        return new{typeof(index),typeof(vertices),typeof(r̃),typeof(D),typeof(J),typeof(x̃ʳ)}(index,vertices, r̃,x̃,n̂, D,lift, J,x̃ʳ,r̃ˣ)
    end
end

### exampleeeee
# function nFaces(::Element2D{N}) where N
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
