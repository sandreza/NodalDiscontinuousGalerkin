include("../src/utils.jl")

abstract type AbstractElement2D end
"""
Element2D(index,vertices, r̃,x̃,n̂, D,lift,fmask)

# Description

    initialize 2D element struct

# Arguments

-   `index`: element number in global map
-   `vertices`: view of vertices this element has
-   `r̃`: ideal coordinates of GL points
-   `x̃`: physical coordinates of GL points
-   `n̂`: normal vectors along the faces
-   `D`: tuple of derivative matrices
-   `lift`: lift matrix
-   `fmask`: matrix of indices of GL points along each face


# Return Values:

    return a properly initiliazed Element2D object

"""
struct Element2D{S, T, U, V, W, X, Y} <: AbstractElement2D
    # identifying features
    index::S
    vertices::T

    # GL points and normals
    x::U # physical coordinates
    nGL::S # number of points

    # boundary information
    fmask::V
    nBP::S # number of points on the boundary
    n̂::U # normal vectors

    # matrices for computation
    D::W
    lift::U

    # geometric factors
    J::X
    xʳ::Y
    rˣ::Y
    volume::X

    function Element2D(index,vertices, x̃, fmask,n̂,Jˢ, D,lift)
        # number of points on the boundary
        nFPᵏ,nFaces = size(fmask)
        nBP = nFPᵏ * nFaces

        # partial derivatives of x
        nGL,nDim = size(x̃)
        x̃ʳ = zeros(nGL, 2, 2)
        r̃ˣ = similar(x̃ʳ)
        J = zeros(nGL)

        # compute the derivates component wise
        xʳ = D[1] * x̃[:,1]
        xˢ = D[2] * x̃[:,1]
        yʳ = D[1] * x̃[:,2]
        yˢ = D[2] * x̃[:,2]

        # save partials as jacobian matrix, inverse, and determinant
        for i in 1:nGL
            𝒥 = [ [xʳ[i] xˢ[i]]; [yʳ[i] yˢ[i]]]
            x̃ʳ[i,:,:] = 𝒥
            r̃ˣ[i,:,:] = inv(𝒥)
            J[i] = -det(𝒥)
        end

        # volume of element
        volume = @. Jˢ / J[fmask][:]

        #### add nodes⁻ and nodes⁺ as struct members

        return new{typeof(index),typeof(vertices),typeof(x̃),typeof(fmask),typeof(D),typeof(J),typeof(x̃ʳ)}(index,vertices, x̃,nGL, fmask,nBP,n̂, D,lift, J,x̃ʳ,r̃ˣ,volume)
    end
end

### exampleeeee
# function nFaces(::Element2D{N}) where N
#     return N
# end




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
