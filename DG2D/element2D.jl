include("../src/utils.jl")

abstract type AbstractFace end
abstract type AbstractFace2D <: AbstractFace end
"""
Face()

# Description

    initialize a face struct

#

"""
struct Face2D{S, T, U, V, W} <: AbstractFace2D
    # identifying features
    index::S
    mask::T # indices of local GL points

    # number of GL points
    nGL::S

    # indices of global GL points
    i⁻::T # interior
    i⁺::T # exterior
    isBoundary::U

    # normals and lift operator for this face
    nˣ::V
    nʸ::V
    C::V  # compactness, or surface-area-to-volume ratio
    ∮::W

    function Face2D(index, mask, C,nˣ,nʸ,∮)
        nGL = length(mask)
        isBoundary = [false]

        # default assignment
        i⁻ = similar(mask)
        i⁺ = similar(mask)

        return new{typeof(index),typeof(mask),typeof(isBoundary),typeof(nˣ),typeof(∮)}(index,mask, nGL,i⁻,i⁺,isBoundary, nˣ,nʸ,C,∮)
    end
end

abstract type AbstractElement end
abstract type AbstractElement2D <: AbstractElement end
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
struct Element2D{S, T, U, V, W, X, Y, Z} <: AbstractElement2D
    # identifying features
    index::S
    vertices::T

    # volume information
    nGL::S # number of points
    x::U   # physical coordinates
    iⱽ::V  # global indices of GL points

    # boundary information
    faces::W  # Array of Face structs

    # geometric factors
    J::X   # magnitude of the jacobian
    rˣ::Y  # jacobian matrix from ideal to physical space
    D::Z   # differentiation matrices
    M::U   # mass matrix
    M⁺::U  # inverse of mass matrix

    function Element2D(index,vertices, x̃,D,M, fmasks,nˣ,nʸ,Jˢ,∮)
        # indices of GL points
        nGL,nDim = size(x̃)
        iⱽ = collect(Int, 1:nGL)

        # partial derivatives of x
        x̃ʳ = zeros(nGL, nDim, nDim)
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
            J[i] = det(𝒥)
        end

        # inverse of mass Matrix
        M⁺ = inv(M)

        # construct faces
        nBP = 0
        faces = Face2D[]
        for (f, fmask) in enumerate(fmasks)
            BPᶠ = (nBP + 1):(nBP + length(fmask))
            nBP += length(fmask)

            C = @. Jˢ[BPᶠ] / J[fmask]

            face = Face2D(f, fmask, C, nˣ[BPᶠ], nʸ[BPᶠ], ∮[:, BPᶠ])

            push!(faces, face)
        end

        return new{typeof(index),typeof(vertices),typeof(x̃),typeof(iⱽ),typeof(faces),typeof(J),typeof(r̃ˣ),typeof(D)}(index,vertices, nGL,x̃,iⱽ, faces, J,r̃ˣ,D,M,M⁺)
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
