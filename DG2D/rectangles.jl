include("../utils.jl")

"""
rectmesh2D(xmin, xmax, ymin, ymax, K, L)

# Description

    Generates a 2D mesh of uniform squares

# Arguments

-   `xmin`: smallest value of first dimension
-   `xmax`:  largest value of first dimension
-   `ymin`: smallest value of second dimension
-   `ymax`:  largest value of second dimension
-   `K`: number of divisions in first dimension
-   `L`: number of divisions in second dimension

# Return Values: VX, EtoV

-   `VX`: vertex values | an Array of size K+1
-   `EtoV`: element to node connectivity | a Matrix of size Kx2

# Example

"""
function rectmesh2D(xmin, xmax, ymin, ymax, K, L)
    # 1D arrays
    vx,mapx = unimesh1D(xmin, xmax, K)
    vy,mapy = unimesh1D(ymin, ymax, L)

    # construct array of vertices
    vertices = Any[] # need to find a way to not make type Any
    for x in vx
        for y in vy
            push!(vertices, (x,y))
        end
    end
    # v = reshape(v, K+1, L+1)

    # construct element to vertex map
    EtoV = Int.(ones(K*L, 4))
    j = 1
    for l in 1:L
        for k in 1:K
            EtoV[j,1] = Int(k + (L+1) * (l-1))
            EtoV[j,3] = Int(k + (L+1) * l)

            EtoV[j,2] = Int(EtoV[j,1] + 1)
            EtoV[j,4] = Int(EtoV[j,3] + 1)

            j += 1
        end
    end

    return vertices,EtoV
end

"""
rectangle(k, vmap, EtoV)

# Description

    initialize rectangle struct

# Arguments

-   `k`: element number in global map
-   `vmap`: array of vertices
-   `EtoV`: element to vertex map

# Return Values: x

    return index and vertices

"""
struct rectangle{T, S, U}
    index::T
    vertices::S

    xmin::U
    xmax::U
    ymin::U
    ymax::U

    function rectangle(k, vmap, EtoV)
        index = k
        vertices = view(EtoV, k, :)

        xmin = vmap[vertices[1]][1]
        ymin = vmap[vertices[1]][2]
        xmax = vmap[vertices[end]][1]
        ymax = vmap[vertices[end]][2]

        return new{typeof(index),typeof(vertices),typeof(xmin)}(index, vertices, xmin,xmax,ymin,ymax)
    end
end

"""
phys2ideal(x, y, Dᵏ)

# Description

    Converts from physical rectangle Dᵏ to ideal [-1,1]⨂[-1,1] square for legendre interpolation

# Arguments

-   `x`: first physical coordinate
-   `y`: second physical coordinate
-   `Dᵏ`: element to compute in

# Return Values

-   `r`: first ideal coordinate
-   `s`: second ideal coordinate

# Example

"""
function phys2ideal(x, y, Dᵏ)
    r = 2 * (x - Dᵏ.xmin) / (Dᵏ.xmax - Dᵏ.xmin) - 1
    s = 2 * (y - Dᵏ.ymin) / (Dᵏ.ymax - Dᵏ.ymin) - 1

    return r,s
end

"""
ideal2phys(r, s, Dᵏ)

# Description

    Converts from ideal [-1,1]⨂[-1,1] square to physical rectangle Dᵏ

# Arguments

-   `r`: first ideal coordinate
-   `s`: second ideal coordinate
-   `Dᵏ`: element to compute in

# Return Values

-   `x`: first physical coordinate
-   `y`: second physical coordinate

# Example

"""
function ideal2phys(r, s, Dᵏ)
    x = 1//2 * (r + 1) * (Dᵏ.xmax - Dᵏ.xmin) + Dᵏ.xmin
    y = 1//2 * (s + 1) * (Dᵏ.ymax - Dᵏ.ymin) + Dᵏ.ymin

    return x,y
end

"""
vandermonde_square(r, s)

# Description

    Return 2D vandermonde matrix using squares

# Arguments

-   `r`: first coordinate
-   `s`: second coordinate

# Return Values

-   `V`: vandermonde matrix evaluated at (r,s) using squares

# Example

"""
function vandermonde_square(α, β, N, M)
    # get GL nodes in each dimnesion
    r = jacobiGL(α, β, N)
    s = jacobiGL(α, β, M)

    # construct 1D vandermonde matrices
    Vⁿ = vandermonde(r, 0, 0, N)
    Vᵐ = vandermonde(s, 0, 0, M)

    # construct identity matrices
    Iⁿ = Matrix(I, N, N)
    Iᵐ = Matrix(I, M, M)

    # compute 2D vandermonde matrix
    V = kron(Iᵐ, Vⁿ) * kron(Vᵐ, Iⁿ)'
    return V
end
