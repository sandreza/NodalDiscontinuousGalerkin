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
vandermondeSQ(N, M)

# Description

    Return 2D vandermonde matrix using squares evaluated at tensor product of NxM GL points on an ideal [-1,1]⨂[-1,1] square

# Arguments

-   `N`: polynomial order in first coordinate
-   `M`: polynomial order in second coordinate

# Return Values

-   `V`: the 2D vandermonde matrix

# Example

"""
function vandermondeSQ(N, M)
    # get GL nodes in each dimnesion
    r = jacobiGL(0, 0, N)
    s = jacobiGL(0, 0, M)

    # construct 1D vandermonde matrices
    Vr = vandermonde(r, 0, 0, N)
    Vs = vandermonde(s, 0, 0, M)

    # construct identity matrices
    Iⁿ = Matrix(I, N+1, N+1)
    Iᵐ = Matrix(I, M+1, M+1)

    # compute 2D vandermonde matrix
    V = kron(Iᵐ, Vr) * kron(Vs, Iⁿ)
    return V
end

"""
dmatricesSQ(N, M)

# Description

    Return the 2D differentiation matrices evaluated at tensor product of NxM GL points on an ideal [-1,1]⨂[-1,1] square

# Arguments

-   `N`: polynomial order in first coordinate
-   `M`: polynomial order in second coordinate

# Return Values

-   `Dr`: the differentiation matrix wrt first coordinate
-   `Ds`: the differentiation matrix wrt to second coordinate

# Example

"""
function dmatricesSQ(N, M)
    # get GL nodes in each dimnesion
    r = jacobiGL(0, 0, N)
    s = jacobiGL(0, 0, M)

    # construct 1D vandermonde matrices
    Dr = dmatrix(r, 0, 0, N)
    Ds = dmatrix(s, 0, 0, M)

    # construct identity matrices
    Iⁿ = Matrix(I, N+1, N+1)
    Iᵐ = Matrix(I, M+1, M+1)

    # compute 2D vandermonde matrix
    Dr = kron(Iᵐ, Dr)
    Ds = kron(Ds, Iⁿ)

    return Dr,Ds
end

"""
dvandermondeSQ(N, M)

# Description

    Return gradient matrices of the 2D vandermonde matrix evaluated at tensor product of NxM GL points on an ideal [-1,1]⨂[-1,1] square

# Arguments

-   `N`: polynomial order in first coordinate
-   `M`: polynomial order in second coordinate

# Return Values

-   `Vr`: gradient of vandermonde matrix wrt to first coordinate
-   `Vs`: gradient of vandermonde matrix wrt to second coordinate

# Example

"""
function dvandermondeSQ(N, M)
    # get 2D vandermonde matrix
    V = vandermondeSQ(N, M)

    # get differentiation matrices
    Dr,Ds = dmatricesSQ(N, M)

    # calculate using definitions
    Vr = Dr * V
    Vs = Ds * V

    return Vr,Vs
end

"""
LiftSQ(N, M)

# Description

    Return the 2D lift matrix evaluated at tensor product of NxM GL points on an ideal [-1,1]⨂[-1,1] square

# Arguments

-   `N`: polynomial order in first coordinate
-   `M`: polynomial order in second coordinate

# Return Values

-   `lift`: the 2D lift matrix

# Example

"""
function LiftSQ(N,M)
    # get 2D vandermonde matrix
    V = vandermondeSQ(N,M)

    # number of GL points in each dimension
    n = N+1
    m = M+1

    # empty matrix
    ℰ = spzeros(n*m, 2*(n+m))

    # starting column number for each face
    rl = 1
    sl = 1+m
    rh = 1+m+n
    sh = 1+m+n+m

    # element number
    k = 0

    # fill matrix for bounds on boundaries
    # += syntax used for debugging, easily shows if multiple statements assign to the same entry
    for i in 1:N+1
        for j in 1:M+1
            k += 1

            # check if on rmin
            if i == 1
                ℰ[k, rl] += 1
                rl += 1
            end

            # check if on smax
            if j == 1
                ℰ[k, sl] += 1
                sl += 1
            end

            # check if on rmax
            if i == N+1
                ℰ[k, rh] += 1
                rh += 1
            end

            # check if on smax
            if j == M+1
                ℰ[k, sh] += 1
                sh += 1
            end
        end
    end

    # compute lift (ordering because ℰ is sparse)
    lift = V * (V' * ℰ)

    return lift
end
