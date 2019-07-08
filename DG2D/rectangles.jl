include("element2D.jl")

"""
rectangle(k, N, M, vmap, EtoV)

# Description

    create a rectangular element

# Arguments

-   `k`: element number in global map
-   `EtoV`: element to vertex map
-   `N`: polynomial order along first axis within element
-   `M`: polynomial order along second axis within element
-   `vmap`: array of vertices

# Return Values:

-   `rect`: a rectangular element object initialized with proper index, vertices, grid points, and geometric factors

"""
function rectangle(k, EtoV, N, M, vmap)
    vertices = view(EtoV, k, :)

    # GL points in each dimension
    a = jacobiGL(0, 0, N)
    b = jacobiGL(0, 0, M)

    # differentiation and lift matrices through tensor products
    Dʳ,Dˢ = dmatricesSQ(a, b)
    lift = liftSQ(a, b)

    # arrays of first,second coordinate of GL tensor product
    r = []
    s = []
    for i in a
        for j in b
            push!(r, i)
            push!(s, j)
        end
    end

    # get min and max values of physical coordinates
    xmin = vmap[vertices[2]][1]
    ymin = vmap[vertices[2]][2]
    xmax = vmap[vertices[end]][1]
    ymax = vmap[vertices[end]][2]

    # create physical coordinates of GL points
    x = @. (xmax - xmin) * (r + 1) / 2
    y = @. (ymax - ymin) * (s + 1) / 2

    # construct element
    rect = Element2D{4}(k,vertices, r,s, x,y, Dʳ,Dˢ,lift)

    return rect
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
function vandermondeSQ(r, s)
    # get order of GL points
    N = length(r) - 1
    M = length(s) - 1

    # construct 1D vandermonde matrices
    Vʳ = vandermonde(r, 0, 0, N)
    Vˢ = vandermonde(s, 0, 0, M)

    # construct identity matrices
    Iⁿ = Matrix(I, N+1, N+1)
    Iᵐ = Matrix(I, M+1, M+1)

    # compute 2D vandermonde matrix
    V = kron(Iᵐ, Vʳ) * kron(Vˢ, Iⁿ)
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

-   `Dʳ`: the differentiation matrix wrt first coordinate
-   `Dˢ`: the differentiation matrix wrt to second coordinate

# Example

"""
function dmatricesSQ(r, s)
    # get order of GL points
    N = length(r) - 1
    M = length(s) - 1

    # construct 1D vandermonde matrices
    Dʳ = dmatrix(r, 0, 0, N)
    Dˢ = dmatrix(s, 0, 0, M)

    # construct identity matrices
    Iⁿ = Matrix(I, N+1, N+1)
    Iᵐ = Matrix(I, M+1, M+1)

    # compute 2D vandermonde matrix
    Dʳ = kron(Iᵐ, Dʳ)
    Dˢ = kron(Dˢ, Iⁿ)

    return Dʳ,Dˢ
end

"""
dvandermondeSQ(N, M)

# Description

    Return gradient matrices of the 2D vandermonde matrix evaluated at tensor product of NxM GL points on an ideal [-1,1]⨂[-1,1] square

# Arguments

-   `N`: polynomial order in first coordinate
-   `M`: polynomial order in second coordinate

# Return Values

-   `Vʳ`: gradient of vandermonde matrix wrt to first coordinate
-   `Vˢ`: gradient of vandermonde matrix wrt to second coordinate

# Example

"""
function dvandermondeSQ(r, s)
    # get 2D vandermonde matrix
    V = vandermondeSQ(r, s)

    # get differentiation matrices
    Dʳ,Dˢ = dmatricesSQ(r, s)

    # calculate using definitions
    Vʳ = Dʳ * V
    Vˢ = Dˢ * V

    return Vʳ,Vˢ
end

"""
liftSQ(N, M)

# Description

    Return the 2D lift matrix evaluated at tensor product of NxM GL points on an ideal [-1,1]⨂[-1,1] square

# Arguments

-   `N`: polynomial order in first coordinate
-   `M`: polynomial order in second coordinate

# Return Values

-   `lift`: the 2D lift matrix

# Example

"""
function liftSQ(r,s)
    # get 2D vandermonde matrix
    V = vandermondeSQ(r,s)

    # number of GL points in each dimension
    n = length(r)
    m = length(s)

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
    for i in 1:n
        for j in 1:m
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
            if i == n
                ℰ[k, rh] += 1
                rh += 1
            end

            # check if on smax
            if j == m
                ℰ[k, sh] += 1
                sh += 1
            end
        end
    end

    # compute lift (ordering because ℰ is sparse)
    lift = V * (V' * ℰ)

    return lift
end
