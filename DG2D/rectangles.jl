include("element2D.jl")

"""
rectangle(k, EtoV, N, M, vmap)

# Description

    create a rectangular element

# Arguments

-   `k`: element number in global map
-   `EtoV`: element to vertex map
-   `N`: polynomial order along first axis within element
-   `M`: polynomial order along second axis within element
-   `vmap`: array of vertices

# Return Values:

-   `rect`: a properly initiliazed Element2D object

"""
function rectangle(index, EtoV, N, M, vmap)
    vertices = view(EtoV, index, :)
    nFaces = length(vertices)

    # GL points in each dimension
    a = jacobiGL(0, 0, N)
    b = jacobiGL(0, 0, M)

    n = length(a)
    m = length(b)

    # get normals
    n̂,Jˢ = normalsSQ(n, m)

    # differentiation and lift matrices through tensor products
    D = dmatricesSQ(a, b)
    lift = liftSQ(a, b)

    # get min and max values of physical coordinates
    xmin = vmap[vertices[2]][1]
    ymin = vmap[vertices[2]][2]
    xmax = vmap[vertices[end]][1]
    ymax = vmap[vertices[end]][2]

    # arrays of first,second coordinate of GL tensor product
    # both ideal and physical coordinates are saved
    r̃ = zeros(n*m, 2)
    x̃ = similar(r̃)
    let k = 0
        for i in 1:n
            for j in 1:m
                k += 1

                r̃[k, :] = [a[i] b[j]]

                x = (xmax - xmin) * (a[i] + 1) / 2 + xmin
                y = (ymax - ymin) * (b[j] + 1) / 2 + ymin
                x̃[k, :] = [x y]
            end
        end
    end

    # build face mask for element
    fmask = fmaskSQ(r̃[:,1], r̃[:,2])

    # construct element
    rect = Element2D(index,vertices, r̃,x̃, fmask,n̂,Jˢ, D,lift)

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
function liftSQ(r, s)
    # get 2D vandermonde matrix
    V = vandermondeSQ(r,s)

    # number of GL points in each dimension
    n = length(r)
    m = length(s)

    # empty matrix
    ℰ = spzeros(Int, n*m, 2*(n+m))

    # starting column number for each face
    rl = 1
    sl = 1+m
    rh = 1+m+n
    sh = 1+m+n+m

    # fill matrix for bounds on boundaries
    # += syntax used for debugging, easily shows if multiple statements assign to the same entry
    let k = 0 # element number
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
    end

    # compute lift (ordering because ℰ is sparse)
    lift = V * (V' * ℰ)

    return lift
end

"""
normalsSQ(n, m)

# Description

    Return the normals for the 2D ideal square

# Arguments

-   `n`: number of GL points along the first axis
-   `m`: number of GL points along the second axis

# Return Values

-   `nˣ`: first coordinate of the normal vector
-   `nʸ`: second coordinate of the normal vector

# Example

"""

function normalsSQ(n, m)
    # empty vectors of right length
    n̂  = zeros(n+m+n+m, 2)
    Jˢ = ones(n+m+n+m) # squares don't need to worry about this

    # ending index for each face
    nf1 = m
    nf2 = m+n
    nf3 = m+n+m
    nf4 = m+n+m+n

    # set values of normals
    for i in 1:nf4
        if     i < nf1+1
            n̂[i, :] = [ 0 -1] # normal is (0, -1) along first face
        elseif i < nf2+1
            n̂[i, :] = [-1  0] # normal is (-1, 0) along second face
        elseif i < nf3+1
            n̂[i, :] = [ 0  1] # normal is (0, 1) along third face
        else
            n̂[i, :] = [ 1  0] # normal is (1, 0) along fourth face
        end
    end

    return n̂,Jˢ
end

"""
fmaskSQ(r, s)

# Description

-   Mask of GL nodes on the faces of a square element

# Arguments

-   `r` : first index of real coordinates
-   `s` : second index of real coordinates


#Return

-   `fmask`: matrix of indices of GL points along each face

"""
function fmaskSQ(r, s)
    fmask1 = findall( abs.( r .+ 1) .< eps(1.0) )'
    fmask2 = findall( abs.( s .+ 1) .< eps(1.0) )'
    fmask3 = findall( abs.( r .- 1) .< eps(1.0) )'
    fmask4 = findall( abs.( s .- 1) .< eps(1.0) )'
    fmask = Array([fmask1; fmask2; fmask3; fmask4]')
    return fmask
end
