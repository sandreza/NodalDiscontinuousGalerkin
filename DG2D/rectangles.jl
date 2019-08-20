include("element2D.jl")

"""
rectangle(k, EtoV, N, M, vmap)

# Description

    create a rectangular element

# Arguments

-   `k`: element number in global map
-   `vertices`: local vertex indices
-   `N`: polynomial order along first axis within element
-   `M`: polynomial order along second axis within element
-   `vmap`: physical coordinates of global vertices

# Return Values:

-   `rect`: a properly initiliazed Element2D object

"""
function rectangle(index, vertices, N, M, vmap)
    # GL points in each dimension
    a = jacobiGL(0, 0, N)
    b = jacobiGL(0, 0, M)

    n = length(a)
    m = length(b)

    # differentiation and lift matrices through tensor products
    D = dmatricesSQ(a, b)
    V = vandermondeSQ(a, b)
    M = inv(V * V')

    # get min and max values of physical coordinates
    xmin = vmap[vertices[1], 1]
    ymin = vmap[vertices[1], 2]
    xmax = vmap[vertices[3], 1]
    ymax = vmap[vertices[3], 2]

    # arrays of first,second coordinate of GL tensor product
    # both ideal and physical coordinates are saved
    r̃ = zeros(n*m, 2)
    x̃ = similar(r̃)
    let k = 0
        for j in 1:m
            for i in 1:n
                k += 1

                r̃[k,:] = [a[i] b[j]]

                x = (xmax - xmin) * (a[i] + 1) / 2 + xmin
                y = (ymax - ymin) * (b[j] + 1) / 2 + ymin
                x̃[k,:] = [x y]
            end
        end
    end

    # build face mask for element
    fmask = fmaskSQ(r̃[:,1], r̃[:,2])

    # lift matrix and normals
    ∮ = liftSQ(a, b, fmask)
    nˣ,nʸ,Jˢ = normalsSQ(n, m, xmax-xmin, ymax-ymin)

    # construct element
    rect = Element2D(index,vertices, x̃,D,M, fmask,nˣ,nʸ,Jˢ,∮)

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
function liftSQ(r, s, fmask)
    # number of GL points in each dimension
    n = length(r)
    m = length(s)

    # empty matrix
    ∮ = spzeros(n*m, 2*(n+m))

    # get 1D mass matrices matrices,
    # need the minus 1 to get it to be the correct size
    Vʳ = vandermonde(r, 0, 0, n-1)
    Vˢ = vandermonde(s, 0, 0, m-1)

    Mʳ = inv(Vʳ * Vʳ')
    Mˢ = inv(Vˢ * Vˢ')

    # ending index for each face
    nf1 = n
    nf2 = n+m
    nf3 = n+m+n
    nf4 = n+m+n+m

    # fill ∮ matrix with mass matrices
    @. ∮[fmask[1],     1:nf1] = Mʳ
    @. ∮[fmask[2], nf1+1:nf2] = Mˢ
    @. ∮[fmask[3], nf2+1:nf3] = Mʳ
    @. ∮[fmask[4], nf3+1:nf4] = Mˢ

    return ∮
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

function normalsSQ(n, m, xwidth, ywidth)
    # empty vectors of right length
    nˣ = zeros(n+m+n+m)
    nʸ = zeros(n+m+n+m)
    Jˢ = ones(n+m+n+m) # squares don't need to worry about this

    # ending index for each face
    nf1 = n
    nf2 = n+m
    nf3 = n+m+n
    nf4 = n+m+n+m

    # set values of normals
    @. nʸ[    1:nf1] = -1 * xwidth/2 # normal is ( 0, -1) along first face
    @. nˣ[nf1+1:nf2] =  1 * ywidth/2 # normal is (-1,  0) along second face
    @. nʸ[nf2+1:nf3] =  1 * xwidth/2 # normal is ( 0,  1) along third face
    @. nˣ[nf3+1:nf4] = -1 * ywidth/2 # normal is ( 1,  0) along fourth face

    @. Jˢ = sqrt(nˣ^2 + nʸ^2)
    @. nˣ *= 1/Jˢ
    @. nʸ *= 1/Jˢ

    return nˣ,nʸ,Jˢ
end

"""
normalsRE(n, m)

# Description

    Return the normals for the 2D ideal square

# Arguments

-   `xmin`: minimum of x on rectangle
-   `xmax`: maximum of x on rectangle
-   `ymin`: minimum of y on rectangle
-   `ymax`: maximum of y on rectangle

-   `n`: number of GL points along the first axis
-   `m`: number of GL points along the second axis

# Return Values

-   `nˣ`: first coordinate of the normal vector
-   `nʸ`: second coordinate of the normal vector

# Example

"""

function normalsRE(n, m, xmin, xmax, ymin, ymax)
    # empty vectors of right length
    n̂  = zeros(n+m+n+m, 2)
    Jˢ = ones(n+m+n+m) # yes they do
    @. Jˢ[1:n] = xmax - xmin
    @. Jˢ[n+1:n+m] = ymax - ymin
    @. Jˢ[n+m+1:n+m+n] = xmax - xmin
    @. Jˢ[n+m+n+1:n+m+n+m] = ymax - ymin
    # ending index for each face
    nf1 = n
    nf2 = n+m
    nf3 = n+m+n
    nf4 = n+m+n+m

    # set values of normals
    @. n̂[    1:nf1, 2] = -1 # normal is ( 0, -1) along first face
    @. n̂[nf1+1:nf2, 1] =  1 # normal is (-1,  0) along second face
    @. n̂[nf2+1:nf3, 2] =  1 # normal is ( 0,  1) along third face
    @. n̂[nf3+1:nf4, 1] = -1 # normal is ( 1,  0) along fourth face

    return n̂, Jˢ
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
    fmask1 = Array(findall( abs.( s .+ 1) .< eps(1.0) ))'
    fmask2 = Array(findall( abs.( r .- 1) .< eps(1.0) ))'
    fmask3 = Array(reverse(findall( abs.( s .- 1) .< eps(1.0) )))'
    fmask4 = Array(reverse(findall( abs.( r .+ 1) .< eps(1.0) )))'
    fmasks = [fmask1', fmask2', fmask3', fmask4']
    return fmasks
end
