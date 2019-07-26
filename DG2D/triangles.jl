include("../src/utils.jl")
include("element2D.jl")

using SparseArrays

"""
triangle(k, EtoV, N, vmap)

# Description

    create a triangular element

# Arguments

-   `k`: element number in global map
-   `vertices`: local vertex indices
-   `N`: polynomial order along each axis within element
-   `vmap`: physical coordinates of global vertices

# Return Values:

-   `tri`: a properly initiliazed Element2D object

"""
function triangle(index, vertices, N, vmap)
    # GL points
    a,b = nodes2D(N)
    r,s = xytors(a,b)

    # build reference elements and matrices
    V = vandermonde2D(N, r, s)
    D = dmatrices2D(N, r, s, V)
    M = inv(V * V')

    # create fmask
    fmask = create_fmask(r, s)
    ∮     = lift_tri(N, fmask, r, s, V)

    # get physical vertices
    x1 = vmap[vertices[1], 1]
    y1 = vmap[vertices[1], 2]
    x2 = vmap[vertices[2], 1]
    y2 = vmap[vertices[2], 2]
    x3 = vmap[vertices[3], 1]
    y3 = vmap[vertices[3], 2]

    # create physical GL points
    x̃ = zeros(length(r), 2)
    for (k, (r̃,s̃)) in enumerate(zip(r,s))
        x =  0.5 * ( -(r̃ + s̃) * x1 + (1 + r̃) * x2 + (1 + s̃) * x3)
        y =  0.5 * ( -(r̃ + s̃) * y1 + (1 + r̃) * y2 + (1 + s̃) * y3)

        x̃[k,:] = [x y]
    end

    # construct normals
    nx, ny, Jˢ = normals2D(x̃[:,1], x̃[:,2], D[1], D[2], fmask, N+1, 1)
    Jˢ = reshape(Jˢ, length(Jˢ))

    n̂ = zeros(length(nx), 2)
    for (i, (x,y)) in enumerate(zip(nx,ny))
        n̂[i,:] = [x y]
    end

    # construct element
    tri = Element2D(index,vertices, x̃,D,M, fmask,n̂,Jˢ,∮)

    return tri
end

"""
rstoab(r,s)

# Description

Changes from (r,s) coordinate to (a,b) coordinates.

# Inputs

- `r`: a vector
- `s`: a vector

# Return

- `a`: a vector: first coordinate pair
- `b`: a vector; second coordinate pair

# Example

a, b = rstoab([-0.5, 0.5], [0.5, -0.5])

"""
function rstoab(r,s)
    nGL = length(r);
    a = zeros(nGL)
    for i in 1:nGL
        if (s[i] ≈ 1)
            a[i] = -1
        else
            a[i] = 2*(1+r[i]) / (1-s[i]) -1
        end
    end
    b = s
    return a, b
end

function warp_factor(n, rout)
    LGLr = jacobiGL(0, 0, n)
    req = collect(range(-1,1, length = n+1))
    veq = vandermonde(req, 0, 0, n)
    nr = length(rout)
    pmat = vandermonde(rout, 0, 0, n)'
    lmat = veq' \ pmat
    warp = lmat' * (LGLr - req)
    zerof = (abs.(rout) .< (1.0 -  1.0e-10) )
    sf = @. 1.0 - (zerof *rout )^2
    warp = warp ./ sf + warp .* (zerof .- 1)
    return warp
end

"""
nodes2D(n)

# Description

- Returns the interpolation nodes for the 2D GL points

# Return: x, y

- `x`: x-coordinate of GL points
- `y`: y-coordinate of GL points

"""
function nodes2D(n)
    # constants
    α° = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258]

    if n < 16
        α = α°[n]
    else
        α = 5/3
    end
    # number of GL points
    nGL = Int((n+1) * (n+2) / 2)

    # create equidistributed nodes on equilateral triangle
    L1 = zeros(nGL)
    L2 = zeros(nGL)
    L3 = zeros(nGL)
    let sk = 1
        for j ∈ 1:(n+1)
            for m ∈ 1:(n+2-j)
                L1[sk] = (j-1)/n
                L3[sk] = (m-1)/n
                sk += 1
            end
        end
    end
    @. L2 = 1.0 - L1 - L3
    x = -L2 + L3
    y = (-L2 - L3 + 2 .* L1 ) / sqrt(3)

    # compute blending function at each node for each edge
    blend1 = 4 * L2 .* L3
    blend2 = 4 * L3 .* L1
    blend3 = 4 * L1 .* L2

    # Amount of warp for each node for each edge
    warpf1 = warp_factor(n, L3 - L2)
    warpf2 = warp_factor(n, L1 - L3)
    warpf3 = warp_factor(n, L2 - L1)

    # combine blend and warp
    warp1 = @. blend1 * warpf1 * (1 + (α * L1 )^2 )
    warp2 = @. blend2 * warpf2 * (1 + (α * L2 )^2 )
    warp3 = @. blend3 * warpf3 * (1 + (α * L3 )^2 )

    # accumulate deformations associated with each edge
    x = x + 1*warp1 + cos(2π/3) * warp2 + cos(4π/3) * warp3
    y = y + 0*warp1 + sin(2π/3) * warp2 + sin(4π/3) * warp3

    return x,y
end


"""
xytors(x,y)

# Description

-    maps the (x,y) coordinates of an equilateral triangle to the standard right triangle

# Arguments

-    `x`: x-coordinate values in the equilateral triangle, an array
-    `y`: y-coordinate values in the equilateral triangle, an array

# Outputs: r,s

-    `r`: r-coordinate values in the standard right triangle, an array
-    `s`: s-coordinate values in the standard right triangle, an array

"""
function xytors(x,y)
    # compute lambda vectors
    L1 = @. (sqrt(3.0) * y + 1) / 3.0
    L2 = @. (-3.0 * x - sqrt(3.0) * y + 2.0) / 6.0
    L3 = @. ( 3.0 * x - sqrt(3.0) * y + 2.0) / 6.0

    # compute right triangle coordinates
    r = -L2 + L3 - L1
    s = -L2 - L3 + L1

    return r,s
end

"""
simplex2DP

# Description

- Evaluate the Jacobi polynomials at nodal values on the a,b grid

# Inputs

- `a`: first coordinate
- `b`: second coordinate
- `i`: jacobi polynomial parameter, Legendre hard coded
- `j`: jacobi polynomial parameter, Legendre hard coded

# Outputs

- `p`: value of the jacobi polynomial

"""
function simplex2DP(a,b, i::Int, j::Int)
    h1 = @. jacobi(a, 0, 0, i)
    h2 = @. jacobi(b, 2*i + 1, 0, j)
    p = @. sqrt(2.0) * h1 * h2 * (1-b)^i
    return p
end

"""
gradsimplex2DP

# Description

- Evaluate the derivative Jacobi polynomials at nodal values on the a,b grid

# Inputs

- `a`: first coordinate
- `b`: second coordinate
- `i`: jacobi polynomial parameter, Legendre hard coded
- `j`: jacobi polynomial parameter, Legendre hard coded

# Outputs

- `dmodedr`: value of the derivative of the jacobi polynomial
- `dmodeds`: value of the derivative of the jacobi polynomial

"""
function gradsimplex2DP(a,b, id::Int, jd::Int)
    fa = @. jacobi(a, 0, 0, id)
    dfa = @. djacobi(a,0,0, id)
    fb = @. jacobi(b, 2*id + 1, 0, jd)
    dfb = @. djacobi(b, 2*id + 1, 0, jd)

    #r-derivative
    # ∂ / ∂_r = ∂a / ∂_r + ∂b /  ∂_r = (2 / (1-b)) ∂ / ∂_a
    dmodedr = dfa .* fb
    if id>0
        @. dmodedr *= ( 0.5 * (1-b) )^(id-1)
    end
    #s-derivative
    # ∂ / ∂_s = (1+a/2)/(1-b/2) ∂ / ∂_a + ∂ / ∂_b
    dmodeds = @. dfa * fb * (1+a) / 2
    if id>0
        @. dmodeds *= ( 0.5 * (1-b) )^(id-1)
    end

    tmp = @. dfb * ( (0.5 * (1-b)) )^id
    if id>0
        @. tmp -= 0.5 * id * fb * ( (0.5 * (1-b)) )^(id-1)
    end
    @. dmodeds += fa * tmp

    #normalize
    @. dmodedr *= 2^(id + 0.5)
    @. dmodeds *= 2^(id + 0.5)
    return dmodedr, dmodeds
end


"""
vandermonde2D(n,r,s)

# Description

    Matrix to convert from modal to nodal basis. Depends on Simplex2DP

# Arguments

- `n`:
- `r`:
- `s`:

# Outputs

- `V2D`: 2D vandermonde matrix

"""
function vandermonde2D(n,r,s)
    nGL = Int((n+1) * (n+2) / 2);
    V2D = zeros(length(r), nGL)
    a,b = rstoab(r,s)
    let sk = 1
        for i ∈ 0:n
            for j ∈ 0:(n - i)
                V2D[:, sk] = simplex2DP(a,b,i,j)
                sk += 1
            end
        end
    end
    return V2D
end

"""
dvandermonde2D(n,r,s)

# Description

- gradient of model basis (i,j) at (r,s) at order N

# Arguments

- `n`:
- `r`:
- `s`:

# Outputs

- `V2Dr`: 2D partial derivative vandermonde matrix
- `V2Ds`: 2D partial derivative vandermonde matrix

"""
function dvandermonde2D(n,r,s)
    nGL = Int((n+1) * (n+2) / 2);
    V2Dr = zeros(length(r), nGL)
    V2Ds = zeros(length(r), nGL)
    a,b = rstoab(r,s)
    let sk = 1
        for i ∈ 0:n
            for j ∈ 0:(n - i)
                V2Dr[:, sk], V2Ds[:, sk] = gradsimplex2DP(a,b,i,j)
                sk += 1
            end
        end
    end
    return V2Dr, V2Ds
end

"""
dmatrices2D(n,r,s)

# Description

- gradient of nodal basis (i,j) at (r,s) at order N

# Arguments

- `n`:
- `r`:
- `s`:
- `V`: vandermonde matrix in 2D

# Outputs

- `∂r`: partial derivative
- `∂s`: partial derivative

"""
function dmatrices2D(n, r, s, V)
     Vr, Vs = dvandermonde2D(n,r,s)
     ∂r = Vr / V
     ∂s = Vs / V
    return ∂r, ∂s
end


"""
lift_tri(n, fmask, r, s, V)

# NEEDS TO BE TESTED

# Description

- lift operation for computing surface integrals

# Arguments

- `n` :
- `fmask`:
- `r`:
- `s`:
- `V`: vandermonde matrix in 2D

# Outputs

- `lift`: the lift operator


"""
function lift_tri(n, fmask, r, s, V)
    nGL = Int( (n+1) * (n+2) /2 )
    nFP = n+1
    nFaces = 3 #it's a triangle
    ∮ = zeros(nGL, nFaces*nFP)

    #face 1
    edge1_mask = fmask[:,1]
    faceR = r[edge1_mask];
    v = vandermonde(faceR, 0, 0, n)
    mass_edge_1 = inv(v * v')
    ∮[edge1_mask, 1:nFP] = mass_edge_1

    #face 2
    edge2_mask = fmask[:,2]
    faceR = r[edge2_mask];
    v = vandermonde(faceR, 0, 0, n)
    mass_edge_2 = inv(v * v')
    ∮[edge2_mask, (nFP+1):(2*nFP)] = mass_edge_2

    #face 3
    edge3_mask = fmask[:,3]
    faceS = s[edge3_mask];
    v = vandermonde(faceS, 0, 0, n)
    mass_edge_3 = inv(v * v')
    ∮[edge3_mask, (2*nFP+1):(3*nFP)] = mass_edge_3

    return ∮
end


"""

geometricfactors2D(x, y, Dr, Ds)

# NEEDS TO BE TESTED

# Description

- Metric elements for local mappings of elements

# Arguments

- `x`:
- `y`:
- `Dr`:
- `Ds`:

# Outputs : rx, sx, ry, sy, J

- `rx`:
- `sx`:
- `ry`:
- `sy`:
- `J`: jacobian

"""
function geometricfactors2D(x, y, Dr, Ds)
    xr = Dr * x; xs = Ds * x; yr = Dr * y; ys = Ds * y;
    J = - xs .* yr + xr .* ys; #determinant
    rx = ys ./ J; sx = - yr ./ J; ry = - xs ./ J; sy = xr ./ J;
    return rx, sx, ry, sy, J
end

"""
normals2D(x, y, Dr, Ds, fmask, nFP, K)

# NOT TESTED

# Description

- Returns unit normal vector along each element

# Arguments

- lots

# Outputs

- `nx`: normal along x-direction
- `ny`: normal along y-direction
- `sJ`: amplitude of unnormalized vector

"""

function normals2D(x, y, Dr, Ds, fmask, nFP, K)
    xr = Dr * x; xs = Ds * x; yr = Dr * y; ys = Ds * y;
    J = - xs .* xr + xr .* ys; #determinant

    #interpolate geometric factors to face nodes
    fxr = xr[fmask[:],:]
    fxs = xs[fmask[:],:]
    fyr = yr[fmask[:],:]
    fys = ys[fmask[:],:]

    #build normals
    nx = zeros(3*nFP, K) #3 edges on a triangle
    ny = zeros(3*nFP, K) #3 edges on a triangle
    fid1 = collect(         1:nFP     )
    fid2 = collect(   (nFP+1):(2*nFP) )
    fid3 = collect( (2*nFP+1):(3*nFP) )

    # face 1

    nx[fid1, :] =  fyr[fid1, :]
    ny[fid1, :] = -fxr[fid1, :]

    # face 2

    nx[fid2, :] =  fys[fid2, :] - fyr[fid2, :]
    ny[fid2, :] = -fxs[fid2, :] + fxr[fid2, :]

    # face 3

    nx[fid3, :] = -fys[fid3, :]
    ny[fid3, :] =  fxs[fid3, :]

    #normalize
    sJ = @. sqrt(nx^2 + ny^2)
    @. nx /= sJ
    @. ny /= sJ
    return nx, ny, sJ
end


"""
triangle_connect2D(EtoV)

# Description

- different connectivity

# Arguments

-  `EtoV`: element to vertices map

# Output

-  `EtoE`: element to element map
-  `EtoF`: element to face map

# Comments


"""
function triangle_connect2D(EtoV)
    nFaces = 3
    K = size(EtoV, 1)
    number_nodes = maximum(EtoV)

    #first  k-rows correspond to face 1
    #second k-rows correspond to face 2
    #third  k-rows correspond to face 3
    fnodes = [EtoV[:,[1,2]]; EtoV[:,[2,3]]; EtoV[:,[3,1]];]
    fnodes = sort(fnodes, dims = 2 ) .- 1

    #default element to element and element to faces connectivity
    EtoE = Int.( collect(1:K) * ones(1, nFaces) )
    EtoF = Int.( ones(K,1) * collect(1:nFaces)' )

    id =  @. fnodes[:,1] * number_nodes + fnodes[:,2] + 1;
    spNodeToNode = Int.( [id collect(1:nFaces*K) EtoE[:] EtoF[:]] )

    # check
    sorted = sortslices(spNodeToNode, dims = 1)
    indices = findall(sorted[1:(end-1),1] .== sorted[2:end,1])

    #make links reflexive
    matchL = [sorted[indices,:] ; sorted[indices .+ 1 , :]]
    matchR = [sorted[indices .+ 1,:] ; sorted[indices , :]]

    # insert matches
    @. EtoE[matchL[:,2]] = matchR[:,3]
    @. EtoF[matchL[:,2]] = matchR[:,4]
    EtoE = Int.(EtoE)
    EtoF = Int.(EtoF)

    return EtoE, EtoF
end


"""
connect2D(EtoV)

# Description

    Build connectivity maps for an arbitrary Mesh2D
    (currently assumes all elements have same number of faces)

# Arguments

-   `EtoV`: element to vertex map

# Return Values

-   `EtoE`: element to element map
-   `EtoF`: element to face map

"""
function connect2D(EtoV)
    nFaces = 3
    K = size(EtoV, 1)

    total_faces = nFaces * K

    # list of local face to local vertex connections
    # face convention implied here
    # face 1 is associated with vn[1,2]
    # face 2 is associated with vn[2,3]
    # ...
    # face nFaces is associated with vn[nFaces, 1]
    vn = zeros(Int, nFaces, 2)
    for i in 1:nFaces
        j = i % nFaces + 1
        vn[i,:] = [i,j]
    end

    # build global face to node sparse array
    # this is done by placing two vertices at each face row
    FtoV = spzeros(Int, total_faces, maximum(EtoV))
    let sk = 1
        for k in 1:K
            for face in 1:nFaces
                @. FtoV[sk, EtoV[k, vn[face,:] ] ] = 1;
                sk += 1
            end
        end
    end

    # the code below is just a way to say that we have vertices matching

    # global face to global face sparse array
    FtoF = FtoV * FtoV' - 2I # gotta love julia

    #find complete face to face connections
    #this just says that if two vertices match then its the same face
    faces1, faces2 = findnz(FtoF .== 2)

    # convert face global number to element and face numbers
    element1 = @. floor(Int, (faces1 - 1) / nFaces ) + 1
    element2 = @. floor(Int, (faces2 - 1) / nFaces ) + 1

    face1 = @. mod((faces1 - 1) , nFaces ) + 1
    face2 = @. mod((faces2 - 1) , nFaces ) + 1

    # Rearrange into Nelement x Nfaces sized arrays
    #ind = diag( LinearIndices(ones(Int, K, nFaces))[element1,face1] ) # this line is a terrible idea.

    # fixed version, just needed to convert to linear indices : /
    ind = element1 + (face1 .- 1) .* K

    #note that a convection has been assumed for the faces here
    EtoE = collect(Int, 1:K) * ones(Int, 1, nFaces)
    EtoF = ones(Int, K, 1) * collect(Int, 1:nFaces)'



    # each row is an element,
    # each index in the row is the neighbor index
    EtoE[ind] = copy(element2)
    # each row is an element
    # each column is a face
    # the entries in a column j link the face of an element i
    # to the face of the element EtoE[i,j]
    EtoF[ind] = copy(face2)

    return EtoE, EtoF
end


"""
connect_periodic_2D(VX, VY, EtoV)

# Description

-    Build connectivity maps for an arbitrary Mesh2D
-    (currently assumes all elements have same number of faces)
-    for periodicity, its easy to mody for only one side periodic

# Arguments

-   `EtoV`: element to vertex map

# Return Values

-   `EtoE`: element to element map
-   `EtoF`: element to face map

"""
function connect_periodic_2D(VX, VY, EtoV)
    ax = minimum(VX)
    bx = maximum(VX)
    ay = minimum(VY)
    by = maximum(VY)
    xperiod = bx - ax
    yperiod = by - ay
    # build periodic index converter
    conv = collect(1:length(VX))
    #build association map to create EtoVp (periodic version)
    minindx = findall( VX .≈ ax )
    minindy = findall( VY .≈ ay )
    maxindx = findall( VX .≈ bx )
    maxindy = findall( VY .≈ by )

    #match up appropriate vertices, does not generalize to 3D
    leftface = sortslices([VY[minindx] minindx], dims = 1)[:,2]
    rightface = sortslices([VY[maxindx] maxindx], dims = 1)[:,2]
    bottomface = sortslices([VX[minindy] minindy], dims = 1)[:,2]
    topface = sortslices([VX[maxindy] maxindy], dims = 1)[:,2]

    nFaces = 3
    K = size(EtoV, 1)

    total_faces = nFaces * K

    # list of local face to local vertex connections
    # face convention implied here
    # face 1 is associated with vn[1,2]
    # face 2 is associated with vn[2,3]
    # ...
    # face nFaces is associated with vn[nFaces, 1]
    vn = zeros(Int, nFaces, 2)
    for i in 1:nFaces
        j = i % nFaces + 1
        vn[i,:] = [i,j]
    end

    # build global face to node sparse array
    # this is done by placing two vertices at each face row
    FtoV = spzeros(Int, total_faces, maximum(EtoV))
    let sk = 1
        for k in 1:K
            for face in 1:nFaces
                @. FtoV[sk, EtoV[k, vn[face,:] ] ] = 1;
                sk += 1
            end
        end
    end

    # now make correction for periodic case
    # should only need to loop over elements on boundary

    nFacesTotal = nFaces*K
    nVX = length(leftface) # == length(bottomface[:,2])
    for i in 1:nFacesTotal
        for k in 1:(nVX-1)
            vecL = Int.( leftface[k:k+1])
            vecR = Int.(rightface[k:k+1])

            # check if face i is vecL
            if sum(FtoV[i, vecL]) == 2
                # identify left face with right face
                @. FtoV[i, vecL] = 0
                @. FtoV[i, vecR] = 1
                dropzeros!(FtoV)
            end

            vecB = Int.(bottomface[k:k+1])
            vecT = Int.(   topface[k:k+1])

            # identify top face with bottom face
            if sum(FtoV[i, vecB]) == 2
                @. FtoV[i, vecB] = 0
                @. FtoV[i, vecT] = 1
                dropzeros!(FtoV)
            end
        end
    end

    # global face to global face sparse array
    FtoF = FtoV * FtoV' - 2I # gotta love julia

    #find complete face to face connections
    #this just says that if two vertices match then its the same face
    faces1, faces2 = findnz(FtoF .== 2)

    # convert face global number to element and face numbers
    element1 = @. floor(Int, (faces1 - 1) / nFaces ) + 1
    element2 = @. floor(Int, (faces2 - 1) / nFaces ) + 1

    face1 = @. mod((faces1 - 1) , nFaces ) + 1
    face2 = @. mod((faces2 - 1) , nFaces ) + 1

    # Rearrange into Nelement x Nfaces sized arrays
    #ind = diag( LinearIndices(ones(Int, K, nFaces))[element1,face1] ) # this line is a terrible idea.

    # fixed version, just needed to convert to linear indices : /
    ind = element1 + (face1 .- 1) .* K

    #note that a convection has been assumed for the faces here
    EtoE = collect(Int, 1:K) * ones(Int, 1, nFaces)
    EtoF = ones(Int, K, 1) * collect(Int, 1:nFaces)'



    # each row is an element,
    # each index in the row is the neighbor index
    EtoE[ind] = copy(element2)
    # each row is an element
    # each column is a face
    # the entries in a column j link the face of an element i
    # to the face of the element EtoE[i,j]
    EtoF[ind] = copy(face2)

    return EtoE, EtoF
end



"""

buildmaps2D(K, np, nfp, nFaces, fmask, EtoE, EtoF, x, y, VX, VY)

# Description

- connectivity matrices for element to elements and elements to face

# Arguments

-   `K`: number of elements
-   `np`: number of points within an element
-   `nfp`: 1
-   `nFaces`: 2
-   `fmask`: an element by element mask to extract edge values
-   `EtoE`: element to element connectivity
-   `EtoF`: element to face connectivity
-   `x`: Guass Lobatto points along x-direction
-   `y`: Guass Lobatto points along y-direction
-   `VX`: vertex stuff
-   `VY`: vertex stuff

# Return Values: vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO

-   `vmapM`: vertex indices, (used for interior u values)
-   `vmapP`: vertex indices, (used for exterior u values)
-   `vmapB`: vertex indices, corresponding to boundaries
-   `mapB`: use to extract vmapB from vmapM

"""
function buildmaps2D(K, np, nfp, nFaces, fmask, EtoE, EtoF, EtoV, x, y, VX, VY)
    # number volume nodes consecutively
    nodeids = reshape(collect(1:(K*np)), np, K)
    vmapM = zeros(Int, nfp, nFaces, K)
    vmapP = zeros(Int, nfp, nFaces, K)
    mapM = collect(1: (K*nfp*nFaces) )'
    mapP = copy( reshape(mapM, nfp, nFaces, K) )
    # find index of face nodes wrt volume node ordering
    for k1 in 1:K
        for f1 in 1:nFaces
            vmapM[:, f1, k1] = nodeids[fmask[:,f1], k1]
        end
    end

    let one1 = ones(1, nfp)
    for k1 = 1:K
        for f1 = 1:nFaces
            # find neighbor
            k2 = EtoE[k1, f1]
            f2 = EtoF[k1, f1]

            # reference length of edge
            v1 = EtoV[k1,f1]
            v2 = EtoV[k1, 1 + mod(f1,nFaces)]
            refd = @. sqrt( (VX[v1] - VX[v2])^2  + (VY[v1] - VY[v2])^2       )

            # find volume node numbers of left and right nodes
            vidM = vmapM[:, f1, k1]
            vidP = vmapM[:, f2, k2]

            x1 = x[vidM]; y1 = y[vidM]
            x2 = x[vidP]; y2 = y[vidP]

            # may need to reshape before multiplication
            x1 = x1 * one1
            y1 = y1 * one1
            x2 = x2 * one1
            y2 = y2 * one1

            # compute distance matrix
            D = @. (x1 - x2' )^2 + (y1 - y2' )^2
            mask = @. D < eps(refd)

            #find linear indices
            m,n = size(D)
            d = collect(1:(m*n))
            idM =  @. Int( floor( (d[mask[:]]-1) / m) + 1 )
            idP =  @. Int( mod( (d[mask[:]]-1) , m) + 1 )
            vmapP[idM, f1, k1] = vidP[idP]
            @. mapP[idM, f1, k1] = idP + (f2-1)*nfp + (k2-1)*nFaces*nfp
        end
    end
    end
        # reshape arrays
        vmapP = Int.( reshape(vmapP, length(vmapP)) )
        vmapM = Int.( reshape(vmapM, length(vmapM)) )

        # Create list of boundary nodes
        mapB = Int.( collect(1:length(vmapP))[vmapP .== vmapM] )
        vmapB = Int.( vmapM[mapB] )

        return vmapM, vmapP, vmapB, mapB
end

"""

build_periodic_maps2D(K, np, nfp, nFaces, fmask, EtoE, EtoF, x, y, VX, VY)

# Description

- connectivity matrices for element to elements and elements to face. assumes periodic

# Arguments

-   `K`: number of elements
-   `np`: number of points within an element (polynomial degree + 1)
-   `nfp`: 1
-   `nFaces`: 2
-   `fmask`: an element by element mask to extract edge values
-   `EtoE`: element to element connectivity
-   `EtoF`: element to face connectivity
-   `x`: Guass Lobatto points along x-direction
-   `y`: Guass Lobatto points along y-direction
-   `VX`: vertex stuff
-   `VY`: vertex stuff

# Return Values: vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO

-   `vmapM`: vertex indices, (used for interior u values)
-   `vmapP`: vertex indices, (used for exterior u values)
-   `vmapB`: vertex indices, corresponding to boundaries
-   `mapB`: use to extract vmapB from vmapM

"""
function build_periodic_maps2D(K, np, nfp, nFaces, fmask, EtoE, EtoF, EtoV, x, y, VX, VY)

    xperiod = maximum(x)-minimum(x)
    yperiod = maximum(y)-minimum(y)
    # number volume nodes consecutively
    nodeids = reshape(collect(1:(K*np)), np, K)
    vmapM = zeros(Int, nfp, nFaces, K)
    vmapP = zeros(Int, nfp, nFaces, K)
    mapM = collect(1: (K*nfp*nFaces) )'
    mapP = copy( reshape(mapM, nfp, nFaces, K) )
    # find index of face nodes wrt volume node ordering
    for k1 in 1:K
        for f1 in 1:nFaces
            vmapM[:, f1, k1] = nodeids[fmask[:,f1], k1]
        end
    end

    let one1 = ones(1, nfp)
    for k1 = 1:K
        for f1 = 1:nFaces
            # find neighbor
            k2 = EtoE[k1, f1]
            f2 = EtoF[k1, f1]

            # reference length of edge
            v1 = EtoV[k1,f1]
            v2 = EtoV[k1, 1 + mod(f1,nFaces)]
            refd = @. sqrt( ( (VX[v1] - VX[v2]))^2  + ( (VY[v1] - VY[v2]))^2       )

            # find volume node numbers of left and right nodes
            # size is number of face points
            vidM = vmapM[:, f1, k1]
            vidP = vmapM[:, f2, k2]

            # get the elements on the face
            x1 = x[vidM]; y1 = y[vidM]
            x2 = x[vidP]; y2 = y[vidP]

            # may need to reshape before multiplication
            x1 = x1 * one1
            y1 = y1 * one1
            x2 = x2 * one1
            y2 = y2 * one1
            #the above is just a convenient way to do all pairwise
            #calculations in terms of matrices

            # compute distance matrix
            D = @. ( x1 - x2')^2 + ( y1 - y2' )^2
            #to handle periodic boundary conditions
            Dpx = @. ( abs(x1 - x2') - xperiod )^2 + ( abs(y1 - y2') )^2
            Dpy = @. ( abs(x1 - x2') )^2 + ( abs(y1 - y2') - yperiod )^2
            mask = @. ( D < eps(refd) ) | ( Dpx < eps(refd) )  | ( Dpy < eps(refd) )

            #find linear indices
            m,n = size(D)
            d = collect(1:(m*n))

            #should be of size nfp
            idM =  @. Int( floor( (d[mask[:]]-1) / m) + 1 )
            idP =  @. Int( mod( (d[mask[:]]-1) , m) + 1 )
            vmapP[idM, f1, k1] = vidP[idP]
            @. mapP[idM, f1, k1] = idP + (f2-1)*nfp + (k2-1)*nFaces*nfp
        end
    end
    end
        # reshape arrays
        vmapP = Int.( reshape(vmapP, length(vmapP)) )
        vmapM = Int.( reshape(vmapM, length(vmapM)) )

        # Create list of boundary nodes
        mapB = Int.( collect(1:length(vmapP))[vmapP .== vmapM] )
        vmapB = Int.( vmapM[mapB] )

        return vmapM, vmapP, vmapB, mapB
end


#=
for debugging
EtoE = EtoEp
EtoF = EtoFp
vmapP = vmapPp
# use after loop has been run
findall(vmapP .== 0)
one1 = ones(1, nfp)
# the last two indices are f1 and k1 respecitively
# inside the loop
f1 = 1
k1 = 1
scatter(x1,y1)
scatter!(x2,y2, legend = false)
# check idP dimension, should be nfp
=#

#= could also do the mask
mask2 = zeros(Bool,nfp,nfp)
for i in 1:nfp
    for j in 1:nfp
        mask2[i,j] =  ( ( (x1[i] - x2[j])%xperiod )^2 + ( (y1[i] - y2[j])%yperiod )^2 ) < eps(refd)
    end
end
#should give the number length(mask) if all true
sum(mask .== mask2)
=#



"""
global_grid(r, s, EtoV, VX, VY)

# Description

- Create a global grid from the elements

# Arguments

- `r` : ideal coordinates
- `s` : ideal coordinates
- `EtoV` : element to vertices
- `VX` : x-coordinate of vertices
- `VY` : y-coordinate of vertices

# Return x, y

- `x` : x-coordinates of grid
- `y` : y-coordinates of grid

"""
function global_grid(r, s, EtoV, VX, VY)
    #need to transpose so that the matrix multiply works
    va = EtoV[:,1]'
    vb = EtoV[:,2]'
    vc = EtoV[:,3]'
    # global x and y values constructed from ideal coordinates and grid
    x =  0.5 * ( - (r+s) * VX[va] + (1 .+ r)*VX[vb] + (1 .+s)*VX[vc])
    y =  0.5 * ( - (r+s) * VY[va] + (1 .+ r)*VY[vb] + (1 .+s)*VY[vc])
    return x, y
end

"""
create_fmask(r,s)

# Description

- mask to get edge nodes

# Arguments

- `r` : ideal coordinates
- `s` : ideal coordinates

#Return

-  `fmask`: array of boolean values of edge nodes

"""
function create_fmask(r, s)
    fmask1 = findall( abs.( s .+ 1) .< eps(10.0) )'
    fmask2 = findall( abs.( r .+ s ) .< eps(10.0) )'
    fmask3 = findall( abs.( r .+ 1) .< eps(10.0) )'
    fmask = Array([fmask1; fmask2; fmask3]')
    return fmask
end

"""
find_edge_nodes(fmask, x, y)

# Description

- find the values on the grid that correspond to edges

# Arguments

- `fmask` : mask to extract edge values, use the function create_fmask(r,s)
-  `x` : x-coordinates of grid
-  `y` : y-coordinates of grid

# Return

- `edge_x`: x-coordinate of edge values (called Fx in NDG)
- `edge_y`: y-coorsdinate of edge values (called Fy in NDG)

"""
function find_edge_nodes(fmask, x, y)
    edge_x = x[fmask[:], :]
    edge_y = y[fmask[:], :]
    return edge_x, edge_y
end

"""
build_bc_maps(mesh, bctype, bc_label)

# Description

- find the indices corresponding to boundary data

# Arguments

- `mesh` : mesh struct
-  `bctype` : matrix corresponding to boundary conditions (sime size as EToV)
-  `bc_label` : label corresponding to index (this is for the user)

# Return

- `mapT`: all the maps, including inflow and outflow maps
- `vmapT`: all the vmaps, including inflow and outlfow global indices
- `newbc_label`: labels for indices

# Example

- the array mapT[index] corresponds to boundary condition type newbc_label[i]

"""
function build_bc_maps(mesh, bctype, bc_label)
    bct = bctype'
    bnodes = ones(mesh.nfp) *  bct[:]'
    btotal_ind = setdiff(union(bctype)[:],[0])
    newbc_label = bc_label[btotal_ind]
    mapT = []
    vmapT = []
    for j in btotal_ind
        bc = findall(bnodes[:] .== j)
        push!(mapT,bc)
        push!(vmapT, mesh.vmapM[bc])
    end
    return mapT, vmapT, newbc_label
end



struct garbage_triangle3{T, S, U, W, V}
    # inputs
    n::S
    filename::V

    # face stuff
    nfp::S
    nFaces::S
    K::S

    # GL points
    r::U
    s::U
    x::T
    y::T

    # vertex maps
    vmapM::W
    vmapP::W
    vmapB::W
    mapB::W

    # structures for computation
    J::T
    sJ::T
    Dʳ::T
    Dˢ::T
    Drw::T
    Dsw::T
    M::T
    Mi::T
    lift::T
    rx::T
    ry::T
    sx::T
    sy::T
    nx::T
    ny::T
    fscale::T

    function garbage_triangle3(n, filename)
        ℳ = meshreader_gambit2D(filename)

        #get number of points
        nfp = n+1
        np = Int( (n+1)*(n+2)/2 )
        nFaces = 3
        nodetol = 1e-12

        #compute nodal set
        x, y = nodes2D(n)
        r, s = xytors(x,y)

        #build reference elements and matrices
        V = vandermonde2D(n, r, s)
        invV = inv(V)
        Mi = V * V'
        M = invV' * invV
        Dʳ, Dˢ = dmatrices2D(n , r, s, V)

        # build global grid
        x,y = global_grid(r, s, ℳ.EtoV, ℳ.vertices[:,1], ℳ.vertices[:,2])

        # create fmask
        fmask = create_fmask(r, s)
        edge_x, edge_y = find_edge_nodes(fmask, x, y)
        lift = lift_tri(n, fmask, r, s, V)

        rx, sx, ry, sy, J = geometricfactors2D(x, y, Dʳ, Dˢ)

        nx, ny, sJ = normals2D(x, y, Dʳ, Dˢ, fmask, nfp, ℳ.K)
        fscale = sJ ./ J[fmask[:],:]

        #EtoE, EtoF = connect2D(EtoV)
        EtoE, EtoF = triangle_connect2D(ℳ.EtoV)

        vmapM, vmapP, vmapB, mapB = buildmaps2D(ℳ.K, np, nfp, nFaces, fmask, EtoE, EtoF, ℳ.EtoV, x, y, ℳ.vertices[:,1], ℳ.vertices[:,2])
        Vr, Vs = dvandermonde2D(n,r,s)
        Drw = (V*Vr') / (V*V')
        Dsw = (V*Vs') / (V*V')
        return new{typeof(x),typeof(K),typeof(r),typeof(vmapP),typeof(filename)}(n, filename, nfp, nFaces, ℳ.K, r,s, x,y, vmapM,vmapP,vmapB,mapB, J, sJ, Dʳ, Dˢ, Drw, Dsw, M, Mi, lift, rx, ry, sx, sy, nx, ny, fscale)
    end
end


struct periodic_triangle{T, S, U, W, V}
    # inputs
    n::S
    filename::V

    # face stuff
    nfp::S
    nFaces::S
    K::S

    # GL points
    r::U
    s::U
    x::T
    y::T

    # vertex maps
    vmapM::W
    vmapP::W
    vmapB::W
    mapB::W

    # structures for computation
    J::T
    sJ::T
    Dʳ::T
    Dˢ::T
    Drw::T
    Dsw::T
    M::T
    Mi::T
    lift::T
    rx::T
    ry::T
    sx::T
    sy::T
    nx::T
    ny::T
    fscale::T

    function periodic_triangle(n, filename)
        Nv, VX, VY, K, EtoV = meshreader_gambit2D(filename)

        #get number of points
        nfp = n+1
        np = Int( (n+1)*(n+2)/2 )
        nFaces = 3
        nodetol = 1e-12

        #compute nodal set
        x, y = nodes2D(n)
        r, s = xytors(x,y)

        #build reference elements and matrices
        V = vandermonde2D(n, r, s)
        invV = inv(V)
        Mi = V * V'
        M = invV' * invV
        Dʳ, Dˢ = dmatrices2D(n , r, s, V)

        # build global grid
        x,y = global_grid(r, s, EtoV, VX, VY)

        # create fmask
        fmask = create_fmask(r, s)
        edge_x, edge_y = find_edge_nodes(fmask, x, y)
        lift = lift_tri(n, fmask, r, s, V)

        rx, sx, ry, sy, J = geometricfactors2D(x, y, Dʳ, Dˢ)

        nx, ny, sJ = normals2D(x, y, Dʳ, Dˢ, fmask, nfp, K)
        fscale = sJ ./ J[fmask[:],:]

        EtoE, EtoF = connect_periodic_2D(VX, VY, EtoV)
        #EtoE, EtoF = triangle_connect2D(EtoVp)

        vmapM, vmapP, vmapB, mapB = build_periodic_maps2D(K, np, nfp, nFaces, fmask, EtoE, EtoF, EtoV, x, y, VX, VY)
        Vr, Vs = dvandermonde2D(n,r,s)
        Drw = (V*Vr') / (V*V')
        Dsw = (V*Vs') / (V*V')
        return new{typeof(x),typeof(K),typeof(r),typeof(vmapP),typeof(filename)}(n, filename, nfp, nFaces, K, r,s, x,y, vmapM,vmapP,vmapB,mapB, J, sJ, Dʳ, Dˢ, Drw, Dsw, M, Mi, lift, rx, ry, sx, sy, nx, ny, fscale)
    end
end

struct dg_garbage_triangle{T}
    u::T
    u̇::T
    φˣ::T
    φʸ::T
    fˣ::T
    fʸ::T
    fⁿ::T
    """
    dg_triangle(mesh)

    # Description

        initialize dg struct

    # Arguments

    -   `mesh`: a mesh to compute on

    # Return Values:

    -   `u` : the field to be computed
    -   `u̇`: numerical solutions for the field
    -   `φˣ`: x-component of flux
    -   `φʸ`: y-component of flux
    -   `fˣ`: the numerical flux on face in the x-direction for the computation
    -   `fʸ`: the numerical flux on face in the y-direction for the computation
    -   `fⁿ`: the numerical flux on face in the normal direction for the computation

    """
    function dg_garbage_triangle(mesh)
        # set up the solution
        u   = similar(mesh.x)
        u̇   = similar(mesh.x)
        φˣ  = similar(mesh.x)
        φʸ  = similar(mesh.x)
        fˣ  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        fʸ  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        fⁿ  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        return new{typeof(u)}(u, u̇, φˣ, φʸ, fˣ, fʸ, fⁿ)
    end
end

"""
make_periodic(VX, VY, EtoV)

# Description

- makes a grid periodic. note that this won't work with any grid. In particular the domain must be recangular and the vertices must exist on both sides. Does not work for small grids (with small number of elements)

# Inputs

-   `VX`: x-coordinate of vertices
-   `VY`: y-coordinate of vertices
-   `EtoV`: element to vertex connection

# Outputs

-   `EtoVp`: periodic element to vertex connections

"""
function make_periodic_2D(VX,VY, EtoV)
    ax = minimum(VX)
    bx = maximum(VX)
    ay = minimum(VY)
    by = maximum(VY)
    xperiod = bx - ax
    yperiod = by - ay
    # build periodic index converter
    conv = collect(1:length(VX))
    #build association map to create EtoVp (periodic version)
    minindx = findall( VX .≈ ax )
    minindy = findall( VY .≈ ay )
    maxindx = findall( VX .≈ bx )
    maxindy = findall( VY .≈ by )

    bool_corner1 = @. (VX==bx) & (VY==by)
    corner1 = findall(  bool_corner1 )

    bool_corner2 = @. (VX==ax) & (VY==by)
    corner2 = findall(  bool_corner2 )

    bool_corner3 = @. (VX==bx) & (VY==ay)
    corner3 = findall(  bool_corner3 )

    bool_corner4 = @. (VX==ax) & (VY==ay)
    corner4 = findall(  bool_corner4 )

    bool_corners = @. ((VX==ax)|(VX==bx)) & ((VY==ay)|(VY==by))
    corners = findall(bool_corners)

    # potentially safer
    corner1 = findall(  bool_corner1 )
    corner2 = findall(  bool_corner2 )
    corner3 = findall(  bool_corner3 )
    corner4 = findall(  bool_corner4 )

    minindt = [minindx; minindy]
    maxindt = [maxindx; maxindy]

    for i in 1:length(minindt)
        ii = minindt[i]
        for j in 1:length(maxindt)
            jj = maxindt[j]
            xtrue = (VX[jj]-VX[ii])%xperiod ≈ 0
            ytrue = (VY[jj]-VY[ii])%yperiod ≈ 0
            if xtrue && ytrue
                conv[ii] = jj
            end
        end
    end

    # make sure there is only one corner

    for i in 1:length(VX)
        if conv[i] in corners
            conv[i] = corners[1]
        end
    end
    EtoVp = copy(EtoV)
    for i in 1:length(EtoVp)
        EtoVp[i] = conv[ EtoVp[i] ] # convert to appropriate vertex
    end
    return EtoVp
end


"""
plot_mesh

# Description

- plots the mesh

"""
function plot_mesh(mesh)
    p1 = scatter(mesh.x, mesh.y, legend=false)
    # plot boundary of triangles
    scatter!(mesh.x[mesh.vmapM] , mesh.y[mesh.vmapM], color = "black", legend = false)
    #plot boundary of domain
    scatter!(mesh.x[mesh.vmapB] , mesh.y[mesh.vmapB], color = "yellow", legend = false)
    display(plot(p1))
end
