include("../utils.jl")
include("mesh2D.jl")

using SparseArrays

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
    np = length(r);
    a = zeros(np)
    for i in 1:np
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
    α° = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258];
    if n < 16
        α = α°[n]
    else
        α = 5/3
    end
    np = Int((n+1) * (n+2) / 2);
    L1 = zeros(np)
    L2 = zeros(np)
    L3 = zeros(np)
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
    blend1 = 4 * L2 .* L3
    blend2 = 4 * L3 .* L1
    blend3 = 4 * L1 .* L2
    warpf1 = warp_factor(n, L3 - L2)
    warpf2 = warp_factor(n, L1 - L3)
    warpf3 = warp_factor(n, L2 - L1)
    warp1 = @. blend1 * warpf1 * (1 + (α * L1 )^2 )
    warp2 = @. blend2 * warpf2 * (1 + (α * L2 )^2 )
    warp3 = @. blend3 * warpf3 * (1 + (α * L3 )^2 )
    x = x + 1*warp1 + cos(2π/3) * warp2 + cos(4π/3) * warp3
    y = y + 0*warp1 + sin(2π/3) * warp2 + sin(4π/3) * warp3
    return x, y
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
    L1 = @. (sqrt(3.0) * y + 1) / 3.0
    L2 = @. (-3.0 * x - sqrt(3.0) * y + 2.0) / 6.0
    L3 = @. ( 3.0 * x - sqrt(3.0) * y + 2.0) / 6.0
    r = -L2 + L3 - L1
    s = -L2 -L3 + L1
    return r, s
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
    np = Int((n+1) * (n+2) / 2);
    V2D = zeros(length(r), np)
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
    np = Int((n+1) * (n+2) / 2);
    V2Dr = zeros(length(r), np)
    V2Ds = zeros(length(r), np)
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
    np = Int( (n+1) * (n+2) /2 )
    nfp = n+1
    nfaces = 3 #it's a triangle
    ℰ = zeros(np, nfaces*nfp)

    #face 1
    edge1_mask = fmask[:,1]
    faceR = r[edge1_mask];
    v = vandermonde(faceR, 0, 0, n)
    mass_edge_1 = inv(v * v')
    ℰ[edge1_mask, 1:nfp] = mass_edge_1

    #face 2
    edge2_mask = fmask[:,2]
    faceR = r[edge2_mask];
    v = vandermonde(faceR, 0, 0, n)
    mass_edge_2 = inv(v * v')
    ℰ[edge2_mask, (nfp+1):(2*nfp)] = mass_edge_2

    #face 3
    edge3_mask = fmask[:,3]
    faceS = s[edge3_mask];
    v = vandermonde(faceS, 0, 0, n)
    mass_edge_3 = inv(v * v')
    ℰ[edge3_mask, (2*nfp+1):(3*nfp)] = mass_edge_3

    #compute lift
    lift = V * (V' * ℰ)
    return lift
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
normals2D(x, y, Dr, Ds, fmask, nfp, K)

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

function normals2D(x, y, Dr, Ds, fmask, nfp, K)
    xr = Dr * x; xs = Ds * x; yr = Dr * y; ys = Ds * y;
    J = - xs .* xr + xr .* ys; #determinant

    #interpolate geometric factors to face nodes
    fxr = xr[fmask[:],:]
    fxs = xs[fmask[:],:]
    fyr = yr[fmask[:],:]
    fys = ys[fmask[:],:]

    #build normals
    nx = zeros(3*nfp, K) #3 edges on a triangle
    ny = zeros(3*nfp, K) #3 edges on a triangle
    fid1 = collect(         1:nfp     )
    fid2 = collect(   (nfp+1):(2*nfp) )
    fid3 = collect( (2*nfp+1):(3*nfp) )

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
connect2D(EToV)

# Description

-

# Arguments

-  `EToV`: element to vertices map

# Output

-  `EToE`: element to element map
-  `EToF`: element to face map

# Comments

- The changes from the 1D are minor. nfaces can probably remain generic

"""
function connect2D(EToV)
    nfaces = 3 # because ... triangles

    #find number of elements and vertices
    K = size(EToV, 1)
    Nv = maximum(EToV)

    # create face to node connectivity matrix
    total_faces = nfaces * K

    # list of local face to local vertex connections
    vn = [[1 2]; [2 3]; [1 3] ]

    # build global face to node sparse array
    SpFToV = spzeros(Int, total_faces, Nv)
    let sk = 1
        for k ∈ 1:K
            for face ∈ 1:nfaces
                @. SpFToV[sk, EToV[k, vn[face,:] ] ] = 1;
                sk += 1
            end
        end
    end

    # global face to global face sparse array
    SpFToF = SpFToV * SpFToV' - 2I #gotta love julia

    #find complete face to face connections
    faces1, faces2 = findnz(SpFToF .== 2)

    # convert face global number to element and face numbers
    element1 = @. Int( floor( (faces1 - 1) / nfaces ) + 1 )
    element2 = @. Int( floor( (faces2 - 1) / nfaces ) + 1 )

    face1 = @. Int( mod( (faces1 - 1) , nfaces ) + 1 )
    face2 = @. Int( mod( (faces2 - 1) , nfaces ) + 1 )

    # Rearrange into Nelement x Nfaces sized arrays
    ind = diag( LinearIndices(ones(K, nfaces))[element1,face1] ) # this line is a terrible idea.
    EToE = collect(1:K) * ones(1, nfaces)
    EToF = ones(K,1) * collect(1:nfaces)'
    EToE[ind] = copy(element2);
    EToF[ind] = copy(face2);
    return EToE, EToF
end


"""
triangle_connect2D(EToV)

# Description

- different connectivity

# Arguments

-  `EToV`: element to vertices map

# Output

-  `EToE`: element to element map
-  `EToF`: element to face map

# Comments


"""
function triangle_connect2D(EToV)
    nfaces = 3
    K = size(EToV, 1)
    number_nodes = maximum(EToV)

    fnodes = [EToV[:,[1,2]]; EToV[:,[2,3]]; EToV[:,[3,1]];]
    fnodes = sort(fnodes, dims = 2 ) .- 1

    #default element to element and element to faces connectivity
    EToE = collect(1:K) * ones(1, nfaces)
    EToF = ones(K,1) * collect(1:nfaces)'

    id =  @. fnodes[:,1] * number_nodes + fnodes[:,2] + 1;
    spNodeToNode = Int.( [id collect(1:nfaces*K) EToE[:] EToF[:]] )

    # check
    sorted = sortslices(spNodeToNode, dims = 1)
    indices = findall(sorted[1:(end-1),1] .== sorted[2:end,1])

    #make links reflexive
    matchL = [sorted[indices,:] ; sorted[indices .+ 1 , :]]
    matchR = [sorted[indices .+ 1,:] ; sorted[indices , :]]

    # insert matches
    @. EToE[matchL[:,2]] = matchR[:,3]
    @. EToF[matchL[:,2]] = matchR[:,4]
    EToE = Int.(EToE)
    EToF = Int.(EToF)

    return EToE, EToF
end


"""

buildmaps2D(K, np, nfp, nfaces, fmask, EToE, EToF, x, y, VX, VY)

# Description

- connectivity matrices for element to elements and elements to face

# Arguments

-   `K`: number of elements
-   `np`: number of points within an element (polynomial degree + 1)
-   `nfp`: 1
-   `nfaces`: 2
-   `fmask`: an element by element mask to extract edge values
-   `EToE`: element to element connectivity
-   `EToF`: element to face connectivity
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
function buildmaps2D(K, np, nfp, nfaces, fmask, EToE, EToF, EToV, x, y, VX, VY)
    # number volume nodes consecutively
    nodeids = reshape(collect(1:(K*np)), np, K)
    vmapM = zeros(Int, nfp, nfaces, K)
    vmapP = zeros(Int, nfp, nfaces, K)
    mapM = collect(1: (K*nfp*nfaces) )'
    mapP = copy( reshape(mapM, nfp, nfaces, K) )
    # find index of face nodes wrt volume node ordering
    for k1 in 1:K
        for f1 in 1:nfaces
            vmapM[:, f1, k1] = nodeids[fmask[:,f1], k1]
        end
    end

    let one1 = ones(1, nfp)
    for k1 = 1:K
        for f1 = 1:nfaces
            # find neighbor
            k2 = EToE[k1, f1]
            f2 = EToF[k1, f1]

            # reference length of edge
            v1 = EToV[k1,f1]
            v2 = EToV[k1, 1 + mod(f1,nfaces)]
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
            @. mapP[idM, f1, k1] = idP + (f2-1)*nfp + (k2-1)*nfaces*nfp
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
global_grid(r, s, EToV, VX, VY)

# Description

- Create a global grid from the elements

# Arguments

- `r` : ideal coordinates
- `s` : ideal coordinates
- `EToV` : element to vertices
- `VX` : x-coordinate of vertices
- `VY` : y-coordinate of vertices

# Return x, y

- `x` : x-coordinates of grid
- `y` : y-coordinates of grid

"""
function global_grid(r, s, EToV, VX, VY)
    #need to transpose so that the matrix multiply works
    va = EToV[:,1]'
    vb = EToV[:,2]'
    vc = EToV[:,3]'
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
    fmask = [fmask1; fmask2; fmask3]'
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
    edge_x = x[fmask[:],:]
    edge_y = y[fmask[:],:]
    return edge_x, edge_y
end


struct garbage_triangle3{T, S, U, W, V}
    # inputs
    n::S
    filename::V

    # face stuff
    nfp::S
    nfaces::S
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
        Nv, VX, VY, K, EToV = meshreader_gambit2D(filename)

        #get number of points
        nfp = n+1
        np = Int( (n+1)*(n+2)/2 )
        nfaces = 3
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
        x,y = global_grid(r, s, EToV, VX, VY)

        # create fmask
        fmask = create_fmask(r, s)
        edge_x, edge_y = find_edge_nodes(fmask, x, y)
        lift = lift_tri(n, fmask, r, s, V)

        rx, sx, ry, sy, J = geometricfactors2D(x, y, Dʳ, Dˢ)

        nx, ny, sJ = normals2D(x, y, Dʳ, Dˢ, fmask, nfp, K)
        fscale = sJ ./ J[fmask[:],:]

        #EToE, EToF = connect2D(EToV)
        EToE, EToF = triangle_connect2D(EToV)

        vmapM, vmapP, vmapB, mapB = buildmaps2D(K, np, nfp, nfaces, fmask, EToE, EToF, EToV, x, y, VX, VY)
        Vr, Vs = dvandermonde2D(n,r,s)
        Drw = (V*Vr') / (V*V')
        Dsw = (V*Vs') / (V*V')
        return new{typeof(x),typeof(K),typeof(r),typeof(vmapP),typeof(filename)}(n, filename, nfp, nfaces, K, r,s, x,y, vmapM,vmapP,vmapB,mapB, J, sJ, Dʳ, Dˢ, Drw, Dsw, M, Mi, lift, rx, ry, sx, sy, nx, ny, fscale)
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
        fˣ  = zeros(mesh.nfp * mesh.nfaces, mesh.K)
        fʸ  = zeros(mesh.nfp * mesh.nfaces, mesh.K)
        fⁿ  = zeros(mesh.nfp * mesh.nfaces, mesh.K)
        return new{typeof(u)}(u, u̇, φˣ, φʸ, fˣ, fʸ, fⁿ)
    end
end
