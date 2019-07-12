include("../src/utils.jl")
include("element2D.jl")
include("rectangles.jl")

using SparseArrays
using LinearAlgebra

abstract type AbstractMesh2D end
abstract type AbstractGrid2D <: AbstractMesh2D end

"""
Mesh2D(K, vertices, EtoV, nFaces)

# Description

    initialize a Mesh2D struct which contains only information about the number of elements, vertices, and the mapping between the two

    No information about GL points is created or stored at this point

# Arguments

-   `K`: number of elements in the grid
-   `vertices`: array of all the vertices in the grid
-   `EtoV`: element to vertex mapping
-   `nFaces`: number of faces of the elements, will be deprecated once connect2D() and buildmaps2D() are flexible enough to handle complex grids


# Return Values:

    return a properly initiliazed Mesh2D object

"""
struct Mesh2D{S, T, U} <: AbstractMesh2D
    K::S
    vertices::T
    EtoV::U
    nFaces::S

    function Mesh2D(K, vertices, EtoV, nFaces)
        return new{typeof(K),typeof(vertices),typeof(EtoV)}(K, vertices, EtoV, nFaces)
    end
end


"""
Grid2D(ℳ::Mesh2D, N::Int) ℳ

# Description

    initialize a Grid2D object with fully constructed elements, GL points, and maps between them. Can be used to perform a computation

# Arguments

-   `ℳ`: a Mesh2D object to fill in with GL points
-   `N`: polynomial order along each face, is the same for each face and each element for now

# Return Values:

    return a properly initiliazed Grid2D object

"""
struct Grid2D{S, T, U, V, W} <: AbstractGrid2D
    # basic grid of vertices
    ℳ::S

    # elements
    Ω::T

    # GL points
    nGL::U # number of points
    x::V # physical coordinates

    # maps maps maps
    vmap⁻::W
    vmap⁺::W
    vmapᴮ::W
    mapᴮ::W

    function Grid2D(ℳ::Mesh2D, N::Int)
        # get elements and GL nodes
        Ω,x̃ = makenodes2D(ℳ, N)
        nGL = length(x̃[:,1])


        # make a facemask
        fmask = Ω[1].fmask # same as with nGLᵏ, should eventually be done on an element by element basis
        nFP,nFaces = size(fmask)
        nGLᵏ = length(Ω[1].x[:, 1]) # get number of GL points from first element, buildmaps2D currently needs all elements to map to the same ideal element

        # build the boundary maps
        vmap⁻,vmap⁺,vmapᴮ,mapᴮ = buildmaps2D(ℳ, nFP, nGLᵏ, fmask, x̃)

        return new{typeof(ℳ),typeof(Ω),typeof(nGL),typeof(x̃),typeof(vmap⁻)}(ℳ, Ω, nGL,x̃, vmap⁻,vmap⁺,vmapᴮ,mapᴮ)
    end
end

"""
meshreader_gambit2D(filename)

# Description

- Reads the .neu files from NDG

# Arguments

- `filename`: a string that contains the path to the file

# Return : Nv, VX, VY, K, EtoV

- `Nv` :   number of vertices
- `VX` :   x coordinate of vertex
- `VY` :   y coordinate of vertex
- `K`  :   number of elements
- `EtoV` : element to vertex connection

"""
function meshreader_gambit2D(_filename)
    #open file
    f = open(_filename)
    #get things line by line, lines[1] is the first line of the file
    lines = readlines(f)

    # in the standard .neu format line 7 contains the data for
    # reading the the number of edges and vertices
    # lines[6] contains the data for whats going on
    #println(lines[6])
    #println(lines[7])
    data = lines[7]
    #split up the blank spaces
    dims = split(data)
    #reinterpret the strings as integers
    dims = map(x->parse(Int,x), dims)

    #Now we can extract Nv and K
    Nv = dims[1]
    K  = dims[2]

    #first get node coordinates
    VX = ones(Nv)
    VY = ones(Nv)

    # the lines with data start at lines[10]
    # the lines end with lines[10+Nv]
    for i ∈ 1:Nv
        data = lines[9+i]
        #split up the blank spaces
        dims = split(data)
        #reinterpret the strings as floats
        dims = map(x->parse(Float64,x), dims)
        VX[i] = dims[2]
        VY[i] = dims[3]
        push!(vertices, (dims[2], dims[3]))
    end
    #define EtoV matrix
    EtoV = zeros(Int, K, 3)

    # the lines with data start at lines[11+Nv]
    for k ∈ 1:K
        data = lines[11+Nv+k]
        #split up the blank spaces
        dims = split(data)
        #reinterpret the strings as Ints
        dims = map(x->parse(Int,x), dims)
        #set the elements
        EtoV[k,:] = dims[4:6]
    end

    #close the file
    close(f)
    return Nv, VX, VY, K, EtoV
end

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

# Return Values:

-   `ℳ`: a Mesh2D object with vertices as specified

# Example

"""
function rectmesh2D(xmin, xmax, ymin, ymax, K, L)
    # 1D arrays
    vˣ,mapˣ = unimesh1D(xmin, xmax, K)
    vʸ,mapʸ = unimesh1D(ymin, ymax, L)

    # construct array of vertices
    vertices = Tuple{Float64,Float64}[]
    for x in vˣ
        for y in vʸ
            push!(vertices, (x,y))
        end
    end

    # construct element to vertex map
    EtoV = Int.(ones(K*L, 4))
    j = 0
    for k in 1:K
        for l in 1:L
            j += 1

            EtoV[j,2] = Int(l + (L+1) * (k-1))
            EtoV[j,3] = Int(l + (L+1) * k)

            EtoV[j,1] = Int(EtoV[j,2] + 1)
            EtoV[j,4] = Int(EtoV[j,3] + 1)
        end
    end

    ℳ = Mesh2D(j, vertices, EtoV, 4)

    return ℳ
end


"""
makenodes2D(ℳ::Mesh2D, N::Int)

# Description

    construct the elements and GL points for a 2D mesh

# Arguments

-   `ℳ`: a Mesh2D object to fill in with GL points
-   `N`: polynomial order along each face, is the same for each element for now

# Return Values:

-   `Ω`: array of all the elements in the grid
-   `x`: physical coordinates of GL points

"""
function makenodes2D(ℳ::Mesh2D, N::Int)
    Ω = Element2D[]
    x = Float64[]
    for k in 1:ℳ.K
        # check number of faces, maybe eventually do on an element by element basis
        if ℳ.nFaces == 3
            # build a triangle, needs to be implemented
            return
        elseif ℳ.nFaces == 4
            # build a rectangle
            Ωᵏ = rectangle(k, ℳ.EtoV, N, N, ℳ.vertices)
        else
            # we will never use any other polygon ???
            return
        end

        # fill in arrays
        push!(Ω, Ωᵏ)
        x = cat(x, Ωᵏ.x; dims=1)
    end

    return Ω,x
end

"""
connect2D(EToV)

# Description

    Build connectivity maps for an arbitrary Mesh2D
    (currently assumes all elements have same number of faces)

# Arguments

-   `ℳ`: a Mesh2D object for which to create the maps

# Return Values

-   `EToE`: element to element map
-   `EToF`: element to face map

"""
function connect2D(ℳ::Mesh2D)
    total_faces = ℳ.nFaces * ℳ.K

    # list of local face to local vertex connections
    vn = zeros(Int, ℳ.nFaces, 2)
    for i in 1:ℳ.nFaces
        j = i % ℳ.nFaces + 1
        vn[i,:] = [i,j]
    end

    # build global face to node sparse array
    FtoV = spzeros(Int, total_faces, maximum(ℳ.EtoV))
    let sk = 1
        for k in 1:ℳ.K
            for face in 1:ℳ.nFaces
                @. FtoV[sk, ℳ.EtoV[k, vn[face,:] ] ] = 1;
                sk += 1
            end
        end
    end

    # global face to global face sparse array
    FtoF = FtoV * FtoV' - 2I # gotta love julia

    #find complete face to face connections
    faces1, faces2 = findnz(FtoF .== 2)

    # convert face global number to element and face numbers
    element1 = @. floor(Int, (faces1 - 1) / ℳ.nFaces ) + 1
    element2 = @. floor(Int, (faces2 - 1) / ℳ.nFaces ) + 1

    face1 = @. mod((faces1 - 1) , ℳ.nFaces ) + 1
    face2 = @. mod((faces2 - 1) , ℳ.nFaces ) + 1

    # Rearrange into Nelement x Nfaces sized arrays
    ind = diag( LinearIndices(ones(Int, ℳ.K, ℳ.nFaces))[element1,face1] ) # this line is a terrible idea.
    EtoE = collect(Int, 1:ℳ.K) * ones(Int, 1, ℳ.nFaces)
    EtoF = ones(Int, ℳ.K, 1) * collect(Int, 1:ℳ.nFaces)'
    EtoE[ind] = copy(element2)
    EtoF[ind] = copy(face2)

    return EtoE,EtoF
end

"""
buildmaps2D(ℳ::Mesh2D, _nFP::Int, _nGL::Int, _fmask, _nodes)

# Description

    Build connectivity matrices for computing fluxes

# Arguments

-   `ℳ`: a Mesh2D object to compute the maps for
-   `_nFP`: number of points along each face, currently fixed for each face and element
-   `_nGL`: number of GL points in an element, currently fixed
-   `_fmask`: matrix of indices of GL points along each face, currently fixed for each element
-   `_nodes`: the complete list of physical GL points on the grid

# Return Values

-   `vmap⁻`: vertex indices, (used for interior u values)
-   `vmap⁺`: vertex indices, (used for exterior u values)
-   `vmapᴮ`: vertex indices, corresponding to boundaries
-   `mapᴮ`: use to extract vmapᴮ from vmap⁻

"""
function buildmaps2D(ℳ::Mesh2D, _nFP::Int, _nGL::Int, _fmask, _nodes)
    # create face to node connectivity matrix
    EtoE,EtoF = connect2D(ℳ)

    # number volume nodes consecutively
    nodeids = reshape( collect(Int, 1:(ℳ.K * _nGL)), _nGL, ℳ.K)
    vmap⁻ = zeros(Int, _nFP, ℳ.nFaces, ℳ.K)
    vmap⁺ = zeros(Int, _nFP, ℳ.nFaces, ℳ.K)
    # not actually used ???
    # map⁻  = collect(Int, 1:(_nFP * ℳ.nFaces * ℳ.K))'
    # map⁺  = copy(reshape(map⁻, _nFP, ℳ.nFaces, ℳ.K))

    # find index of interior face nodes wrt volume node ordering
    for k in 1:ℳ.K
        for f in 1:ℳ.nFaces
            vmap⁻[:, f, k] = nodeids[_fmask[:, f], k]
        end
    end

    # find indices of exterior face nodes that match interior face nodes
    let one = ones(1, _nFP)
        for k1 in 1:ℳ.K
            for f1 in 1:ℳ.nFaces
                # find neighbor
                k2 = EtoE[k1, f1]
                f2 = EtoF[k1, f1]

                # find volume node numbers of interior and exterior nodes
                vid⁻ = vmap⁻[:, f1, k1]
                x⁻ = _nodes[:, 1][vid⁻]
                y⁻ = _nodes[:, 2][vid⁻]

                vid⁺ = vmap⁻[:, f2, k2]
                x⁺ = _nodes[:, 1][vid⁺]
                y⁺ = _nodes[:, 2][vid⁺]

                # create distance matrix
                D = zeros(length(x⁻), _nFP)
                for i in 1:_nFP
                    for j in 1:i
                        D[j,i] = (x⁻[i] - x⁺[j])^2 + (y⁻[i] - y⁺[j])^2
                    end
                end
                D = Symmetric(D)

                # reference length of edge
                v1 = ℳ.vertices[ℳ.EtoV[k1,f1]]
                v2 = ℳ.vertices[ℳ.EtoV[k1, 1 + mod(f1, ℳ.nFaces)]]
                refd = @. sqrt((v1[1] - v2[1])^2 + (v1[2] - v2[2])^2)

                # find indices of GL points on the boundary of the element
                # i.e. distance between them is less than the reference difference
                mask = @. sqrt(abs(D)) < eps(refd)

                # convert from matrix mask to linear mask
                m,n = size(D)
                d = collect(Int, 1:(m*n))[mask[:]]

                # find IDs of matching interior and exterior GL nodes
                id⁻ =  @. floor(Int, (d-1)/m) + 1
                id⁺ =  @. mod(d-1, m) + 1

                # save exterior node that interior node maps to
                vmap⁺[id⁻, f1, k1] = vid⁺[id⁺]
                # @. map⁺[id⁻, f1, k1] = id⁺ + (f2-1) * _nFP + (k2-1) * ℳ.nFaces * _nFP
            end
        end
    end

    # reshape arrays
    vmap⁺ = reshape(vmap⁺, length(vmap⁺))
    vmap⁻ = reshape(vmap⁻, length(vmap⁻))

    # Create list of boundary nodes
    mapᴮ = collect(Int, 1:length(vmap⁺))[vmap⁺ .== vmap⁻]
    vmapᴮ = vmap⁻[mapᴮ]

    return vmap⁻, vmap⁺, vmapᴮ, mapᴮ
end
