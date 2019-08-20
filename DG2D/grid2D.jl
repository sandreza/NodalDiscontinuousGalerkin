include("../src/utils.jl")
include("rectangles.jl")
include("triangles.jl")

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
struct Grid2D{S, T, U, V} <: AbstractGrid2D
    # basic grid of vertices
    ℳ::S

    # elements
    Ω::T

    # GL points
    x::U # physical coordinates
    nGL::V # number of points

    function Grid2D(ℳ::Mesh2D, N::Int; periodic::Bool=false)
        # get elements and GL nodes
        Ω,x̃ = makenodes2D(ℳ, N)
        nGL,_ = size(x̃)
        # make a facemask

        # build the boundary maps
        l̃ˣ = abs(maximum(x̃[:,1]) - minimum(x̃[:,1]))
        l̃ʸ = abs(maximum(x̃[:,2]) - minimum(x̃[:,2]))
        buildmaps2D(ℳ, Ω, nGL; lˣ=l̃ˣ, lʸ=l̃ʸ, periodic=periodic)

        return new{typeof(ℳ),typeof(Ω),typeof(x̃),typeof(nGL)}(ℳ,Ω, x̃,nGL)
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
    vertices = zeros(Nv, 2)

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
        vertices[i,:] = [dims[2] dims[3]]
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

    # close the file
    close(f)

    # make mesh object
    ℳ = Mesh2D(K, vertices, EtoV, 3)

    return ℳ
end

"""
meshreader_gambit_bc_2D(filename)

# Description

- Reads the .neu files from NDG, includes bc data

# Arguments

- `filename`: a string that contains the path to the file

# Return : Nv, VX, VY, K, EtoV

- `Nv` :   number of vertices
- `VX` :   x coordinate of vertex
- `VY` :   y coordinate of vertex
- `K`  :   number of elements
- `EtoV` : element to vertex connection
- `bctype`: matrix that contains which vertices for boundary data
- `bc_name`: for knowing what the indices actually mean

# more detail
if bctype[element_number, vertex] = index
then it is a bc_name[index] bundary condition

"""
function meshreader_gambit_bc_2D(_filename)
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
    vertices = zeros(Nv, 2)

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
        vertices[i,:] = [dims[2] dims[3]]
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

    #first find the line numbers
    bc_name = ["In", "Out", "Wall", "Cyl"] #not exhaustive
    bc_line = []
    bc_label = []
    for j in 1:length(lines)
        for bc_type ∈ bc_name
            if  occursin(bc_type, lines[j])
                push!(bc_line, j)
                push!(bc_label, bc_type)
            end
        end
    end

    #then find the data for the line numbers
    bctype = zeros(Int, size(EtoV))
    for j in bc_line
        bc_type_index = findall(j .== bc_line)
        for i in 1:(length(lines)-j)
            if occursin("ENDOFSECTION", lines[j+i])
                break
            end
            bcs = split(lines[j+i])
            bcs = map(x->parse(Int,x), bcs)
            bctype[bcs[1], bcs[3]] = bc_type_index[1]
        end
    end

    close(f)
    return Nv, VX, VY, K, EtoV, bctype, bc_name
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
    vertices = zeros((K+1)*(L+1), 2)
    let i = 0
        for y in vʸ
            for x in vˣ
                i += 1
                vertices[i,:] = [x y]
            end
        end
    end

    # construct element to vertex map
    EtoV = ones(Int, K*L, 4)
    j = 0
    for l in 1:L
        for k in 1:K
            j += 1

            EtoV[j,1] = Int(k + (K+1) * (l-1))
            EtoV[j,4] = Int(k + (K+1) * l)

            EtoV[j,3] = Int(EtoV[j,4] + 1)
            EtoV[j,2] = Int(EtoV[j,1] + 1)
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
        vertices = view(ℳ.EtoV, k, :)
        nFaces = length(vertices)
        if nFaces == 3
            # build a triangle
            Ωᵏ = triangle(k, vertices, N, ℳ.vertices)
        elseif nFaces == 4
            # build a rectangle
            Ωᵏ = rectangle(k, vertices, N, N, ℳ.vertices)
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
function connect2D(ℳ::Mesh2D; periodic::Bool=false)
    nFacesTotal = ℳ.nFaces * ℳ.K

    # list of local face to local vertex connections
    vn = zeros(Int, ℳ.nFaces, 2)
    for i in 1:ℳ.nFaces
        j = i % ℳ.nFaces + 1
        vn[i,:] = [i,j]
    end

    # build global face to node sparse array
    FtoV = spzeros(Int, nFacesTotal, maximum(ℳ.EtoV))
    let sk = 1
        for k in 1:ℳ.K
            for face in 1:ℳ.nFaces
                @. FtoV[sk, ℳ.EtoV[k, vn[face,:] ] ] = 1;
                sk += 1
            end
        end
    end

    # make periodic
    if periodic == true
        VX = ℳ.vertices[:,1]
        VY = ℳ.vertices[:,2]

        # find min and max values for each axis
        xmin = minimum(VX)
        ymin = minimum(VY)
        xmax = maximum(VX)
        ymax = maximum(VY)

        # find vertices on each boundary
        x⁻ = findall( VX .≈ xmin )
        y⁻ = findall( VY .≈ ymin )
        x⁺ = findall( VX .≈ xmax )
        y⁺ = findall( VY .≈ ymax )

        #match up appropriate vertices, does not generalize to 3D
        facesᴸ = sortslices([VY[x⁻] x⁻], dims = 1)[:,2]
        facesᴿ = sortslices([VY[x⁺] x⁺], dims = 1)[:,2]
        facesᴮ = sortslices([VX[y⁻] y⁻], dims = 1)[:,2]
        facesᵀ = sortslices([VX[y⁺] y⁺], dims = 1)[:,2]

        # loop over faces
        nFacesᴴ = length(facesᴸ)
        nFacesⱽ = length(facesᴮ)
        for i in 1:nFacesTotal

            # enforce periodicity in the horizontal direction
            for k in 1:(nFacesᴴ-1)
                faceᴸ = Int.(facesᴸ[k:k+1])
                faceᴿ = Int.(facesᴿ[k:k+1])

                # check if face i is faceᴸ
                if sum(FtoV[i, faceᴸ]) == 2
                    # identify left face with right face
                    @. FtoV[i, faceᴸ] = 0
                    @. FtoV[i, faceᴿ] = 1
                    dropzeros!(FtoV)
                end
            end

            # enforce periodicity in the vertical direction
            for k in 1:(nFacesⱽ-1)
                faceᴮ = Int.(facesᴮ[k:k+1])
                faceᵀ = Int.(facesᵀ[k:k+1])

                # check if face i is faceᴮ
                if sum(FtoV[i, faceᴮ]) == 2
                    # identify top face with bottom face
                    @. FtoV[i, faceᴮ] = 0
                    @. FtoV[i, faceᵀ] = 1
                    dropzeros!(FtoV)
                end
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

    connectivity = Array{Tuple{Int,Int}}[]
    for k in 1:ℳ.K
        push!(connectivity, Tuple{Int,Int}[])
        for f in 1:ℳ.nFaces
            E = EtoE[k,f]
            F = EtoF[k,f]
            push!(connectivity[k], (E,F))
        end
    end

    return connectivity
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
function buildmaps2D(ℳ::Mesh2D, Ω::Array{Element2D}, nGL::Int; lˣ=-1, lʸ=-1, periodic::Bool=false)
    # create face to node connectivity matrix
    connectivity = connect2D(ℳ, periodic=periodic)

    # find indices of interior face nodes
    nodes = collect(Int, 1:nGL)
    let nGLᵏ = 0
        for k in 1:ℳ.K
            # get element
            Ωᵏ = Ω[k]

            # get number of GL points in element k
            xᵏ = (nGLᵏ + 1):(nGLᵏ + Ωᵏ.nGL)
            nGLᵏ += Ωᵏ.nGL

            # extract indices of GL points in element k
            @. Ωᵏ.iⱽ = nodes[xᵏ]

            # save global indices of GL on face f
            for face in Ωᵏ.faces
                @. face.i⁻ = Ωᵏ.iⱽ[face.mask]
            end
        end
    end

    # find indices of exterior nodes
    for k in 1:ℳ.K
        # get element
        Ωᵏ = Ω[k]

        for h in 1:length(Ωᵏ.faces)
            # find neighbor and matching face
            j,g = connectivity[k][h]

            # get neighbor
            Ωʲ = Ω[j]

            # get faces
            f⁻ = Ωᵏ.faces[h]
            f⁺ = Ωʲ.faces[g]

            # get coordinates of GL points on the faces
            x⁻ = Ωᵏ.x[f⁻.mask, :]
            x⁺ = Ωʲ.x[f⁺.mask, :]

            # create distance matrix
            nFP⁻,_ = size(x⁻)
            nFP⁺,_ = size(x⁺)
            D = zeros(Int, nFP⁺, nFP⁻)
            for j in 1:nFP⁺
                for i in 1:nFP⁻
                    exact = (x⁻[i, 1] ≈ x⁺[j, 1]) && (x⁻[i, 2] ≈ x⁺[j, 2])
                    periodicˣ = (abs(x⁻[i, 1] - x⁺[j, 1]) ≈ lˣ) && (x⁻[i, 2] ≈ x⁺[j, 2])
                    periodicʸ = (abs(x⁻[i, 2] - x⁺[j, 2]) ≈ lʸ) && (x⁻[i, 1] ≈ x⁺[j, 1])
                    if periodic
                        D[i,j] = (exact || periodicˣ || periodicʸ)
                    else
                        D[i,j] = exact
                    end
                end
            end

            # convert from matrix mask to linear mask
            m,n = size(D)
            d = collect(Int, 1:(m*n))[Bool.(D[:])]

            # find IDs of matching interior and exterior GL nodes
            id⁻ =  @. floor(Int, (d-1)/m) + 1
            id⁺ =  @. mod(d-1, m) + 1

            # save exterior node that interior node maps to
            @. f⁻.i⁺[id⁻] = f⁺.i⁻[id⁺]

            # create list of boundary nodes
            mapᴮ = collect(Int, 1:length(f⁻.i⁻))[f⁻.i⁺ .== f⁻.i⁻]
            x = !isempty(mapᴮ)
            @. f⁻.isBoundary = [x]
        end
    end

    return nothing
end
