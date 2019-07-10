include("../src/utils.jl")


using SparseArrays
using LinearAlgebra

abstract type AbstractGrid2D end
abstract type AbstractMesh2D <: AbstractGrid2D end

struct Grid2D{S, T, U} <: AbstractGrid2D
    K::S
    vertices::T
    EtoV::U
    nFaces::S

    function Grid2D(K, vertices, EtoV, nFaces)
        return new{typeof(K),typeof(vertices),typeof(EtoV)}(K, vertices, EtoV, nFaces)
    end
end



struct Mesh2D{S, T, U, V} <: AbstractMesh2D
    # basic grid of vertices
    grid::S

    # elements
    elements::T

    # GL points
    nodes::U

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
function meshreader_gambit2D(filename)
    #open file
    f = open(filename)
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
    vertices = Tuple{Float64,Float64}[]
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
    return Nv, VX, VY, vertices, K, EtoV
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
    vertices = Tuple{Float64,Float64}[]
    for x in vx
        for y in vy
            push!(vertices, (x,y))
        end
    end
    # v = reshape(v, K+1, L+1)

    # construct element to vertex map
    EtoV = Int.(ones(K*L, 4))
    j = 0
    for l in 1:L
        for k in 1:K
            j += 1

            EtoV[j,2] = Int(k + (L+1) * (l-1))
            EtoV[j,3] = Int(k + (L+1) * l)

            EtoV[j,1] = Int(EtoV[j,2] + 1)
            EtoV[j,4] = Int(EtoV[j,3] + 1)
        end
    end

    grid = Grid2D(j, vertices, EtoV, 4)

    return grid
end


"""
makenodes2D()
"""
function makenodes2D(𝒢::Grid2D, N::Int)
    Ω = Element2D[]
    for k in 1:𝒢.K
        # check number of faces, maybe eventually do on an element by element basis
        if 𝒢.nFaces == 3
            # build a triangle
            return
        elseif 𝒢.nFaces == 4
            Ωᵏ = rectangle(k, 𝒢.EtoV, N, N, 𝒢.vertices)
        else
            return
        end

        push!(Ω, Ωᵏ)
    end

    return Ω
end


"""
buildmaps2D(EtoV)

# Description

- Make element-to-element and element-to-face maps for a 2D grid

# Arguments

-   `𝒢`: a 2D grid object
-   `x,y`: physical coordinates of all the GL points for the mesh
-

# Output

# Comments

- The changes from the 1D are minor. nFaces can probably remain generic

"""
function buildmaps2D(𝒢::Grid2D, nFP::Int, nodes, fmask)
    # create face to node connectivity matrix
    total_faces = 𝒢.nFaces * 𝒢.K

    # list of local face to local vertex connections
    vn = zeros(Int, 𝒢.nFaces, 2)
    for i in 1:𝒢.nFaces
        j = i % 𝒢.nFaces + 1
        vn[i, :] = [i,j]
    end

    # build global face to node sparse array
    FtoV = spzeros(Int, total_faces, maximum(𝒢.EtoV))
    let sk = 1
        for k in 1:𝒢.K
            for face in 1:𝒢.nFaces
                @. FtoV[sk, 𝒢.EtoV[k, vn[face,:] ] ] = 1;
                sk += 1
            end
        end
    end

    # global face to global face sparse array
    FtoF = FtoV * FtoV' - 2I #gotta love julia

    #find complete face to face connections
    faces1, faces2 = findnz(FtoF .== 2)

    # convert face global number to element and face numbers
    element1 = @. floor(Int, (faces1 - 1) / 𝒢.nFaces ) + 1
    element2 = @. floor(Int, (faces2 - 1) / 𝒢.nFaces ) + 1

    face1 = @. mod((faces1 - 1) , nFaces ) + 1
    face2 = @. mod((faces2 - 1) , nFaces ) + 1

    # Rearrange into Nelement x Nfaces sized arrays
    ind = diag( LinearIndices(ones(Int, 𝒢.K, 𝒢.nFaces))[element1,face1] ) # this line is a terrible idea.
    EtoE = collect(Int, 1:𝒢.K) * ones(Int, 1, 𝒢.nFaces)
    EtoF = ones(Int, 𝒢.K, 1) * collect(Int, 1:𝒢.nFaces)'
    EtoE[ind] = copy(element2)
    EtoF[ind] = copy(face2)
    ### end connect2D from book

    ### start buildmaps2D from book
    # number volume nodes consecutively
    nodeids = reshape( collect(Int, 1:(𝒢.K * length(nodes))), length(nodes), 𝒢.K)
    vmapM = zeros(Int, nFP, 𝒢.nFaces, 𝒢.K)
    vmapP = zeros(Int, nFP, 𝒢.nFaces, 𝒢.K)
    mapM  = collect(Int, 1:(nFP * 𝒢.nFaces * 𝒢.K))'
    mapP  = copy(reshape(mapM, nFP, 𝒢.nFaces, 𝒢.K))

    # find index of face nodes wrt volume node ordering
    for k in 1:𝒢.K
        for f in 1:𝒢.nFaces
            vmapM[:, f, k] = nodeids[fmask[:, f], k]
        end
    end

    let one = ones(1, nFP)
        for k1 in 1:𝒢.K
            for f1 in 1:𝒢.nFaces
                # find neighbor
                k2 = EtoE[k1, f1]
                f2 = EtoF[k1, f1]

                # reference length of edge
                v1 = 𝒢.EtoV[k1,f1]
                v2 = 𝒢.EtoV[k1, 1 + mod(f1, 𝒢.nFaces)]
                refd = @. sqrt((𝒢.vertices[v1][1] - 𝒢.vertices[v2][1])^2 + (𝒢.vertices[v1][2] - 𝒢.vertices[v2][2])^2)

                # find volume node numbers of left and right nodes
                vidM = vmapM[:, f1, k1]
                x⁻ = nodes[vidM]

                vidP = vmapM[:, f2, k2]
                x⁺ = nodes[vidP]

                # create distance matrix
                D = Symmetric(zeros(length(x⁻), nFP))
                for i in 1:nFP
                    for j in 1:i
                        D[i,j] = (x⁻[i][1] - x⁺[j][1])^2 + (x⁻[i][2] - x⁺[j][2])^2
                    end
                end

                mask = @. D < eps(refd)

                # find linear indices
                m,n = size(D)
                d = collect(Int, 1:(m*n))
                idM =  @. floor(Int, (d[mask[:]]-1) / m) + 1
                idP =  @. mod( (d[mask[:]]-1), m) + 1
                vmapP[idM, f1, k1] = vidP[idP]
                @. mapP[idM, f1, k1] = idP + (f2-1) * nFP + (k2-1) * 𝒢.nFaces * nFP
            end
        end
    end

    # reshape arrays
    vmapP = reshape(vmapP, length(vmapP))
    vmapM = reshape(vmapM, length(vmapM))

    # Create list of boundary nodes
    mapB = collect(Int, 1:length(vmapP))[vmapP .== vmapM]
    vmapB = vmapM[mapB]

    return vmapM, vmapP, vmapB, mapB
end
