include("../src/utils.jl")

using SparseArrays

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
    println(lines[6])
    println(lines[7])
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
            EtoV[j,2] = Int(k + (L+1) * (l-1))
            EtoV[j,3] = Int(k + (L+1) * l)

            EtoV[j,1] = Int(EtoV[j,2] + 1)
            EtoV[j,4] = Int(EtoV[j,3] + 1)

            j += 1
        end
    end

    return vertices,EtoV
end

"""
connect2D(EtoV)

# Description

- Make element-to-element and element-to-face maps for a 2D grid

# Arguments

-   `nfaces`: number of faces on each element (constant across the grid for now)
-   `EtoV`: element to vertices map

# Output

-   `EtoE`: element to element map
-   `EtoF`: element to face map

# Comments

- The changes from the 1D are minor. nfaces can probably remain generic

"""
function connect2D(nfaces::Int, EtoV)
    #find number of elements and vertices
    K = size(EtoV, 1)
    Nv = maximum(EtoV)

    # create face to node connectivity matrix
    total_faces = nfaces * K

    # list of local face to local vertex connections
    vn = zeros(Int, nfaces, 2)
    for i in 1:nfaces
        j = i % nfaces + 1
        vn[i, :] = [i,j]
    end

    # build global face to node sparse array
    FtoV = spzeros(Int, total_faces, Nv)
    let sk = 1
        for k in 1:K
            for face in 1:nfaces
                @. FtoV[sk, EtoV[k, vn[face,:] ] ] = 1;
                sk += 1
            end
        end
    end

    # global face to global face sparse array
    FtoF = FtoV * FtoV' - 2I #gotta love julia

    #find complete face to face connections
    faces1, faces2 = findnz(FtoF .== 2)

    # convert face global number to element and face numbers
    element1 = @. floor(Int, (faces1 - 1) / nfaces ) + 1
    element2 = @. floor(Int, (faces2 - 1) / nfaces ) + 1

    face1 = @. mod((faces1 - 1) , nfaces ) + 1
    face2 = @. mod((faces2 - 1) , nfaces ) + 1

    # Rearrange into Nelement x Nfaces sized arrays
    ind = diag( LinearIndices(ones(Int, K, nfaces))[element1,face1] ) # this line is a terrible idea.
    EtoE = collect(Int, 1:K) * ones(Int, 1, nfaces)
    EtoF = ones(Int, K,1) * collect(Int, 1:nfaces)'
    EtoE[ind] = copy(element2);
    EtoF[ind] = copy(face2);
    return EtoE, EtoF
end
