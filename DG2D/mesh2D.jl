include("../DG1D/mesh.jl")


"""
meshreader_gambit2D(filename)

# Description

- Reads the .neu files from NDG

# Arguments

- `filename`: a string that contains the path to the file

# Return : Nv, VX, VY, K, EToV

- `Nv` :   number of vertices
- `VX` :   x coordinate of vertex
- `VY` :   y coordinate of vertex
- `K`  :   number of elements
- `EToV` : element to vertex connection

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
    #define EToV matrix
    EToV = zeros(Int, K, 3)

    # the lines with data start at lines[11+Nv]
    for k ∈ 1:K
        data = lines[11+Nv+k]
        #split up the blank spaces
        dims = split(data)
        #reinterpret the strings as Ints
        dims = map(x->parse(Int,x), dims)
        #set the elements
        EToV[k,:] = dims[4:6]
    end

    #close the file
    close(f)
    return Nv, VX, VY, K, EToV
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
