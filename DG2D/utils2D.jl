"""
∇(u, Ωᵏ)

# Description

    Compute gradient of u wrt physical grid

# Arguments

-   `u`: scalar to take gradient of
-   `Ωᵏ`: element to compute in

# Return Values

-   `uˣ`: first component of the gradient
-   `uʸ`: second component of the gradient

"""
function ∇(u, Ωᵏ)
    # compute partial derivatives on ideal grid
    uʳ = Ωᵏ.Dʳ * u
    uˢ = Ωᵏ.Dˢ * u

    # compute partial derivatives on physical grid
    uˣ = @. Ωᵏ.rˣ * uʳ + Ωᵏ.sˣ * uˢ
    uʸ = @. Ωᵏ.rʸ * uʳ + Ωᵏ.sʸ * uˢ

    return uˣ,uʸ
end

"""
∇⨀(x, y, Ωᵏ)

# Description

    Compute the divergence of u=(x,y) wrt physical grid

# Arguments

-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ωᵏ`: element to compute in

# Return Values

-   `∇⨀u`: the divergence of u

"""
function ∇⨀(x, y, Ωᵏ)
    # compute partial derivatives on ideal grid
    xʳ = Ωᵏ.Dʳ * x
    xˢ = Ωᵏ.Dˢ * x
    yʳ = Ωᵏ.Dʳ * y
    yˢ = Ωᵏ.Dˢ * y

    # compute gradient on physical grid
    ∇⨀u = @. Ωᵏ.rˣ * xʳ + Ωᵏ.sˣ * xˢ + Ωᵏ.rʸ * yʳ + Ωᵏ.sʸ * yˢ

    return ∇⨀u
end

"""
∇⨂(x, y, Ωᵏ)

# Description

    Compute the curl of u=(x,y) wrt physical grid

# Arguments

-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ωᵏ`: element to compute in

# Return Values

-   `∇⨂u`: the curl of u

"""
function ∇⨂(x, y, Ωᵏ)
    # compute partial derivatives on ideal grid
    xʳ = Ωᵏ.Dʳ * x
    xˢ = Ωᵏ.Dˢ * x
    yʳ = Ωᵏ.Dʳ * y
    yˢ = Ωᵏ.Dˢ * y

    # compute gradient on physical grid
    ∇⨂u = @. Ωᵏ.rˣ * yʳ + Ωᵏ.sˣ * yˢ - Ωᵏ.rʸ * xʳ - Ωᵏ.sʸ * xˢ

    return ∇⨂u
end


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
