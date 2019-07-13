# This file here is mostly for debugging
#=
Some general comments:

EToV matrix : element to vertex matrix, each row corresponds to an element, each column, j,  corresponds to the vertex (VX[j], VY[j])
=#
using Plots
#n = 10
#n = 3
n = 5
periodic = false
#FileName = "Maxwell025.neu"
FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"


include("../src/utils.jl")
include("../DG2D/triangles.jl")


filename = filepath*FileName
#test1 = garbage_triangle3(n, filename)
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

if periodic

    ax = minimum(VX)
    bx = maximum(VX)
    ay = minimum(VY)
    by = maximum(VY)
    xperiod = bx - ax
    yperiod = by - ay
    # build periodic index converter
    conv = collect(1:length(VX))
    #build association map to create EToVp (periodic version)
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

    # relabel indices
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

    EToVp = copy(EToV)
    for i in 1:length(EToVp)
        EToVp[i] = conv[ EToVp[i] ] # convert to appropriate vertex
    end
    scatter(VX,VY)
    scatter!(VX[conv],VY[conv],legend=false)

end

p1 = scatter(VX, VY, legend=false)
display(p1)


# Basically gonna be the same as startup2D
#n = 8

nfp = n+1
np = Int( (n+1)*(n+2)/2 )
nfaces = 3
nodetol = 1e-12

#compute nodal set
x, y = nodes2D(n)
r, s = xytors(x,y)

#plot them
theme(:default)
p1 = scatter(x,y, legend = false)
p2 = scatter(r,s, legend = false)
display(plot(p1,p2))

#build reference elements and matrices
V = vandermonde2D(n, r, s)
invV = inv(V)
mass_matrix = invV' * invV
Dr, Ds = dmatrices2D(n , r, s, V)

x,y = global_grid(r, s, EToV, VX, VY)

fmask = create_fmask(r, s)
edge_x, edge_y = find_edge_nodes(fmask, x, y)
lift = lift_tri(n, fmask, r, s, V)

rx, sx, ry, sy, J = geometricfactors2D(x, y, Dr, Ds)

nx, ny, sJ = normals2D(x, y, Dr, Ds, fmask, nfp, K)
Fscale = sJ ./ J[fmask[:],:]

EToE, EToF = connect2D(EToV)

EToEp, EToFp = connect2D(EToVp) #periodic map

EToE, EToF = triangle_connect2D(EToV)

vmapM, vmapP, vmapB, mapB = buildmaps2D(K, np, nfp, nfaces, fmask, EToE, EToF, EToV, x, y, VX, VY)
vmapMp, vmapPp, vmapBp, mapBp = build_periodic_maps2D(K, np, nfp, nfaces, fmask, EToEp, EToFp, EToVp, x, y, VX, VY)

Vr, Vs = dvandermonde2D(n,r,s)
Drw = (V*Vr') / (V*V')
Dsw = (V*Vs') / (V*V')



p1 = scatter(VX, VY, legend=false)
p2 = scatter(x, y, legend=false)
display(plot(p1,p2))


scatter(x[vmapMp[1:4]],y[vmapMp[1:4]])
scatter!(x[vmapPp[1:4]],y[vmapPp[1:4]])
