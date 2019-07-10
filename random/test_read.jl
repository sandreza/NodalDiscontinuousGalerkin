# This file here is mostly for debugging

#n = 10
#n = 3
n = 3
#FileName = "Maxwell025.neu"
FileName = "Maxwell2.neu"
filepath = "./DG2D/grids/"


include("../utils.jl")
include("../DG2D/triangles.jl")

filename = filepath*FileName
test1 = garbage_triangle2(n, filename)
#open file
f = open(filepath*FileName)
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

#EToE, EToF = connect2D(EToV)
EToE, EToF = triangle_connect2D(EToV)

vmapM, vmapP, vmapB, mapB = buildmaps2D(K, np, nfp, nfaces, fmask, EToE, EToF, EToV, x, y, VX, VY)
Vr, Vs = dvandermonde2D(n,r,s)
Drw = (V*Vr') / (V*V')
Dsw = (V*Vs') / (V*V')



p1 = scatter(VX, VY, legend=false)
p2 = scatter(x, y, legend=false)
display(plot(p1,p2))
