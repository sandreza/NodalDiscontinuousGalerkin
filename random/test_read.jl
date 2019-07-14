# This file here is mostly for debugging
#=
Some general comments:

EtoV matrix : element to vertex matrix, each row corresponds to an element, each column, j,  corresponds to the vertex (VX[j], VY[j])
=#
using Plots
#n = 10
#n = 3
n = 3
periodic = false
#FileName = "Maxwell025.neu"
FileName = "Maxwell2.neu"
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

if periodic
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

    leftface = sort([VY[minindx] minindx], dims = 1)
    rightface = sort([VY[maxindx] maxindx], dims = 1)
    bottomface = sort([VX[minindy] minindy], dims = 1)
    topface = sort([VX[maxindy] maxindy], dims = 1)

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

    EtoVp = copy(EtoV)
    for i in 1:length(EtoVp)
        EtoVp[i] = conv[ EtoVp[i] ] # convert to appropriate vertex
    end
    scatter(VX,VY)
    scatter!(VX[conv],VY[conv],legend=false)

    for i in 1:length(VX)
        println("vertex ")
        println((VX[i],VY[i]))
        println("gets converted to")
        println((VX[conv[i]],VY[conv[i]]))
    end
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

x,y = global_grid(r, s, EtoV, VX, VY)

fmask = create_fmask(r, s)
edge_x, edge_y = find_edge_nodes(fmask, x, y)
lift = lift_tri(n, fmask, r, s, V)

rx, sx, ry, sy, J = geometricfactors2D(x, y, Dr, Ds)

nx, ny, sJ = normals2D(x, y, Dr, Ds, fmask, nfp, K)
Fscale = sJ ./ J[fmask[:],:]

EtoE, EtoF = connect2D(EtoV)

EtoVp = make_periodic_2D(VX,VY, EtoV)

EtoEp, EtoFp = connect2D(EtoVp) #periodic map
EtoEp, EtoFp = connect_periodic_2D(VX, VY, EtoV) #periodic map

EtoE, EtoF = triangle_connect2D(EtoV)

vmapM, vmapP, vmapB, mapB = buildmaps2D(K, np, nfp, nfaces, fmask, EtoE, EtoF, EtoV, x, y, VX, VY)

vmapMp, vmapPp, vmapBp, mapBp = build_periodic_maps2D(K, np, nfp, nfaces, fmask, EtoEp, EtoFp, EtoV, x, y, VX, VY)

Vr, Vs = dvandermonde2D(n,r,s)
Drw = (V*Vr') / (V*V')
Dsw = (V*Vs') / (V*V')



p1 = scatter(VX, VY, legend=false)
p2 = scatter(x, y, legend=false)
display(plot(p1,p2))

#=
scatter(x[vmapMp[1:4]],y[vmapMp[1:4]])
scatter!(x[vmapPp[1:4]],y[vmapPp[1:4]])
=#
sum( vmapPp .== 0)

getpair(i) = (VX[i], VY[i])


#=
 e1 = @. mod(vmapB-1, np) + 1
 e2 = @. div(vmapB-1, np) + 1

j = 4
ind1 = 1 + (j-1)*nfp
ind2 = j * nfp
scatter!(x[vmapB[ind1:ind2]],y[vmapB[ind1:ind2]])
=#

###
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

leftface = sortslices([VY[minindx] minindx], dims = 1)
rightface = sortslices([VY[maxindx] maxindx], dims = 1)
bottomface = sortslices([VX[minindy] minindy], dims = 1)
topface = sortslices([VX[maxindy] maxindy], dims = 1)

nFaces = 3
K = size(EtoV, 1)

total_faces = nFaces * K

# list of local face to local vertex connections
# face convention implied here
# face 1 is associated with vn[1,2]
# face 2 is associated with vn[2,3]
# ...
# face nfaces is associated with vn[nfaces, 1]
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

for i in 1:(nfaces*K)
    for k in 1:(length(leftface[:,2])-1)
        vecL = Int.(leftface[k:k+1 , 2])
        vecR = Int.(rightface[k:k+1, 2])
        # identify left face with right face
        if sum(FtoV[i, vecL])==2
            @. FtoV[i, vecL] = 0
            @. FtoV[i,vecR] = 1
            dropzeros!(FtoV)
        end
        vecB = Int.(bottomface[k:k+1 , 2])
        vecT = Int.(topface[k:k+1, 2])
        # identify top face with bottom face
        if sum(FtoV[i, vecB])==2
            @. FtoV[i, vecB] = 0
            @. FtoV[i,vecT] = 1
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

###
