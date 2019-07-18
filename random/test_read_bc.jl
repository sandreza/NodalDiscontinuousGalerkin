using Plots
#n = 10
#n = 3
n = 8
periodic = false
#FileName = "Maxwell025.neu"

FileName = "pvortex4A01.neu"
#FileName = "pvortex2A025.neu"
#FileName = "pvortex3A025.neu"
#FileName = "pvortex4A01.neu"
filepath = "./DG2D/grids/"

include("../src/utils.jl")
include("../DG2D/triangles.jl")

filename = filepath*FileName
test1 = garbage_triangle3(n, filename)
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

bc_line = []
bc_label = []
bc_name = ["In", "Out", "Wall", "Cyl"]
for j in 1:length(lines)
    for bc_type ∈ bc_name
        if  occursin(bc_type, lines[j])
            push!(bc_line, j)
            push!(bc_label, bc_type)
        end
    end
end

#close the file
bctype = zeros(Int,size(EtoV))
for j in bc_line
    bc_type_index = findall(j .== bc_line)
    for i in 1:(length(lines)-j)
        if occursin("ENDOFSECTION", lines[j+i])
            break
        end
        bcs = split(lines[j+i])
        display(bcs)
        bcs = map(x->parse(Int,x), bcs)
        bctype[bcs[1], bcs[3]] = bc_type_index[1]
    end
end


leftface = findall(VX .== minimum(VX))
rightface = findall(VX .== maximum(VX))
bottomface = findall(VY .== minimum(VY))
topface = findall(VY .== maximum(VY))
boundary = union(leftface,rightface, topface, bottomface)
bc1 = [EtoV[i] in boundary for i in 1:length(EtoV)]
boundary_element = [EtoV[i] == true for i in 1:length(EtoV[:,1])]


close(f)

mesh = test1

plot_mesh(mesh)
Nv, VX, VY, K, EtoV, bctype, bc_name = meshreader_gambit_bc_2D(filename)

mapT, vmapT, newbc_label = build_bc_maps(mesh, bctype, bc_name)


color_list = ["red", "blue"]
p = []
for j in 1:length(mapT)
    p1 = scatter(mesh.x[vmapT[j]], mesh.y[vmapT[j]], color = color_list[j], title = newbc_label[j], legend = false)
    push!(p,p1)
end
display(plot(p...))
