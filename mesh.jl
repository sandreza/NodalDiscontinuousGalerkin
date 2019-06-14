
using SparseArrays #for connectivity matrix

"""
unimesh1D(xmin, xmax, K)

#Description

    Generates a uniform 1D mesh

#Arguments

    xmin: smallest value of array

    xmax: largest values of array

    K: number of elements in an array

#Return Values: VX, EToV

    VX: vertex values | an Array of size K+1

    EToV: element to node connectivity | a Matrix of size Kx2

#Example
xmin = -1
xmax =  1
K    =  4
VX, EToV = unimesh1D(xmin, xmax, K)

"""
function unimesh1D(xmin, xmax, K)
    VX = collect(0:K) ./ K .* (xmax - xmin) .+ xmin
    EToV = Int.(ones(K, 2))
    for i = 1:K
        EToV[i,1] = Int(i)
        EToV[i,2] = Int(i+1)
    end
    return VX, EToV
end

"""
gridvalues1D(xmin, xmax, K)

#Description

    Generates physical gridpoints with each element

#Arguments

    VX: vertex values | an Array of size K+1

    EToV: element to node connectivity | a Matrix of size Kx2

    r: LGL nodes in reference element | an array

#Return Values: x

    x: physical coordinates of solution

#Example (uses dg_utils.jl as well)

xmin = 0
xmax = 2π
K = 4
#call functions
VX, EToV = unimesh1D(xmin, xmax, K)
r = jacobiGL(0, 0, 4)
x = gridvalues1D(VX, EToV, r)
#x[:,1] is the physical coordinates within the first element
#for plotting
f(x) = sin(x)
plot(x, f.(x))
#scatter(x,f.(x)) tends to work better
"""
function gridvalues1D(VX, EToV, r)
    va = view(EToV,:,1)
    vb = view(EToV,:,2)
    x = ones(length(r),1) * (VX[va]') .+ 0.5 .* (r .+ 1 ) *((VX[vb] -VX[va])')
    return x
end

"""
edgevalues1D(r, x)

#Description

    calculates edge values

#Arguments

    r: GL points

    x:  physical coordinates of solution on each element

#Return Values: x

    fx: face values of x

#Example | dg_utils.jl

r = jacobiGL(0, 0, 4)

x = gridvalues1D(VX, EToV, r)

fx = edgevalues1D(r,x)

#the locations of the edges in element 1 is fx[:, 1]


"""
function edgevalues1D(r, x)
    f1 = abs.(r .+ 1) .< eps(1.0)*10^5
    f2 = abs.(r .- 1) .< eps(1.0)*10^5
    fx1 = x[f1,:]
    fx2 = x[f2,:]
    fx = [fx1; fx2]
    return fx
end

"""
normals1D(K)

#Description

    calculates face normals

#Arguments

    K: number of elements

#Return Values: nx

    nx: face normals along each grid

#Example

"""
function normals1D(K)
    nx  = ones(2,K)
    @. nx[1,:] *= -1
    return nx
end

"""
geometric_factors(x, Dr)

#Description

    computes the geometric factors for local mappings of 1D elements

#Arguments

    x: physical coordinates of solution for each element

    Dr:

#Return Values: rx, J

    rx: inverse jacobian

    J: jacobian (in 1D a scalar)

#Example

"""
function geometric_factors(x, Dr)
    J = Dr * x
    rx = 1 ./ J #for 1D
    return rx, J
end

"""
connect1D(EToV)

#Description

    builds global connectivity arrays for 1D

#Arguments

    EToV: element to node connectivity | a Matrix of size Kx2

#Return Values: EToE, EToF

    EToE:
    EToF:

#Example

"""
function connect1D(EToV)
    nfaces = 2 #for 1d elements
    K = size(EToV,1)
    total_faces = nfaces * K
    Nv = K+1
    vn = [1, 2]
    SpFToV = Int.(spzeros(total_faces, Nv))
    let sk = 1
    for k = 1:K
        for faces = 1:nfaces
            SpFToV[ sk, EToV[k, vn[faces]]] = 1;
            sk += 1
        end
    end
    end
    SpFToF = SpFToV * (SpFToV') - sparse(I, total_faces, total_faces)
    (faces1,faces2) = findnz(SpFToF)

    element1 = @. Int(floor((faces1-1)/nfaces ) + 1)
    face1 = @. Int(mod((faces1-1),nfaces) + 1)
    element2 = @. Int(floor((faces2-1)/nfaces ) + 1)
    face2 = @. Int(mod((faces2-1),nfaces) + 1)

    #the line below is a terrible idea.
    ind = diag(LinearIndices(ones(K, nfaces))[element1,face1])
    EToE = collect(1:K) * ones(1, nfaces)
    EToF = ones(K,1) * (collect(1:nfaces)' )
    EToE[ind] = copy(element2);
    EToF[ind] = face2;
    return EToE, EToF
end

"""
buildmaps1D(K, np, nfp, nfaces, EToE, EToF, x)

#Description

    connectivity matrices for element to elements and elements to face

#Arguments
    K: number of elements
    np: number of points within an element (polynomial degree + 1)
    nfp: 1
    nfaces: 2
    fmask: an element by element mask to extract edge values
    EToE: element to element connectivity
    EToF: element to face connectivity
    x: Guass lobatto points

#Return Values: vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO

    vmapM: vertex indices, (used for interior u values)

    vmapP: vertex indices, (used for exterior u values)

    vmapB: vertex indices, corresponding to boundaries

    mapB: use to extract vmapB from vmapM

    mapI: Index of left boundary condition

    mapO: Index of right boundary condition

#Example | uses dg_utils.jl

K = 3
n = 3; α = 0; β = 0; xmin = 0; xmax = 2π;
np = n + 1
nfp = 1
nfaces = 2

r = jacobiGL(α, β, n)

VX, EToV = unimesh1D(xmin, xmax, K)
EToE, EToF = connect1D(EToV)
x = gridvalues1D(VX, EToV, r)
fx = edgevalues1D(r,x)
#build fmask
fmask1 = @. abs(r+1) < eps(1.0);
fmask2 = @. abs(r-1) < eps(1.0);
fmask  = (fmask1, fmask2)

vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO = buildmaps1D(K, np, nfp, nfaces, fmask, EToE, EToF, x)
"""
function buildmaps1D(K, np, nfp, nfaces, fmask, EToE, EToF, x)
    nodeids = reshape(collect(1:(K*np)), np, K)
    vmapM = zeros(nfp, nfaces, K)
    vmapP = zeros(nfp, nfaces, K)
    for k1 in 1:K
        for f1 in 1:nfaces
            vmapM[:,f1,k1] = nodeids[fmask[f1], k1]
        end
    end

    for k1 = 1:K
        for f1 = 1:nfaces
            k2 = Int.(EToE[k1,f1])
            f2 = Int.(EToF[k1,f1])
            vidM = Int.(vmapM[:,f1,k1])
            vidP = Int.(vmapM[:,f2,k2])
            x1 = x[vidM]
            x2 = x[vidP]

            D = (x1 .- x2) .^ 2
            if D[1] < eps(1.0)*10^5
                vmapP[:,f1,k1] = vidP
            end
        end
    end

    #reshape arrays
    vmapP = Int.( reshape(vmapP,length(vmapP)) )
    vmapM = Int.( reshape(vmapM,length(vmapM)) )

    mapB = Int.( collect(1:length(vmapP))[vmapP .== vmapM] )
    vmapB = Int.( vmapM[mapB] )

    #inflow and outflow maps
    mapI = 1
    mapO = K * nfaces
    vmapI = 1
    vmapO = K*np
    return vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO
end

#make the grid periodic
function make_periodic1D(vmapP,u)
    vmapP[1] =length(u)
    vmapP[end] = 1
    return nothing
end
