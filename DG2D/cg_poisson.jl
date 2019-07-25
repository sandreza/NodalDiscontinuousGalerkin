#Note that nothing in this file is correct


"""
cg_Δ!(Δu, u, params, t)


# Description

- computes the second derivative

# Arguments

- `Δu` : the laplacian of u
- `u` :  the thing we want to take laplacian of
- `ι` : struct for temporary variables
- `mesh` : the mesh struct with all the grid information


"""
function cg_Δ!(Δu, u, ι, mesh)
    # computing gradient
    ∇!(ι.φˣ, ι.φʸ, u, mesh)
    # compute divergence of flux, gets the Laplacian
    ∇⨀!(Δu, ι.φˣ, ι.φʸ, mesh)
    return nothing
end

# need this to enforce boundary conditions, faster than setdiff!
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end


#builds the matrix (one column at a time) (sparse matrix)
function cg_poisson_setup(ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
    L = spzeros(length(mesh.x), length(mesh.x))
    @. ι.u = 0.0
    Δq = copy(ι.u)
    q =  copy(ι.u)
    mapC = copy(mesh.vmapM)
    mapD = copy(mesh.vmapP)
    mapB = copy(mesh.vmapB)
    @. q = 0
    # first construct second derivative matrix
    for j in 1:length(mesh.x)
        q[j] = 1.0
        dg_poisson!(Δq, q, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
        q[j] = 0.0
        L[:,j] =  Δq[:]
    end
    dropzeros!(L)
    # partition boundary condition nodes to enforce continuity and derivative
    for j in 1:length(mesh.vmapM)
        if mesh.vmapM[j] ∈ mapC
            # remove both the node and its neighbor
            remove!(mapC, mesh.vmapM[j])
            remove!(mapC, mesh.vmapP[j])
            # remove only the node
            remove!(mapD, mesh.vmapM[j])
        end
    end
    # find nodes to enforce continuity
    mapC = setdiff(mesh.vmapM, mapD)
    # find nodes to enforce derivative, get rid of duplicates
    mapD = union(mapD)
    # enforce continuity across nodes, modify row by row
    for j in mapC
        @. L[j,:] = 0
        L[j,j] = 1.0
        L[j, mesh.vmapP[j]] = -1.0
    end
    dropzeros!(L)
    for j in mapD
        @. L[j,:] = 0
        L[j,j] = normal_derivative[j]
        L[j, mesh.vmapM[j]] = -normal_derivative[vmapP[j]]
    end
    # enfore boundary condtions
    LL = sparse( mesh.J .* (mesh.M * L) )
    return LL, mapC, mapD
end
