# builds the 2D helmholtz matrix
"""
solveHelmholtz!(Hu, u, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)


# Description

- Evaluate the helmholtz operator

# Arguments

- `Hu` : helholtz operator acting on u
- `u` :  the thing we want to take laplacian of
- `ι` : struct for temporary variables
- `params`: any penalty parameters that we would like to include
- `mesh` : the mesh struct with all the grid information
- `bc_u!`: function that computes boundary conditions
- `bc` : boundary condition tuple with indices
- `bc_φ!`: function that computes derivative boundary conditions
- `dbc` : boundary condition tuple with indices


"""
function solveHelmholtz!(Δu, u, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
    # unpack parameters
    τ = params[1]
    γ = params[2]
    # Form q-flux differences at faces
    @. ι.fⁿ[:] =  u[mesh.vmapM] - (u[mesh.vmapM] + u[mesh.vmapP])/2

    # Choose boundary condition type, dirichlet
    @. ι.u = u
    bc_u!(ι, mesh, bc)

    # get the numerical flux jump in the x and y directions
    @. ι.fˣ = mesh.nx * ι.fⁿ
    @. ι.fʸ = mesh.ny * ι.fⁿ

    liftx = mesh.lift * (mesh.fscale .* ι.fˣ )
    lifty = mesh.lift * (mesh.fscale .* ι.fʸ )

    # lhs of the semi-discerte PDE, ∇⋅(q) = f , q  = ∇u, qˣ = ∂ˣu, qʸ = ∂ʸu
    #first get ∇q + flux terms
    ∇!(ι.φˣ, ι.φʸ, u, mesh)
    @.  ι.φˣ -= liftx
    @.  ι.φʸ -= lifty

    # Form field differences at faces for x and y partial derivatives
    @. ι.fˣ[:] = ι.φˣ[mesh.vmapM] - (ι.φˣ[mesh.vmapP] + ι.φˣ[mesh.vmapM] )/2
    @. ι.fʸ[:] = ι.φʸ[mesh.vmapM] - (ι.φʸ[mesh.vmapP] + ι.φʸ[mesh.vmapM] )/2

    #enfore boundary conditions for flux (neumann)
    bc_φ!(ι, mesh, dbc)

    #modify with τ, remember fⁿ is field differences at face points
    # which are are overwriting here

    #compute surface term
    @. ι.fⁿ = (mesh.nx * ι.fˣ + mesh.ny * ι.fʸ + τ * ι.fⁿ)

    # compute divergence of flux, volume term
    ∇⨀!(ι.u̇, ι.φˣ, ι.φʸ, mesh)

    # combine the terms
    lift = mesh.lift * (mesh.fscale .* ι.fⁿ)

    tmp =  mesh.J .* ( mesh.M *  ι.u̇)

    @. Δu = tmp

    return nothing
end


# builds the affine operator (one column at a time) (sparse matrix)
# here Δ[u] = L[u] + b (b is where the boundary conditions go as a forcing term)
function constructHelmholtzOperator(ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
    L = spzeros(length(mesh.x), length(mesh.x))
    @. ι.u = 0.0
    Δq = copy(ι.u)
    q =  copy(ι.u)
    b = copy(ι.u)
    @. q = 0
    @. b = 0

    # affine part of operator
    solveHelmholtz!(b, q, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)

    #now construct linear part
    for i in 1:length(mesh.x)
        q[i] = 1.0
        solveHelmholtz!(Δq, q, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
        @. L[:,i] = Δq[:] - b[:]
        q[i] = 0.0
    end
    dropzeros!(L)
    return L, b
end


"""
computeTau(mesh)

# Description

- Computes the tau parameter in NDG

# Arguments

- `mesh` : mesh struct

# Returns

- `τ` : the value of τ at every grid point. (in the code could be either)
"""
function computeTau(mesh)
    matP = mesh.J[mesh.vmapP] ./ mesh.sJ[:]
    matM = mesh.J[mesh.vmapM] ./ mesh.sJ[:]
    hmin = zeros(length(matP))
    for i in 1:length(matP)
        matP[i] < matM[i] ? hmin[i] = 2 * matP[i] : hmin[i] = 2 * matM[i]
    end
    np = (mesh.n + 1) * (mesh.n + 2) / 2
    return reshape(np ./ hmin, mesh.nfp * mesh.nFaces, mesh.K)
end
