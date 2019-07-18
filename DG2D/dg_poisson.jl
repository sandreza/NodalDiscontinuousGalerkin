# builds the 2D poisson matrix
"""
dg_poisson!(Δu, u, params, t)


# Description

- Evaluate the right hand side for poisson's equation

# Arguments

- `Δu` : the laplacian of u
- `u` :  the thing we want to take laplacian of
- `ι` : struct for temporary variables
- `params`: any penalty parameters that we would like to include
- `mesh` : the mesh struct with all the grid information
- `bc_function!`: function that computes boundary conditions
- `bc` : boundary condition tuple with indices


"""
function dg_poisson!(Δu, u, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
    # unpack parameters
    τ = params[1]

    # Form field differences at faces
    @. ι.fⁿ[:] =  (u[mesh.vmapM] - u[mesh.vmapP])

    # Choose boundary condition type, dirichlet
    bc_u!(ι.fⁿ, u, bc)

    # get the numerical flux jump in the x and y directions
    @. ι.fˣ = mesh.nx * ι.fⁿ / 2
    @. ι.fʸ = mesh.ny * ι.fⁿ / 2

    liftx = mesh.lift * (mesh.fscale .* ι.fˣ )
    lifty = mesh.lift * (mesh.fscale .* ι.fʸ )
    # lhs of the semi-discerte PDE, ∇⋅(q) = f , q  = ∇u, qˣ = ∂ˣu, qʸ = ∂ʸu
    #first get ∇q + flux terms
    ∇!(ι.φˣ, ι.φʸ, u, mesh)
    @.  ι.φˣ -= liftx
    @.  ι.φʸ -= lifty

    # Form field differences at faces for x and y partial derivatives
    @. ι.fˣ[:] = ι.φˣ[mesh.vmapM] - ι.φˣ[mesh.vmapP]
    @. ι.fʸ[:] = ι.φʸ[mesh.vmapM] - ι.φʸ[mesh.vmapP]

    #enfore boundary conditions for flux (neumann)
    bc_φ!(ι.fˣ, ι.fʸ,ι.φˣ, ι.φʸ, dbc)

    #modify with τ, remember fⁿ is field differences at face points
    # which are are overwriting here

    #compute surface term
    @. ι.fⁿ = (mesh.nx * ι.fˣ + mesh.ny * ι.fʸ + τ * ι.fⁿ)/ 2

    # compute divergence of flux, volume term
    ∇⨀!(ι.u̇, ι.φˣ, ι.φʸ, mesh)

    # combine the terms
    tmp =  mesh.J .* ( mesh.M * (ι.u̇ - mesh.lift * (mesh.fscale .* ι.fⁿ) ))
    @. Δu = tmp
    return nothing
end



#builds the matrix (one column at a time) (sparse matrix)
function poisson_setup(ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
    L = spzeros(length(mesh.x), length(mesh.x))
    @. ι.u = 0.0
    Δq = copy(ι.u)
    q =  copy(ι.u)
    @. q = 0
    for i in 1:length(mesh.x)
        q[i] = 1.0
        dg_poisson!(Δq, q, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
        @. L[:,i] = Δq[:]
        q[i] = 0.0
    end
    dropzeros!(L)
    return L
end
