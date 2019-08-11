"""
dg_stokes_bc!(𝒮ϕ, ϕ, ns, params, mesh, bc_ϕ!, bc, bc_φ!, dbc)


# Description

- Evaluate the right hand side for (helmholtz) stoke's equation

# Arguments

- `𝒮ϕ` : the stokes operator acting on ϕ = (u,v,p)
- `ϕ` :  the thing we want to take laplacian of
- `ns` : navier-stokes struct for temporary variables
- `params`: any penalty parameters that we would like to include
- `mesh` : the mesh struct with all the grid information
- `bc_ϕ!`: function that computes dirichlet boundary conditions
- `bc` : boundary condition tuple with indices
- `bc_φ!`: function that computes neumann boundary conditions
- `dbc` : boundary condition tuple with indices

# return

- nothing


"""
function dg_stokes_bc!(𝒮ϕ, ϕ, ns, params, mesh, bc!, bc_ϕ, bc_∇!, dbc_ϕ)
    # unpack parameters
    τ = params[1]
    γ = params[2]

    # unpack boundary conditions
    # u-velocity
    bc_u = bc_ϕ[1]
    dbc_u = dbc_ϕ[1]
    # v-velocity
    bc_v = bc_ϕ[2]
    dbc_v = dbc_ϕ[2]
    # pressure
    bc_p = bc_ϕ[3]
    dbc_p = dbc_ϕ[3]

    # unpack ϕ and 𝒮ϕ
    u = ϕ[:,:, 1]
    v = ϕ[:,:, 2]
    p = ϕ[:,:, 3]

    @. ns.u.ϕ = u
    @. ns.v.ϕ = v
    @. ns.p.ϕ = p

    ℋu = similar(mesh.x)
    ℋv = similar(mesh.x)
    #compute the easy block operators
    dg_helmholtz_bc!(ℋu, u, ns.u, params, mesh, bc!, bc_u, bc_∇!, dbc_u)
    dg_helmholtz_bc!(ℋv, v, ns.v, params, mesh, bc!, bc_v, bc_∇!, dbc_v)

    # compute off diagonal terms, pressure terms, derivative + lift
    ∂ˣp = ∂ˣ_∮(ns.p, mesh, bc!, bc_p)
    ∂ʸp = ∂ʸ_∮(ns.p, mesh, bc!, bc_p)

    # compute divergence condition
    ∂ˣu = ∂ˣ_∮(ns.u, mesh, bc!, bc_u)
    ∂ʸv = ∂ʸ_∮(ns.v, mesh, bc!, bc_v)

    #penalty = mesh.lift * reshape((mesh.fscale[:] .* (τ[:] .* (p[mesh.vmapM] - (p[mesh.vmapM] + p[mesh.vmapP])))/2), size(mesh.fscale))
    #penalty *= 0.0

    𝒮ᵘ = ℋu - mesh.J .*  (mesh.M * ∂ˣp)
    𝒮ᵛ = ℋv - mesh.J .*  (mesh.M * ∂ʸp)
    𝒮ᵖ = mesh.J .*  (mesh.M * ( ∂ˣu + ∂ʸv ))

    @. 𝒮ϕ[:,:,1] = 𝒮ᵘ
    @. 𝒮ϕ[:,:,2] = 𝒮ᵛ
    @. 𝒮ϕ[:,:,3] = 𝒮ᵖ

    return nothing
end

# lift + volume helper functions, uses central difference
function ∂ˣ_∮(ι, mesh, bc_ϕ!, bc)
    # form field differnces at faces
    @. ι.fⁿ[:] =  (ι.ϕ[mesh.vmapM] - ι.ϕ[mesh.vmapP]) / 2 #central flux
    # enforce bc
    bc_ϕ!(ι, mesh, bc)
    # compute normal component in the x-direction
    @. ι.fˣ = mesh.nx * ι.fⁿ
    # compute lift term
    liftx = mesh.lift * (mesh.fscale .* ι.fˣ )
    # compute partial with respect to x
    ∇!(ι.∂ˣ, ι.∂ʸ, ι.ϕ, mesh)
    return ∂ˣϕ =  ι.∂ˣ + liftx
end

function ∂ʸ_∮(ι, mesh, bc_ϕ!, bc)
    # form field differnces at faces
    @. ι.fⁿ[:] =  (ι.ϕ[mesh.vmapM] - ι.ϕ[mesh.vmapP]) / 2 #central flux
    # enforce bc
    bc_ϕ!(ι, mesh, bc)
    # compute normal component in the x-direction
    @. ι.fʸ = mesh.ny * ι.fⁿ
    # compute lift term
    lifty = mesh.lift * (mesh.fscale .* ι.fʸ )
    # compute partial with respect to x
    ∇!(ι.∂ˣ, ι.∂ʸ, ι.ϕ, mesh)
    return ∂ʸϕ =  ι.∂ʸ + lifty
end


function stokes_setup_bc(ϕ, ns, params, mesh, bc!, bc_ϕ, bc_∇!, dbc_ϕ)
    L = spzeros(length(ϕ), length(ϕ))
    @. ϕ = 0.0
    𝒮ϕ = copy(ϕ)
    q =  copy(ϕ)
    b = copy(ϕ)
    @. q = 0
    @. b = 0
    # affine part of operator
    dg_stokes_bc!(b, q, ns, params, mesh, bc!, bc_ϕ, bc_∇!, dbc_ϕ)
    @. q = 0 #just in case
    #now construct linear part
    for i in 1:length(ϕ)
        q[i] = 1.0
        dg_stokes_bc!(𝒮ϕ, q, ns, params, mesh, bc!, bc_ϕ, bc_∇!, dbc_ϕ)
        @. L[:,i] = 𝒮ϕ[:] - b[:]
        q[i] = 0.0
        dropϵzeros!(L)
    end
    return L, b
end


function modify_stokes_operator(L, b)
    m,n = size(L) #divisible by three
    mr = Int(m/3)
    nr = Int(n/3)

    #new operator L
    nL = spzeros(m+1,n+1)
    nb = zeros(n+1)
    @. nL[1:m, 1:n] = L
    # pad with enforcement of Lagrange multipliers to make problem invertible
    @. nL[(2*mr+1):m,n+1] = 1.0
    @. nL[m+1,(2*nr+1):n] = 1.0
    dropϵzeros!(nL)
    # set entries for new b
    @. nb[1:n] = b[:]
    return nL, nb
end
