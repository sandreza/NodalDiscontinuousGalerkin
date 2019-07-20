# module for all the seperate pieces of navier-stokes

#



function step_euler()
end

"""
step_ab!(uᵀ, vᵀ, u¹, u¹, u⁰, v⁰ nu¹, nv¹, nu², nv², ab, dt, γ)

# Description

- computes an adam - bashforth step

# Arguments

- `!uᵀ` : nonlinear explicit term
- `!vᵀ` : nonlinear explicit term
- `u¹` : u-velocity at the current timestep
- `v¹` : v-velocity at the current timestep
- `u⁰` : u-velocity at the previous timestep
- `v⁰` : v-velocity at the previous timestep
- `nu²` : nonlinear u-velocity term at the current timestep
- `nv²` : nonlinear v-velocity term at the current timestep
- `nu¹` : nonlinear u-velocity term at the previous timestep
- `nv¹` : nonlinear v-velocity term at the previous timestep
- `ab` : adam-bashforth weights
- `dt` : timestep size
- `γ` : parameter

"""
function step_ab!(uᵀ, vᵀ, u¹, v¹, u⁰, v⁰, nu¹, nv¹, nu², nv², ab, dt, γ)
        @. uᵀ = ( (ab[1,1] * u¹ + ab[2,1] * u⁰) - dt * (ab[1,2] * nu² + ab[2,2] * nu¹) ) / γ
        @. vᵀ = ( (ab[1,1] * v¹ + ab[2,1] * v⁰) - dt * (ab[1,2] * nv² + ab[2,2] * nv¹) ) / γ
        return nothing
end

"""
ns_flux(φux, φuy, φvx, φvy, u, v)

# Description

- Compute advection term in Navier-Stokes equations. used in advection

# Arguments

- `φux`: flux for the u velocity in the x-direction
- `φuy`: flux for the u velocity in the y-direction
- `φvx`: flux for the v velocity in the x-direction
- `φvy`: flux for the v velocity in the y-direction
- `u`  : velocity field in the x-direction
- `v`  : velocity field in the y-direction

# Return : nothing

"""
function ns_flux!(φux, φuy, φvx, φvy, u, v)
        @. φux = u^2
        @. φuy = u * v
        @. φvx = φuy
        @. φvy = v^2
end

"""
nonlinear(nu¹, nv¹, nu², nv², φux, φuy, φvx, φvy, 𝒢)

# Description

- computes the nonlinear term and saves the old one.

# Arguments

- `nu¹` : nonlinear term in the u-velocity field at the previous timestep
- `nv¹` : nonlinear term in the v-velocity field at the previous timestep
- `nu²` : nonlinear term in the u-velocity field at the current timestep
- `nv²` : nonlinear term in the v-velocity field at the current timestep
- `φux`: flux for the u velocity in the x-direction
- `φuy`: flux for the u velocity in the y-direction
- `φvx`: flux for the v velocity in the x-direction
- `φvy`: flux for the v velocity in the y-direction
"""
function nonlinear!(nu¹, nv¹, nu², nv², φux, φuy, φvx, φvy, 𝒢)
        # save old values
        @. nu¹ = nu²
        @. nv¹ = nu²
        #compute new nonlinear term, may be a good idea to split advection here
        ∇⨀!(nu², φux, φuy, 𝒢)
        ∇⨀!(nv², φvx, φvy, 𝒢)
end

"""
face_velocity!(u⁻, v⁻, u⁺, v⁺, u, v, mesh)

# Description

- compute the velocity field on the face

# Arguments

- `!u⁻`: u-velocity on the face interior to the node
- `!v⁻`: v-velocity on the face interior to the node
- `!u⁺`: u-velocity on the face exterior to the node
- `!v⁺`: v-velocity on the face exterior to the node
- `u`   : u-velocity at every grid point
- `v`   : v-velocity at every grid point
- `mesh`: mesh struct

"""
function face_velocity!(u⁻, v⁻, u⁺, v⁺, u, v, mesh)
        @. u⁻ = u[mesh.vmapM]
        @. v⁻ = v[mesh.vmapM]
        @. u⁺ = u[mesh.vmapP]
        @. v⁺ = v[mesh.vmapP]
        return nothing
end



"""
normal_face_velocity!(un⁻,un⁺, u⁻, v⁻, u⁺, v⁺, mesh)

# Description

- compute the maximum normal velocity field on the face

# Arguments

- `!n⁻`: normal velocity on the face interior to the node
- `!n⁺`: normal velocity on the face exterior to the node
- `u⁻`: u-velocity on the face interior to the node
- `v⁻`: v-velocity on the face interior to the node
- `u⁺`: u-velocity on the face exterior to the node
- `v⁺`: v-velocity on the face exterior to the node
- `mesh`: mesh struct

"""
function normal_face_velocity!(n⁻, n⁺, u⁻, v⁻, u⁺, v⁺, mesh)
        @. n⁻ = mesh.nx * u⁻ + mesh.ny * v⁻
        @. n⁺ = mesh.nx * u⁺ + mesh.ny * v⁺
        return nothing
end

"""
maximum_face_velocity!(maxv, n⁻, n⁺, mesh)

# Description

- compute the maximum velocity field on the face. allocates memory.

# Arguments

- `!maxv`: maximum face velocity
- `!n⁻`: v-velocity on the face interior to the node
- `!n⁺`: u-velocity on the face exterior to the node
- `!v⁺`: v-velocity on the face exterior to the node
- `u`   : u-velocity at every grid point
- `v`   : v-velocity at every grid point
- `mesh`: mesh struct

"""
function maximum_face_velocity!(maxv, n⁻, n⁺,  mesh)
        maxtmp = [ maximum([n⁻[i] n⁺[i]]) for i in 1:length(a) ]
        # reorder so that we can just take the max along a given dimension
        # the output will be the maximum along each face in linear ordering
        maxtmp = reshape(maxtmp, mesh.nfp, mesh.nfaces *  mesh.k)
        # duplicate values on the face
        maxtmp = ones(mesh.nfp, 1) * maximum(maxtmp, dims =  1)
        maxtmp = reshape(maxtmp, mesh.nfp * mesh.nfaces, mesh.K)
        @. maxv = maxtmp
        return nothing
end


"""
face_flux(φux, φuy, φvx, φvy, u, v)

# Description

- calculate the flux on the faces

# Arguments

- `φux⁻`: interior flux for the u velocity in the x-direction
- `φuy⁻`: interior flux for the u velocity in the y-direction
- `φvx⁻`: interior flux for the v velocity in the x-direction
- `φvy⁻`: interior flux for the v velocity in the y-direction
- `φux⁺`: exterior flux for the u velocity in the x-direction
- `φuy⁺`: exterior flux for the u velocity in the y-direction
- `φvx⁺`: exterior flux for the v velocity in the x-direction
- `φvy⁺`: exterior flux for the v velocity in the y-direction
- `u⁻`: interior  u velocity on a face
- `v⁻`: interior  v velocity on a face
- `u⁺`: exterior  u velocity on a face
- `v⁺`: exterior  v velocity on a face
- `mesh`  : velocity field in the x-direction

# Return : nothing

"""
function face_flux!(φux⁻, φuy⁻, φvx⁻, φvy⁻, φux⁺, φuy⁺, φvx⁺, φvy⁺, u⁻, v⁻, u⁺, v⁺, mesh)
        # interior face
        @. φux⁻ = u⁻[mesh.vmapM] * u⁻[mesh.vmapM]
        @. φuy⁻ = u⁻[mesh.vmapM] * v⁻[mesh.vmapM]
        @. φvx⁻ = v⁻[mesh.vmapM] * u⁻[mesh.vmapM]
        @. φvy⁻ = v⁻[mesh.vmapM] * v⁻[mesh.vmapM]

        # exterior face
        @. φux⁺ = u⁺[mesh.vmapM] * u⁺[mesh.vmapM]
        @. φuy⁺ = u⁺[mesh.vmapM] * v⁺[mesh.vmapM]
        @. φvx⁺ = v⁺[mesh.vmapM] * u⁺[mesh.vmapM]
        @. φvy⁺ = v⁺[mesh.vmapM] * v⁺[mesh.vmapM]
end

"""
ns_rusonov_flux!(sφu, sφv, ux⁻, φuy⁻, φvx⁻, φvy⁻, φux⁺, φuy⁺, φvx⁺, φvy⁺, u⁻, v⁻, u⁺, v⁺, mesh)

# Description

- Calulate the total flux on the faces

# Arguments

- `sφu` : surface flux for the u-velocity
- `sφv` : surface flux for the v-velocity
- `φux⁻`: interior flux for the u velocity in the x-direction
- `φuy⁻`: interior flux for the u velocity in the y-direction
- `φvx⁻`: interior flux for the v velocity in the x-direction
- `φvy⁻`: interior flux for the v velocity in the y-direction
- `φux⁺`: exterior flux for the u velocity in the x-direction
- `φuy⁺`: exterior flux for the u velocity in the y-direction
- `φvx⁺`: exterior flux for the v velocity in the x-direction
- `φvy⁺`: exterior flux for the v velocity in the y-direction
- `u⁻`: interior  u velocity on a face
- `v⁻`: interior  v velocity on a face
- `u⁺`: exterior  u velocity on a face
- `v⁺`: exterior  v velocity on a face
- `maxv`: maximum velocity on a face
- `mesh`  : velocity field in the x-direction
"""
function ns_rusonov_flux!(sφu, sφv, ux⁻, φuy⁻, φvx⁻, φvy⁻, φux⁺, φuy⁺, φvx⁺, φvy⁺, u⁻, v⁻, u⁺, v⁺, maxv, mesh)
        # yes the signs are correct on the last entry
        @. sφu = - mesh.nx * ( φux⁻ - φux⁺) - mesh.ny * ( φuy⁻ - φuy⁺) - maxv * (u⁺ - u⁻)
        @. sφv = - mesh.nx * ( φvx⁻ - φvx⁺) - mesh.ny * ( φvy⁻ - φvy⁺) - maxv * (v⁺ - v⁻)

        @. sφu *= 0.5
        @. sφv *= 0.5
        return nothing
end


"""
explicit_nonlinear_rhs!(nu, nv, sφu, sφv, mesh)

# Description

- computes the nonlinear term for navier-stokes. the other terms are implicit or pressure

# Arguments

- `nu`   : nonlinear volume term for the u-velocity
- `nv`   : nonlinear volume term for the v-velocity
- `sφu`  : surface flux for the u-velocity
- `sφv`  : surface flux for the v-velocity
- `mesh` : mesh struct

"""
function explicit_nonlinear_rhs!(nu, nv, sφu, sφv, mesh)
        liftu = mesh.lift * (mesh.fscale .* sφu)
        liftv = mesh.lift * (mesh.fscale .* sφv)
        @. nu += liftu
        @. nv += liftv
        return nothing
end


"""
enforce_bc!(uf⁺, vf⁺, bc_u, bc_v, mapT)

# Description

- enforce boundary conditions by utilizing exterior face nodes

# Arguments

- `!uf⁺`: u-velocity on the face exterior to the node
- `!vf⁺`: v-velocity on the face exterior to the node
- `bc_u`: boundary condition for the u-velocity
- `bc_v`: boundary condition for the v-velocity
- `mapT`: a tuple of arrays that correspond to boundary conditions
"""
function enforce_bc!(uf⁺, vf⁺, bc_u, bc_v, mapT)
        for i in 1:length(mapT)
                uf⁺[mapT[i]] = bc_u[mapT[i]]
                vf⁺[mapT[i]] = bc_v[mapT[i]]
        end
        return nothing
end


function pressure_solve()
end

function viscous_step()
end

"""
pearson_vortex!(u, v, 𝒢, t)

# Description

- An exact solution to the Navier-Stokes equations

# Arguments

- `u` : velocity field component in the x-direction
- `v` : veloctiy field component in the y-direction
- `p` : pressure field
- `𝒢` : grid struct
- `t` : time
"""
function pearson_vortex!(u, v, p, 𝒢, ν, t)
        @.  u = -sin(2 * pi * 𝒢.y ) * exp( - ν * 4 * pi^2 * t)
        @.  v =  sin(2 * pi * 𝒢.x ) * exp( - ν * 4 * pi^2 * t)
        @.  p = -cos(2 * pi * 𝒢.x ) * cos(2 * pi * 𝒢.y) * exp( - ν * 8 * pi^2 * t)
end
