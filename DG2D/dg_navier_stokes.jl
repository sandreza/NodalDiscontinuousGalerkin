# module for all the seperate pieces of navier-stokes

#



function step_euler()
end

function step_ab()
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
nonlinear()

# Description

- computes the nonlinear term and saves the old one.

# Arguments

- `nu¹` : nonlinear term in the u-velocity field at the previous timestep
- `nv¹` : nonlinear term in the v-velocity field at the previous timestep
- `nu²` : nonlinear term in the u-velocity field at the current timestep
- `nv²` : nonlinear term in the v-velocity field at the current timestep
"""
function nonlinear!(nu¹, nv¹, nu², nv², φux, φuy, φvx, φvy, 𝒢)
        # save old values
        @. nu¹ = nu²
        @. nv¹ = nu²
        #compute new nonlinear term, may be a good idea to split advection here
        ∇⨀!(nu², φux, φuy, 𝒢)
        ∇⨀!(nv², φvx, φvy, 𝒢)
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

function pearson_vortex!(u, v, p, 𝒢, t)
        @.  u = -sin(2 * pi * 𝒢.y ) * exp( - nu * 4 * pi^2 * t)
        @.  v =  sin(2 * pi * 𝒢.x ) * exp( - nu * 4 * pi^2 * t)
        @.  p = -cos(2 * pi * 𝒢.x ) * cos(2 * pi * 𝒢.y) * exp( - nu * 8 * pi^2 * t)
end
