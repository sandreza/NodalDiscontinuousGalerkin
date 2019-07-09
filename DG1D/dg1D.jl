include("mesh.jl")

using Plots
using BenchmarkTools

struct dg{T}
    u::T
    u̇::T
    flux::T

    """
    dg(mesh)

    # Description

        initialize dg struct

    # Arguments

    -   `mesh`: a mesh to compute on

    # Return Values:

    -   `u` : the field to be computed
    -   `u̇`: numerical solutions for the field
    -   `flux`: the numerical flux for the computation

    """
    function dg(mesh)
        # set up the solution
        u    = copy(mesh.x)
        u̇   = copy(mesh.x)
        flux = zeros(mesh.nfp * mesh.nfaces, mesh.K)

        return new{typeof(u)}(u, u̇, flux)
    end
end

# low storage Runge-Kutta coefficients
rk4a = [ 0.0, -567301805773.0/1357537059087.0, -2404267990393.0/2016746695238.0, -3550918686646.0/2091501179385.0, -1275806237668.0/842570457699.0]
rk4b = [ 1432997174477.0/9575080441755.0, 5161836677717.0/13612068292357.0, 1720146321549.0/2090206949498.0, 3134564353537.0/4481467310338.0, 2277821191437.0/14882151754819.0]
rk4c = [ 0.0, 1432997174477.0/9575080441755.0, 2526269341429.0/6820363962896.0, 2006345519317.0/3224310063776.0, 2802321613138.0/2924317926251.0]

"""
rk_solver!(u̇, u, params, t)

# Description

    time stepping with 4th order runge-kutta

# Arguments

-   `u̇ = (Eʰ, Hʰ)`: container for numerical solutions to fields
-   `u  = (E , H )`: container for starting field values
-   `params = (𝒢, E, H, ext)`: mesh, E sol, H sol, and material parameters
-   `t`: time to evaluate at

"""
function rk_solver!(rhs!, u̇, u, params, tspan, dt)
    # Runge-Kutta residual storage
    nsol = length(u)
    res = Any[]
    for iRes in 1:nsol
        push!(res, zeros(size(u[iRes])))
    end

    # store solutions at all times
    Nsteps = ceil(Int, tspan[end] / dt)
    sol = Any[]

    # time step loop
    for tstep in 1:Nsteps
        for iRK in 1:5
            # get numerical solution
            rhs!(u̇, u, params, dt)

            # update solutions
            for iRes in 1:nsol
                res[iRes] = rk4a[iRK] * res[iRes] + dt * u̇[iRes]
                u[iRes] = u[iRes] + rk4b[iRK] * res[iRes]
                # seems to differ from matlab code during this step ???
            end
        end

        uᵗ = similar(u)
        @. uᵗ = u
        push!(sol, uᵗ)

        if (tstep % 10000) == 0
            println( string(tstep, " / ", Nsteps))
        end
    end

    return sol
end
