include("../utils.jl")
include("mesh.jl")

using Plots
using BenchmarkTools


struct mesh{T,S,U,W}
    # inputs
    K::S
    n::S

    # face stuff
    nfp::S
    nfaces::S

    # GL points
    r::U
    x::T

    # vertex maps
    vmapM::W
    vmapP::W
    vmapB::W
    mapB::W

    # inflow/outflow maps
    mapI::S
    mapO::S
    vmapI::S
    vmapO::S

    # structures for computation
    D::T
    M::T
    Mi::T
    lift::T
    rx::T
    normals::T
    fscale::T

    """
    mesh(KK, nn, xmin, xmax)

    # Description

        initialize mesh struct

    # Arguments

        KK: number of elements
        nn: polynomial order
        xmin: lower bound
        xmax: upper bound


    # Return Values: x
        return grid values

    """
    function mesh(KK, nn, xmin, xmax)
        # initialize parameters
        K = KK
        α = 0; β = 0;
        n = nn

        # number of vertices
        np = n+1
        nfp = 1
        nfaces = 2

        # compute Gauss Lobatto grid
        r = jacobiGL(α, β, n)

        # build grid
        VX, EtoV = unimesh1D(xmin, xmax, K)

        # build coordinates of all the nodes
        x = gridvalues1D(VX, EtoV, r)

        # build connectivity matrix
        EtoE, EtoF = connect1D(EtoV)

        # build face masks
        fmask1,fmask2 = fmask1D(r)

        # build connectibity maps
        vmapM,vmapP,vmapB,mapB, mapI,mapO,vmapI,vmapO = buildmaps1D(K, np,nfp,nfaces, fmask1, EtoE,EtoF, x)

        # build differentiation matrix
        D = dmatrix(r, α, β, n)

        # build surface integral terms
        V = vandermonde(r, α, β, n)
        lift = ∮dΩ(V)

        # build mass matrix and inverse of mass matrix
        Mi = V * V'
        M = inv(Mi)

        # calculate geometric factors
        rx,J = geometric_factors(x, D)

        # build surface normals
        normals = normals1D(K)

        # build inverse metric at the surface
        fscale = 1 ./ J[fmask2,:]

        return new{typeof(x),typeof(K),typeof(r),typeof(vmapP)}(K,n, nfp,nfaces, r,x, vmapM,vmapP,vmapB,mapB, mapI,mapO,vmapI,vmapO, D,M,Mi,lift,rx,normals,fscale)
    end
end

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
