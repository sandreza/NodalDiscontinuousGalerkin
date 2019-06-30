include("dg_utils.jl")
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
        D = dmatrix(r, α, β)

        # build surface integral terms
        V = similar(D)
        vandermonde!(V, r, α, β)
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
    uʰ::T
    flux::T

    """
    dg(mesh)

    # Description

        initialize dg struct

    # Arguments

    -   `mesh`: a mesh to compute on

    # Return Values:

    -   `u` : the field to be computed
    -   `uʰ`: numerical solutions for the field
    -   `flux`: the numerical flux for the computation

    """
    function dg(mesh)
        # set up the solution
        u    = copy(mesh.x)
        uʰ   = copy(mesh.x)
        flux = zeros(mesh.nfp * mesh.nfaces, mesh.K)

        return new{typeof(u)}(u, uʰ, flux)
    end
end
