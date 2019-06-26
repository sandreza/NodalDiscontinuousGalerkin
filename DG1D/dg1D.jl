include("dg_utils.jl")
include("mesh.jl")

using Plots
using BenchmarkTools

struct dg{T,S,U,V,W}
    u::T
    du::T
    uʰ::T

    K::S
    n::S

    np::S
    nfp::S
    nfaces::S

    r::U
    x::T

    VX::U
    EtoV::V

    EtoE::T
    EtoF::T

    fx::T
    fmask::W

    vmapM::W
    vmapP::W
    vmapB::W
    mapB::W

    mapI::S
    mapO::S
    vmapI::S
    vmapO::S

    D::T
    lift::T
    rx::T
    nx::T
    fscale::T

    """
    dg(KK, nn, xmin, xmax)

    # Description

        initialize DG struct

    # Arguments

        KK: number of elements
        nn: polynomial order
        xmin: lower bound
        xmax: upper bound


    # Return Values: x
        all the members of the DG struct

    """
    function dg(KK, nn, xmin, xmax)
        # initialize parameters
        K = KK
        α = 0; β = 0;
        n = nn

        # number of vertices
        np = n+1
        nfp = 1
        nfaces = 2

        # comput Gauss Lobatto grid
        r = jacobiGL(α, β, n)

        # build reference element matrices
        D = dmatrix(r, α, β)
        V = similar(D)
        vandermonde!(V, r, α, β)

        # build surface integral terms
        lift = ∮dΩ(V)

        # build grid
        VX, EtoV = unimesh1D(xmin, xmax, K)

        # build coordinates of all the nodes
        x = gridvalues1D(VX, EtoV, r)

        # calculate geometric factors
        rx, J  = geometric_factors(x, D)

        # compute masks for edge nodes
        fmask1,fmask2 = fmask1D(r)
        fx = edgevalues1D(fmask1,x)

        # build surface normals and inverse metric at the surface
        nx = normals1D(K)
        fscale = 1 ./ J[fmask2,:]

        # build connectivity matrix
        EtoE, EtoF = connect1D(EtoV)

        # set up the solution
        u = copy(x)
        du = zeros(nfp*nfaces,K)
        uʰ = copy(x)

        # build connectibity maps
        vmapM,vmapP,vmapB,mapB, mapI,mapO,vmapI,vmapO = buildmaps1D(K, np,nfp,nfaces, fmask1, EtoE,EtoF, x)

        return new{typeof(u),typeof(K),typeof(r),typeof(EtoV),typeof(vmapP)}(u,du,uʰ, K,n, np,nfp,nfaces, r,x, VX,EtoV, EtoE,EtoF, fx,fmask2, vmapM,vmapP,vmapB,mapB, mapI,mapO,vmapI,vmapO, D,lift,rx,nx,fscale)
    end
end
