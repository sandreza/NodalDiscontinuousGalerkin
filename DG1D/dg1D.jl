include("dg_utils.jl")
include("mesh.jl")



using Plots
using BenchmarkTools

#need to figure out the fmask situation

struct dg
    u
    du
    rhsu
    K
    n
    np
    nfp
    nfaces
    r
    x
    VX
    EToV
    EToE
    EToF
    fx
    fmask
    vmapM
    vmapP
    vmapB
    mapB
    mapI
    mapO
    vmapI
    vmapO
    D
    lift
    rx
    nx
    fscale
    function dg(KK, nn, xmin, xmax)
        K = KK
        α = 0; β = 0;
        n = nn
        np = n+1
        nfp = 1
        nfaces = 2
        #local element stuff
        r = jacobiGL(α, β, n)
        D = dmatrix(r, α, β)
        V = similar(D)
        jacobi!(V, r, α, β)
        lift = lift1D(V)

        #grid stuff
        nx = normals1D(K)
        VX, EToV = unimesh1D(xmin, xmax, K)
        EToE, EToF = connect1D(EToV)
        x = gridvalues1D(VX, EToV, r)
        rx, J  = geometric_factors(x, D)
        #field stuff
        u = copy(x)
        du = zeros(nfp*nfaces,K)
        rhsu = copy(x)
        #convenience grid stuff
        fx = edgevalues1D(r,x)
        #build fmask
        fmask1 = @. abs(r+1) < eps(1.0);
        fmask2 = @. abs(r-1) < eps(1.0);
        tmp = collect(1:length(r))
        fmask  = [tmp[fmask1]; tmp[fmask2]]
        fmask2 = (fmask1,fmask2)
        fscale = 1 ./ J[fmask,:]
        vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO = buildmaps1D(K, np, nfp, nfaces, fmask2, EToE, EToF, x)
        return new(u, du, rhsu, K, n,np,nfp,nfaces, r, x, VX, EToV, EToE, EToF, fx, fmask,vmapM, vmapP,vmapB,mapB, mapI,mapO,vmapI, vmapO, D, lift, rx, nx, fscale)
    end
end


struct dg_parametric{T,S,U,V,W}
    u::T
    du::T
    rhsu::T
    K::S
    n::S
    np::S
    nfp::S
    nfaces::S
    r::U
    x::T
    VX::U
    EToV::V
    EToE::T
    EToF::T
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
    function dg_parametric(KK, nn, xmin, xmax)
        K = KK
        α = 0; β = 0;
        n = nn
        np = n+1
        nfp = 1
        nfaces = 2
        #local element stuff
        r = jacobiGL(α, β, n)
        D = dmatrix(r, α, β)
        V = similar(D)
        jacobi!(V, r, α, β)
        lift = lift1D(V)
        #grid stuff
        nx = normals1D(K)
        VX, EToV = unimesh1D(xmin, xmax, K)
        EToE, EToF = connect1D(EToV)
        x = gridvalues1D(VX, EToV, r)
        rx, J  = geometric_factors(x, D)
        #field stuff
        u = copy(x)
        du = zeros(nfp*nfaces,K)
        rhsu = copy(x)
        #convenience grid stuff
        fx = edgevalues1D(r,x)
        #build fmask
        fmask1 = @. abs(r+1) < eps(1.0);
        fmask2 = @. abs(r-1) < eps(1.0);
        tmp = collect(1:length(r))
        fmask  = [tmp[fmask1]; tmp[fmask2]]
        fmask2 = (fmask1,fmask2)
        fscale = 1 ./ J[fmask,:]
        vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO = buildmaps1D(K, np, nfp, nfaces, fmask2, EToE, EToF, x)
        return new{typeof(u),typeof(K),typeof(r),typeof(EToV),typeof(vmapP)}(u, du, rhsu, K, n,np,nfp,nfaces, r, x, VX, EToV, EToE, EToF, fx, fmask,vmapM, vmapP,vmapB,mapB, mapI,mapO,vmapI, vmapO, D, lift, rx, nx, fscale)
    end
end


struct dg_collocation{T,S,U,V,W}
    u::T
    du::T
    rhsu::T
    K::S
    n::S
    np::S
    nfp::S
    nfaces::S
    r::U
    x::T
    VX::U
    EToV::V
    EToE::T
    EToF::T
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
    function dg_collocation(KK, nn, xmin, xmax)
        K = KK
        α = 0; β = 0;
        n = nn
        np = n+1
        nfp = 1
        nfaces = 2
        #local element stuff
        r = jacobiGL(α, β, n)
        D = dmatrix(r, α, β)
        V = similar(D)
        jacobi!(V, r, α, β)
        lift = lift1D_v2(V)
        #grid stuff
        nx = normals1D(K)
        VX, EToV = unimesh1D(xmin, xmax, K)
        EToE, EToF = connect1D(EToV)
        x = gridvalues1D(VX, EToV, r)
        rx, J  = geometric_factors(x, D)
        #field stuff
        u = copy(x)
        du = zeros(nfp*nfaces,K)
        rhsu = copy(x)
        #convenience grid stuff
        fx = edgevalues1D(r,x)
        #build fmask
        fmask1 = @. abs(r+1) < eps(1.0);
        fmask2 = @. abs(r-1) < eps(1.0);
        tmp = collect(1:length(r))
        fmask  = [tmp[fmask1]; tmp[fmask2]]
        fmask2 = (fmask1,fmask2)
        fscale = 1 ./ J[fmask,:]
        vmapM, vmapP, vmapB, mapB, mapI, mapO, vmapI, vmapO = buildmaps1D(K, np, nfp, nfaces, fmask2, EToE, EToF, x)
        return new{typeof(u),typeof(K),typeof(r),typeof(EToV),typeof(vmapP)}(u, du, rhsu, K, n,np,nfp,nfaces, r, x, VX, EToV, EToE, EToF, fx, fmask,vmapM, vmapP,vmapB,mapB, mapI,mapO,vmapI, vmapO, D, lift, rx, nx, fscale)
    end
end
