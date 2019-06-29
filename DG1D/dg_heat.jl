"""
dg_heat!(uʰ, u, params, t)


# Description

    Evaluate the right hand side for the heat equation

# Example

K = 2^2 # number of elements
n = 2^2-1 # polynomial order
println("The degrees of freedom are ")
println( (n+1)*K)

# domain parameters
xmin = 0.0
xmax = 2π

par_i = dg(K, n, xmin, xmax)
par_e = external_params(1.0, 1.0)
periodic = false
params = (par_i, par_e, periodic)

x = par_i.x
u = par_i.u

@. u = sin(par_i.x) # initial condition
uʰ = par_i.uʰ

@btime dg_heat!(uʰ, u, params, t)
scatter!(x,u, leg = false)

"""
function dg_heat!(uʰ, u, params, t)
    # unpack params
    𝒢 = params[1]
    ι = params[2] # internal parameters
    ε = params[3] # external parameters
    periodic = params[4] #case parameter
    q = params[5]  #temporary arrray for allocation, same size as u
    dq = params[6] #temporary array for allocation, same size as dq
    τ = params[7]   #penalty parameter

    # Form field differences at faces
    diffs = reshape( (u[𝒢.vmapM] - u[𝒢.vmapP]), (𝒢.nfp * 𝒢.nfaces, 𝒢.K ))
    #@. ι.flux = 1//2 * diffs * (ε.v * 𝒢.normals - (1 - ε.α) * abs(ε.v * 𝒢.normals))
    @. ι.flux =  diffs / 2

    # Inflow and Outflow boundary conditions
    if !periodic
        uin  = -u[𝒢.vmapI]
        uout = -u[𝒢.vmapO]
        ι.flux[𝒢.mapI]  =  @. (u[𝒢.vmapI] - uin) / 2
        ι.flux[𝒢.mapO]  =  @. (u[𝒢.vmapO] - uout) / 2
    end

    # rhs of the semi-discerte PDE, ∂ᵗu = ∂ˣq, ∂ˣq  = u
    #first solve for q,
    mul!(q, 𝒢.D, u)
    @. q *= 𝒢.rx
    lift = 𝒢.lift * (𝒢.fscale .* 𝒢.normals .* ι.flux )
    @. q -= lift
    # Form field differences at faces for q
    diffs = reshape( (q[𝒢.vmapM] - q[𝒢.vmapP]), (𝒢.nfp * 𝒢.nfaces, 𝒢.K ))
    #@. dq = 1//2 * diffs * (ε.v * 𝒢.normals - (1 - ε.α) * abs(ε.v * 𝒢.normals))
    @. dq = 0 #reset dq
    @. dq = diffs / 2
    #impose neumann boundary conditions for q
    if !periodic
        qin  = q[𝒢.vmapI]
        qout = q[𝒢.vmapO]
        dq[𝒢.mapI]  =  @. (q[𝒢.vmapI] - qin) / 2
        dq[𝒢.mapO]  =  @. (q[𝒢.vmapO] - qout) / 2
    end
    # solve for uʰ
    mul!(uʰ, 𝒢.D, q)
    @. uʰ *=  𝒢.rx
    lift = 𝒢.lift * (𝒢.fscale .* 𝒢.normals .* dq )
    @. uʰ -= lift
    return nothing
end
