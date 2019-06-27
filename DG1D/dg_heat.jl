"""
dg_heat!(uʰ, u, params, t)

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
    ι = params[1] # internal parameters
    ε = params[2] # external parameters
    periodic = params[3] #case parameter
    q = params[4]  #temporary arrray for allocation, same size as u
    dq = params[5] #temporary array for allocation, same size as dq

    # Form field differences at faces
    diffs = reshape( (u[ι.vmapM] - u[ι.vmapP]), (ι.nfp * ι.nfaces, ι.K ))
    #@. ι.du = 1//2 * diffs * (ε.v * ι.nx - (1 - ε.α) * abs(ε.v * ι.nx))
    @. ι.du =  diffs / 2

    # Inflow and Outflow boundary conditions
    if !periodic
        uin  = -u[ι.vmapI]
        uout = -u[ι.vmapO]
        ι.du[ι.mapI]  =  @. (u[ι.vmapI] - uin) / 2
        ι.du[ι.mapO]  =  @. (u[ι.vmapO] - uout) / 2
    end

    # rhs of the semi-discerte PDE, ∂ᵗu = ∂ˣq, ∂ˣq  = u
    #first solve for q
    mul!(q, ι.D, u)
    @. q *= ι.rx
    lift = ι.lift * (ι.fscale .* ι.nx .* ι.du )
    @. q -= lift
    # Form field differences at faces for q
    diffs = reshape( (q[ι.vmapM] - q[ι.vmapP]), (ι.nfp * ι.nfaces, ι.K ))
    #@. dq = 1//2 * diffs * (ε.v * ι.nx - (1 - ε.α) * abs(ε.v * ι.nx))
    @. dq = 0 #reset dq
    @. dq = diffs / 2
    #impose neumann boundary conditions for q
    if !periodic
        qin  = q[ι.vmapI]
        qout = q[ι.vmapO]
        dq[ι.mapI]  =  @. (q[ι.vmapI] - qin) / 2
        dq[ι.mapO]  =  @. (q[ι.vmapO] - qout) / 2
    end
    # solve for uʰ
    mul!(uʰ, ι.D, q)
    @. uʰ *=  ι.rx
    lift = ι.lift * (ι.fscale .* ι.nx .* dq )
    @. uʰ -= lift
    return nothing
end
