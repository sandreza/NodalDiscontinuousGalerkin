include("dg1D.jl")

"""
dg_upwind!(du, u, params, t)

#Example

K = 2^2 #number of elements
n = 2^2-1 #polynomial order
println("The degrees of freedom are ")
println( (n+1)*K)
#domain parameters
xmin = 0.0
xmax = 2π
par_i = dg(K, n, xmin, xmax)
struct extern_params
    v
    α
end
par_e = extern_params(1.0, 1.0) #first is velocity, second is value for α
params = (par_i, par_e)
x = par_i.x
u = par_i.u
@. u = sin(par_i.x) #initial condition
rhsu = par_i.rhsu
j
@btime dg_upwind!(rhsu, u, params, t)
scatter!(x,u, leg = false)

maybe define a function that acts on dg structs?
"""
function dg_upwind!(rhsu, u, params, t)
    ι = params[1] #internal parameters, this is a struct
    ε = params[2] #external parameters, advection velocity and alpha
    #Internal node stuff
    #@.  du = (u[ι.vmapM] - u[ι.vmapP]) * (ε.v * ι.nx - (1 - ε.α) * abs(ε.v * ι.nx) ) /2

    #interface terms, calculate flux
    tmp1 = reshape(u[ι.vmapM] - u[ι.vmapP], ( ι.nfp * ι.nfaces , ι.K ))
    @. ι.du  = tmp1 * (ε.v * ι.nx - (1 - ε.α) * abs(ε.v * ι.nx) ) /2

    #Inflow and Outflow boundary conditions
    uin = -sin(ε.v * t)
    ι.du[ι.mapI] = (u[ι.vmapI] .- uin)
    ι.du[ι.mapI] *= (ε.v .* ι.nx[ι.mapI] .- (1-ε.α) .* abs.(ε.α .* abs.( ε.v .* ι.nx[ι.mapI] ))) ./ 2
    ι.du[ι.mapO] = 0

    #rhs of the semi-discerte PDE, ∂_t u = - ∂_x u
    mul!(rhsu, ι.D, u)
    @. rhsu *= -ε.v * ι.rx
    tmp2 = ι.lift * (ι.fscale .* ι.du )
    @. rhsu += tmp2
    return nothing
end

#for periodic functions
function dg_upwind_p!(rhsu, u, params, t)
    ι = params[1] #internal parameters, this is a struct
    ε = params[2] #external parameters, advection velocity and alpha
    #Internal node stuff
    #@.  du = (u[ι.vmapM] - u[ι.vmapP]) * (ε.v * ι.nx - (1 - ε.α) * abs(ε.v * ι.nx) ) /2

    #interface terms, calculate flux
    tmp1 = reshape(u[ι.vmapM] - u[ι.vmapP], ( ι.nfp * ι.nfaces , ι.K ))
    @. ι.du  = tmp1 * (ε.v * ι.nx - (1 - ε.α) * abs(ε.v * ι.nx) ) /2

    #rhs of the semi-discerte PDE, ∂_t u = - ∂_x u
    mul!(rhsu, ι.D, u)
    @. rhsu *= -ε.v * ι.rx
    tmp2 = ι.lift * (ι.fscale .* ι.du )
    @. rhsu += tmp2
    return nothing
end
