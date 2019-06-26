include("dg1D.jl")

"""
external_params{T,S}

# Description

    struct for external params needed for advection

# Members

    first is velocity
    second is value for α

"""
struct external_params{T,S}
    v::T
    α::S
end

"""
dg_upwind!(uʰ, u, params, t)

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

@btime dg_upwind!(uʰ, u, params, t)
scatter!(x,u, leg = false)

maybe define a function that acts on dg structs?
"""
function dg_upwind!(uʰ, u, params, t)
    # unpack params
    ι = params[1] # internal parameters
    ε = params[2] # external parameters
    periodic = params[3]

    # Form field differences at faces
    diffs = reshape( (u[ι.vmapM] - u[ι.vmapP]), (ι.nfp * ι.nfaces, ι.K ))
    @. ι.du = 1//2 * diffs * (ε.v * ι.nx - (1 - ε.α) * abs(ε.v * ι.nx))

    # Inflow and Outflow boundary conditions
    if !periodic
        uin = -sin(ε.v * t)
        ι.du[ι.mapI]  = @. (u[ι.vmapI] - uin)
        ι.du[ι.mapI] *= @. 1//2 * (ε.v * ι.nx[ι.mapI] - (1-ε.α) * abs(ε.α * abs(ε.v * ι.nx[ι.mapI])))
        ι.du[ι.mapO]  = 0
    end

    # rhs of the semi-discerte PDE, ∂ᵗu = -∂ˣu
    mul!(uʰ, ι.D, u)
    @. uʰ *= -ε.v * ι.rx
    lift = ι.lift * (ι.fscale .* ι.du )
    @. uʰ += lift
    return nothing
end
