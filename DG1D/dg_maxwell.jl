include("dg1D.jl")

"""
material_params{T}

# Description

    struct for material params needed for maxwell's equations

# Members

    ϵ is the electric permittivity
    μ is the magnetic permeability

"""
struct material_params{T}
    ϵ::T
    μ::T
end

"""
dg_maxwell!(uʰ, u, params, t)

# Description

    numerical solution to 1D maxwell's equation

# Arguments

-   `uʰ = (Eʰ, Hʰ)`: container for numerical solutions to fields
-   `u  = (E , H )`: container for starting field values
-   `params = (intE, intH, ext)`: mesh, dg, and material parameters
-   `t`: time to evaluate at

"""
function dg_maxwell!(uʰ, u, params, t)
    # unpack variables
    E  = u[1]
    H  = u[2]
    Eʰ = uʰ[1]
    Hʰ = uʰ[2]

    # unpack params
    int  = params[1] # internal parameters
    intH = params[2] # internal parameters
    ext  = params[3] # external parameters

    # compute impedence
    Z = @. sqrt(ext.μ / ext.ϵ)

    # define field differences at faces
    dE = similar(int.du)
    @. dE[:] = E[int.vmapM] - E[int.vmapP]
    dH = similar(dE)
    @. dH[:] = H[int.vmapM] - H[int.vmapP]

    # define impedances at the faces
    Z⁻ = similar(dE)
    @. Z⁻[:] = Z[int.vmapM]
    Z⁺ = similar(dE)
    @. Z⁺[:] = Z[int.vmapP]
    Y⁻ = similar(dE)
    @. Y⁻ = 1 / Z⁻
    Y⁺ = similar(dE)
    @. Y⁺ = 1 / Z⁺

    # homogenous boundary conditions, Ez = 0
    dE[int.mapB] = E[int.vmapB] + E[int.vmapB]
    dH[int.mapB] = H[int.vmapB] - H[int.vmapB]

    # evaluate upwind fluxes
    fluxE = @. 1/(Z⁻ + Z⁺) * (int.nx * Z⁻ * dH - dE)
    fluxH = @. 1/(Y⁻ + Y⁺) * (int.nx * Y⁻ * dE - dH)

    # compute right hand side of the PDE's
    mul!(Eʰ, int.D, H)
    @. Eʰ *= -int.rx
    liftE = int.lift * (int.fscale .* fluxE)
    @. Eʰ += liftE / ext.ϵ

    mul!(Hʰ, int.D, E)
    @. Hʰ *= -int.rx
    liftH = int.lift * (int.fscale .* fluxH)
    @. Hʰ += liftH / ext.μ

    return nothing
end
