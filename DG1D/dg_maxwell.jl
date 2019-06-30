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
-   `params = (𝒢, E, H, ext)`: mesh, E sol, H sol, and material parameters
-   `t`: time to evaluate at

"""
function dg_maxwell!(uʰ, u, params, t)
    # unpack params
    𝒢   = params[1] # grid parameters
    E   = params[2] # internal parameters for E
    H   = params[3] # internal parameters for H
    ext = params[4] # external parameters

    # unpack variables, no need to
    #=
    @. E.u  = u[1]
    @. H.u  = u[2]
    @. E.uʰ = uʰ[1]
    @. H.uʰ = uʰ[2]
    =#
    # compute impedence
    Z = @. sqrt(ext.μ / ext.ϵ)

    # define field differences at faces
    dE = similar(E.flux)
    @. dE[:] = E.u[𝒢.vmapM] - E.u[𝒢.vmapP]
    dH = similar(H.flux)
    @. dH[:] = H.u[𝒢.vmapM] - H.u[𝒢.vmapP]

    # define impedances at the faces
    Z⁻ = similar(dE)
    @. Z⁻[:] = Z[𝒢.vmapM]
    Z⁺ = similar(dE)
    @. Z⁺[:] = Z[𝒢.vmapP]
    Y⁻ = similar(dE)
    @. Y⁻ = 1 / Z⁻
    Y⁺ = similar(dE)
    @. Y⁺ = 1 / Z⁺

    # homogenous boundary conditions, Ez = 0
    dE[𝒢.mapB] = E.u[𝒢.vmapB] + E.u[𝒢.vmapB]
    dH[𝒢.mapB] = H.u[𝒢.vmapB] - H.u[𝒢.vmapB]

    # evaluate upwind fluxes
    @. E.flux = 1/(Z⁻ + Z⁺) * (𝒢.normals * Z⁻ * dH - dE)
    @. H.flux = 1/(Y⁻ + Y⁺) * (𝒢.normals * Y⁻ * dE - dH)

    # compute right hand side of the PDE's
    mul!(E.uʰ, 𝒢.D, H.u)
    @. E.uʰ *= -𝒢.rx
    liftE = 𝒢.lift * (𝒢.fscale .* E.flux)
    @. E.uʰ += liftE / ext.ϵ

    mul!(H.uʰ, 𝒢.D, E.u)
    @. H.uʰ *= -𝒢.rx
    liftH = 𝒢.lift * (𝒢.fscale .* H.flux)
    @. H.uʰ += liftH / ext.μ

    # pass values back into arguments, no need to
    #=
    @. u[1]  = E.u
    @. u[2]  = H.u
    @. uʰ[1] = E.uʰ
    @. uʰ[2] = H.uʰ
    =#
    return nothing
end
