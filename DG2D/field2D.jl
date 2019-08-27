include("grid2D.jl")

abstract type AbstractField2D end

"""
Field2D(𝒢::Grid2D)

# Description

    Contains all the computational elements necessary to evolve a field in time

# Arguments

-   `mesh`: a mesh to compute on

# Return Values:

-   `u` : the field to be computed
-   `u̇`: numerical solutions for the field
-   `flux`: the numerical flux for the computation

"""
struct Field2D{T} <: AbstractField2D
    # field value and numerical value
    ϕ::T
    ϕ°::T

    # volume contributions to tendency
    𝚽::T

    # physical fluxes
    φˣ::T
    φʸ::T

    # numerical fluxes
    fˣ::T
    fʸ::T

    # jump in the flux
    Δf::T

    # surface contributions to the tendency
    ∮f::T

    # tendency and residual for RK4 methods
    ϕ̇::T
    r::T

    function Field2D(𝒢::Grid2D)
        ϕ  = zeros(𝒢.nGL)
        ϕ° = zeros(𝒢.nGL)

        𝚽  = zeros(𝒢.nGL)

        φˣ = zeros(𝒢.nGL)
        φʸ = zeros(𝒢.nGL)

        fˣ = zeros(𝒢.nGL)
        fʸ = zeros(𝒢.nGL)

        Δf = zeros(𝒢.nGL)
        ∮f = zeros(𝒢.nGL)

        ϕ̇  = zeros(𝒢.nGL)
        r  = zeros(𝒢.nGL)

        return new{typeof(ϕ)}(ϕ,ϕ°, 𝚽, φˣ,φʸ, fˣ,fʸ, Δf,∮f, ϕ̇,r)
    end
end

function computeCentralDifference!(𝑓::Field2D, f::Face2D)
    @. 𝑓.ϕ°[f.i⁻] = 0.5 * (𝑓.ϕ[f.i⁻] + 𝑓.ϕ[f.i⁺])

    return nothing
end

function computeLaxFriedrichsFluxes!(𝑓::Field2D, f::Face2D, C)
    @. 𝑓.fˣ[f.i⁻] += 0.5 * C * f.nˣ * (𝑓.ϕ[f.i⁻] - 𝑓.ϕ[f.i⁺])
    @. 𝑓.fʸ[f.i⁻] += 0.5 * C * f.nʸ * (𝑓.ϕ[f.i⁻] - 𝑓.ϕ[f.i⁺])

    return nothing
end

function computeSurfaceTerms!(ϕ, 𝑓::Field2D, Ωᵏ::Element2D, f::Face2D)
    # compute jump in flux
    @. 𝑓.Δf[f.i⁻] = f.nˣ * (𝑓.φˣ[f.i⁻] - 𝑓.fˣ[f.i⁻]) + f.nʸ * (𝑓.φʸ[f.i⁻] - 𝑓.fʸ[f.i⁻])

    # compute surface terms
    𝑓.∮f[Ωᵏ.iⱽ] = Ωᵏ.M⁺ * f.∮ * (f.C .* 𝑓.Δf[f.i⁻])
    @. ϕ[Ωᵏ.iⱽ] -= 𝑓.∮f[Ωᵏ.iⱽ]

    return nothing
end
