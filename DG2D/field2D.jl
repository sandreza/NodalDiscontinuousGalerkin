include("grid2D.jl")

abstract type AbstractField2D end
abstract type AbstractAuxiliaryField2D <: AbstractField2D end

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
    # field value and tendency
    ϕ::T
    ϕ̇::T

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

    # residual for RK4 methods
    r::T

    function Field2D(𝒢::Grid2D)
        ϕ  = zeros(𝒢.nGL)
        ϕ̇  = zeros(𝒢.nGL)

        𝚽  = zeros(𝒢.nGL)

        φˣ = zeros(𝒢.nGL)
        φʸ = zeros(𝒢.nGL)

        fˣ = zeros(𝒢.nGL)
        fʸ = zeros(𝒢.nGL)

        Δf = zeros(𝒢.nGL)
        ∮f = zeros(𝒢.nGL)

        r  = zeros(𝒢.nGL)

        return new{typeof(ϕ)}(ϕ,ϕ̇, 𝚽, φˣ,φʸ, fˣ,fʸ, Δf,∮f, r)
    end
end

"""
AuxiliaryField2D(𝒢::Grid2D)

# Description

    Contains all the computational elements necessary to compute an auxiliary field

# Arguments

-   `mesh`: a mesh to compute on

# Return Values:

-   `u` : the field to be computed
-   `u̇`: numerical solutions for the field
-   `flux`: the numerical flux for the computation

"""
struct AuxiliaryField2D{T} <: AbstractAuxiliaryField2D
    # physical value and numerical value
    ϕ::T
    ϕ°::T

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

    function AuxiliaryField2D(𝒢::Grid2D)
        ϕ  = zeros(𝒢.nGL)
        ϕ° = zeros(𝒢.nGL)

        φˣ = zeros(𝒢.nGL)
        φʸ = zeros(𝒢.nGL)

        fˣ = zeros(𝒢.nGL)
        fʸ = zeros(𝒢.nGL)

        Δf = zeros(𝒢.nGL)
        ∮f = zeros(𝒢.nGL)

        return new{typeof(ϕ)}(ϕ,ϕ°, φˣ,φʸ, fˣ,fʸ, Δf,∮f)
    end
end

function computeCentralFluxes!(𝑓::Field2D, f::Face2D)
    @. 𝑓.fˣ[f.i⁻] = 0.5 * (𝑓.φˣ[f.i⁻] + 𝑓.φˣ[f.i⁺])
    @. 𝑓.fʸ[f.i⁻] = 0.5 * (𝑓.φʸ[f.i⁻] + 𝑓.φʸ[f.i⁺])

    return nothing
end

function computeLaxFriedrichsFluxes!(𝑓::Field2D, f::Face2D)
    C = maximum(abs.([𝑓.ϕ[f.i⁻]; 𝑓.ϕ[f.i⁺]]))
    @. 𝑓.fˣ[f.i⁻] += C * f.nˣ * (𝑓.ϕ[f.i⁻] - 𝑓.ϕ[f.i⁺])
    @. 𝑓.fʸ[f.i⁻] += C * f.nˣ * (𝑓.ϕ[f.i⁻] - 𝑓.ϕ[f.i⁺])

    return nothing
end

function computeCentralDifference!(𝑓::AuxiliaryField2D, f::Face2D)
    @. 𝑓.ϕ°[f.i⁻] = 0.5 * (𝑓.ϕ[f.i⁻] + 𝑓.ϕ[f.i⁺])

    return nothing
end

function computeSurfaceTerms!(𝑓::Field2D, f::Face2D)
    # compute jump in flux
    @. 𝑓.Δf[f.i⁻] = f.nˣ * (𝑓.φˣ[f.i⁻] - 𝑓.fˣ[f.i⁻]) + f.nʸ * (𝑓.φʸ[f.i⁻] - 𝑓.fʸ[f.i⁻])

    # compute surface terms
    𝑓.∮f[f.iⱽ] = Ωᵏ.M⁺ * f.∮ * (f.C .* 𝑓.Δf[f.i⁻])
    @. 𝑓.ϕ̇[f.iⱽ] += 𝑓.∮f[f.iⱽ]

    return nothing
end

function computeSurfaceTerms!(𝑓::AuxiliaryField2D, f::Face2D)
    # compute jump in flux
    @. 𝑓.Δf[f.i⁻] = f.nˣ * (𝑓.φˣ[f.i⁻] - 𝑓.fˣ[f.i⁻]) + f.nʸ * (𝑓.φʸ[f.i⁻] - 𝑓.fʸ[f.i⁻])

    # compute surface terms
    𝑓.∮f[f.iⱽ] = Ωᵏ.M⁺ * f.∮ * (f.C .* 𝑓.Δf[f.i⁻])
    @. 𝑓.ϕ[f.iⱽ] -= 𝑓.∮f[f.iⱽ]

    return nothing
end
