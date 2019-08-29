include("field2D.jl")
include("flux2D.jl")
include("utils2D.jl")

"""
solveAdvection2D!(fields, params, t)

# Description

    numerical solution to Chorin Navier Stokes equation
    in vector form:
    ∂ᵗu = -∇⋅(ṽu)
    written out component wise for DG formulation:
    ∂ᵗu = -∂ˣ(vˣ * u) - ∂ʸ(vʸ * u)

# Arguments

-   `fields = (u)`: velocity field
-   `params = (𝒢, vˣ, vʸ)`: grid struct and velocities in each direction
-   `t`: time to compute BC at

"""
function solveAdvection2D!(fields, fluxes, auxils, params, t)
    # unpack params
    𝒢  = params[1]
    vˣ = params[2]
    vʸ = params[3]

    u  = fields[1]
    θˣ = auxils[1]
    θʸ = auxils[2]

    φˣ = fluxes[1]
    φʸ = fluxes[2]

    # define physical fluxes
    @. θˣ.ϕ = vˣ .* u.ϕ
    @. θʸ.ϕ = vʸ .* u.ϕ

    # compute volume contributions
    for Ω in 𝒢.Ω
        computePhysicalFlux!(u.φˣ, φˣ, Ω)
        computePhysicalFlux!(u.φʸ, φʸ, Ω)

        # compute volume contributions
        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ω)
        @. u.ϕ̇[Ω.iⱽ] = u.𝚽[Ω.iⱽ]
    end

    # compute surface contributions
    for Ω in 𝒢.Ω
        for f in Ω.faces
            computeCentralDifference!(θˣ, f)
            computeCentralDifference!(θʸ, f)

            # impose BC
            if f.isBoundary[1]
                @. θˣ.ϕ°[f.i⁻] = θˣ.ϕ[f.i⁻]
                @. θʸ.ϕ°[f.i⁻] = θʸ.ϕ[f.i⁻]
            end

            computeNumericalFlux!(u.fˣ, φˣ, f)
            computeNumericalFlux!(u.fʸ, φʸ, f)

            v⁻ = @. abs(f.nˣ * vˣ[f.i⁻] + f.nʸ * vʸ[f.i⁻])
            v⁺ = @. abs(f.nˣ * vˣ[f.i⁺] + f.nʸ * vʸ[f.i⁺])
            C = -maximum([v⁻, v⁺])
            computeLaxFriedrichsFluxes!(u, f, C)

            computeSurfaceTerms!(u.ϕ̇, u, Ω, f)
        end
    end

    return nothing
end
