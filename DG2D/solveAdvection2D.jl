include("field2D.jl")
include("flux2D.jl")
include("utils2D.jl")

"""
solveAdvection2D!(fields, params, t)

# Description

    numerical solution to Chorin Navier Stokes equation
    in vector form:
    ∂ᵗθ = -∇⋅(ṽθ)
    written out component wise for DG formulation:
    ∂ᵗθ = -∂ˣ(vˣ * θ) - ∂ʸ(vʸ * θ)

# Arguments

-   `fields = (θ)`: velocity field
-   `params = (𝒢, vˣ, vʸ)`: grid struct and velocities in each direction
-   `t`: time to compute BC at

"""
function solveAdvection2D!(fields, fluxes, auxils, params, t)
    # unpack params
    𝒢  = params[1]
    vˣ = params[2]
    vʸ = params[3]

    θ  = fields[1]
    θˣ = auxils[1]
    θʸ = auxils[2]

    φˣ = fluxes[1]
    φʸ = fluxes[2]

    # define physical fluxes
    @. θˣ.ϕ = vˣ .* θ.ϕ
    @. θʸ.ϕ = vʸ .* θ.ϕ

    # compute volume contributions
    for Ω in 𝒢.Ω
        computePhysicalFlux!(θ.φˣ, φˣ, Ω)
        computePhysicalFlux!(θ.φʸ, φʸ, Ω)

        # compute volume contributions
        ∇⨀!(θ.𝚽, θ.φˣ, θ.φʸ, Ω)
        @. θ.ϕ̇[Ω.iⱽ] = θ.𝚽[Ω.iⱽ]
    end

    # compute surface contributions
    for Ω in 𝒢.Ω
        for f in Ω.faces
            computeCentralDifference!(θˣ, f)
            computeCentralDifference!(θʸ, f)

            # impose BC
            if f.isBoundary[1]
                @. θˣ.ϕ°[f.i⁻] = 0. # -θˣ.ϕ[f.i⁻]
                @. θʸ.ϕ°[f.i⁻] = 0. # -θʸ.ϕ[f.i⁻]
            end

            computeNumericalFlux!(θ.fˣ, φˣ, f)
            computeNumericalFlux!(θ.fʸ, φʸ, f)

            v⁻ = @. abs(f.nˣ * vˣ[f.i⁻] + f.nʸ * vʸ[f.i⁻])
            v⁺ = @. abs(f.nˣ * vˣ[f.i⁺] + f.nʸ * vʸ[f.i⁺])
            C = -maximum([v⁻, v⁺])
            # computeLaxFriedrichsFluxes!(θ, f, C)

            computeSurfaceTerms!(θ.ϕ̇, θ, Ω, f)
        end
    end

    return nothing
end
