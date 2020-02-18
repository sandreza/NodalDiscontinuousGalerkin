include("field2D.jl")
include("utils2D.jl")

"""
solveBurgers1D!(fields, auxils, params, t)

# Description

    numerical solution to Chorin Navier Stokes equation
    in vector form:
    ∂ᵗu = -∇(u²/2) + ν∇²u
    written out component wise for DG formulation:
    ∂ᵗu = -∂ˣ(u²/2 - νuˣ) - ∂ʸ(u²/2 - νuʸ)
    we are setting the flux in the y-direction to be zero for the 1D case


# Arguments

-   `fields = (u)`: velocity field
-   `auxils = (uˣ, uʸ, u²)`: auxiliary fields for computation
-   `params = (𝒢, ν, α, β)`: grid struct, viscosity, nonlinear switch, and 2D switch
-   `t`: time to compute BC at

"""
function solveBurgers1D!(fields, fluxes, auxils, params, t)
    # unpack params
    𝒢  = params[1] # grid parameters
    ν  = params[2]
    α  = params[3]
    β  = params[4]

    # unpack fields
    u  = fields[1]

    # auxiliary fields
    u² = auxils[1]
    uˣ = auxils[2]
    uʸ = auxils[3]

    φᵘ  = fluxes[1]
    φˣᵤ = fluxes[2]
    φʸᵤ = fluxes[3]

    # compute volume contribution to uˣ and uʸ
    for Ω in 𝒢.Ω
        computePhysicalFlux!(uˣ.φˣ, φᵘ, Ω)
        computePhysicalFlux!(uʸ.φʸ, φᵘ, Ω)

        # compute volume contributions
        ∇!(u.φˣ, u.φʸ, u.ϕ, Ω)
        @. uˣ.ϕ[Ω.iⱽ] = sqrt(ν) * u.φˣ[Ω.iⱽ]
        @. uʸ.ϕ[Ω.iⱽ] = sqrt(ν) * u.φʸ[Ω.iⱽ]
    end

    # compute surface contributions to uˣ, uʸ
    for Ω in 𝒢.Ω
        for f in Ω.faces
            computeCentralDifference!(u, f)

            # impose BC
            if f.isBoundary[1]
                uᴮ = [u⁰(𝒢.x[i,1],t) for i in f.i⁻]
                @. u.ϕ°[f.i⁻] = uᴮ
            end

            computeNumericalFlux!(uˣ.fˣ, φᵘ, f)
            computeNumericalFlux!(uʸ.fʸ, φᵘ, f)

            computeSurfaceTerms!(uˣ.ϕ, uˣ, Ω, f)
            computeSurfaceTerms!(uʸ.ϕ, uʸ, Ω, f)
        end
    end

    # compute u²
    @. u².ϕ = u.ϕ^2

    # compute volume contribution to tendency
    for Ω in 𝒢.Ω
        computePhysicalFlux!(u.φˣ, φˣᵤ, Ω)
        computePhysicalFlux!(u.φʸ, φʸᵤ, Ω)

        # compute volume contributions
        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ω)
        @. u.ϕ̇[Ω.iⱽ] = u.𝚽[Ω.iⱽ]
    end

    # compute surface contributions to tendency
    for Ω in 𝒢.Ω
        for f in Ω.faces
            computeCentralDifference!(uˣ, f)
            computeCentralDifference!(uʸ, f)
            computeCentralDifference!(u², f)

            # impose BC on uˣ, uʸ, and u²
            if f.isBoundary[1]
                uᴮ = [u⁰(𝒢.x[i,1],t) for i in f.i⁻]
                @. uˣ.ϕ°[f.i⁻] = uˣ.ϕ[f.i⁻]
                @. uʸ.ϕ°[f.i⁻] = uʸ.ϕ[f.i⁻]
                @. u².ϕ°[f.i⁻] = uᴮ^2
            end

            # evaluate numerical flux for u
            computeNumericalFlux!(u.fˣ, φˣᵤ, f)
            computeNumericalFlux!(u.fʸ, φʸᵤ, f)

            C = -maximum(abs.(u.ϕ[f.i⁻]))
            computeLaxFriedrichsFluxes!(u, f, C)

            computeSurfaceTerms!(u.ϕ̇, u, Ω, f)
        end
    end

    return nothing
end
