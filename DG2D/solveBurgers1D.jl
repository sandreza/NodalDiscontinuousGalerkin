include("field2D.jl")
include("utils2D.jl")

"""
solveMaxwell!(u̇, u, params)

# Description

    numerical solution to 1D maxwell's equation

# Arguments

-   `u̇ = (Eʰ, Hʰ)`: container for numerical solutions to fields
-   `u  = (E , H )`: container for starting field values
-   `params = (𝒢, E, H, ext)`: mesh, E sol, H sol, and material parameters

"""
function solveBurgers1D!(fields, auxil, params, t)
    # unpack params
    𝒢 = params[1] # grid parameters
    ε = params[2]
    α = params[3]

    # unpack fields
    u  = fields[1]

    # auxiliary fields
    u² = auxil[1]
    uˣ = auxil[2]
    uʸ = auxil[3]

    for Ωᵏ in 𝒢.Ω
        # get volume nodes
        iⱽ = Ωᵏ.iⱽ

        # compute volume contribution to uˣ and uʸ
        ∇!(u.φˣ, u.φʸ, u.ϕ, Ωᵏ)
        @. uˣ.ϕ[iⱽ] = sqrt(ε) * u.φˣ[iⱽ]
        @. uʸ.ϕ[iⱽ] = sqrt(ε) * u.φʸ[iⱽ]

        # define physical fluxes for uˣ and uʸ
        @. uˣ.φˣ[iⱽ] = sqrt(ε) * u.ϕ[iⱽ]
        @. uʸ.φʸ[iⱽ] = sqrt(ε) * u.ϕ[iⱽ]

        # compute surface contributions to uˣ, uʸ
        for f in Ωᵏ.faces
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            computeCentralFluxes!(uˣ, f)
            computeCentralFluxes!(uʸ, f)

            # impose BC
            if f.isBoundary[1]
                uᴮ = [u⁰(𝒢.x[i,1],t) for i in f.i⁻]
                @. uˣ.fˣ[f.i⁻] = sqrt(ε) * uᴮ
                @. uʸ.fʸ[f.i⁻] = sqrt(ε) * uᴮ
            end

            computeSurfaceTerms!(uˣ, Ωᵏ, f)
            computeSurfaceTerms!(uʸ, Ωᵏ, f)
        end

        # compute u²
        @. u².ϕ[iⱽ] = u.ϕ[iⱽ]^2

        # define physical fluxes
        @. u.φˣ[iⱽ] = 0.5 * α * u².ϕ[iⱽ] - sqrt(ε) * uˣ.ϕ[iⱽ]
        @. u.φʸ[iⱽ] = 0.0 # make non-zero for 2D burgers eqn

        # compute volume contributions
        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ωᵏ)
        @. u.ϕ̇[iⱽ] = -u.𝚽[iⱽ]

        # compute surface contributions to tendency
        for f in Ωᵏ.faces
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
            C = maximum(abs.(u.ϕ[f.i⁻]))
            @. u.fˣ[f.i⁻] = 0.5 * α * u².ϕ°[f.i⁻] - sqrt(ε) * uˣ.ϕ°[f.i⁻]
            computeLaxFriedrichsFluxes!(u, f, C)
            @. u.fʸ[f.i⁻] = 0.0 # make non-zero for 2D burgers eqn

            computeSurfaceTerms!(u, Ωᵏ, f)
        end
    end

    return nothing
end
