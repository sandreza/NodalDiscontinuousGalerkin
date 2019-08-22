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
function solveAdvection2D!(U̇, U, params, t)
    # unpack params
    𝒢 = params[1] # grid parameters
    α = params[2]
    vˣ = params[3]
    vʸ = params[4]
    u = params[end]

    @. u.ϕ = U

    # perform calculations over elements
    for Ωᵏ in 𝒢.Ω
        # get volumes nodes
        iⱽ = Ωᵏ.iⱽ

        # define physical fluxes
        @. u.φˣ[iⱽ] = vˣ[iⱽ] .* u.ϕ[iⱽ]
        @. u.φʸ[iⱽ] = vʸ[iⱽ] .* u.ϕ[iⱽ]

        # compute volume contributions
        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ωᵏ)
        @. u.ϕ̇[iⱽ] = -u.𝚽[iⱽ]

        # compute surface contributions
        for f in Ωᵏ.faces
            # evaluate numerical fluxes
            v⁻ = @. abs(f.nˣ * vˣ[f.i⁻] + f.nʸ * vʸ[f.i⁻])
            v⁺ = @. abs(f.nˣ * vˣ[f.i⁺] + f.nʸ * vʸ[f.i⁺])
            C = maximum([v⁻, v⁺])
            computeCentralFluxes!(u, f)
            computeLaxFriedrichsFluxes!(u, f, C)

            # impose BC
            if f.isBoundary[1]
                @. u.fˣ[i⁻] = u.φˣ[i⁻]
                @. u.fʸ[i⁻] = u.φʸ[i⁻]
            end

            computeSurfaceTerms!(u, Ωᵏ, f)
        end
    end

    @. U̇ = u.ϕ̇

    return nothing
end
