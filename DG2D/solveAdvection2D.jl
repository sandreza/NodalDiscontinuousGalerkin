include("field2D.jl")
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
function solveAdvection2D!(fields, params, t)
    # unpack params
    𝒢  = params[1]
    vˣ = params[2]
    vʸ = params[3]

    u = fields[1]

    # compute volume contributions
    for Ωᵏ in 𝒢.Ω
        # get volumes nodes
        iⱽ = Ωᵏ.iⱽ

        # define physical fluxes
        @. u.φˣ[iⱽ] = vˣ[iⱽ] .* u.ϕ[iⱽ]
        @. u.φʸ[iⱽ] = vʸ[iⱽ] .* u.ϕ[iⱽ]

        # compute volume contributions
        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ωᵏ)
        @. u.ϕ̇[iⱽ] = -u.𝚽[iⱽ]
    end

    # compute surface contributions
    for Ωᵏ in 𝒢.Ω
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

    return nothing
end
