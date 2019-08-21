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
function solveBurgers1D!(fields, params)
    # unpack params
    𝒢 = params[1] # grid parameters
    ε = params[2]

    # unpack fields
    u  = fields[1]
    u² = fields[2]
    uˣ = fields[3]
    uʸ = fields[4]

    for Ωᵏ in 𝒢.Ω
        # get volume nodes
        iⱽ = Ωᵏ.iⱽ

        # compute volume contribution to uˣ and uʸ
        ∇!(u.φˣ, u.φʸ, u.ϕ, Ωᵏ)
        @. uˣ.ϕ[iⱽ] = sqrt(ε) * u.φˣ[iⱽ]
        @. uʸ.ϕ[iⱽ] = sqrt(ε) * u.φʸ[iⱽ]

        # define physical fluxes for uˣ and uʸ
        @. uˣ.φˣ[iⱽ] = sqrt(ε) * u.ϕ[iⱽ]
        @. uˣ.φʸ[iⱽ] = 0.0
        @. uʸ.φˣ[iⱽ] = 0.0
        @. uʸ.φʸ[iⱽ] = sqrt(ε) * u.ϕ[iⱽ]

        # compute surface contributions to uˣ, uʸ
        for f in Ωᵏ.faces
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            # evaluate numerical fluxes
            @. uˣ.fˣ[i⁻] = 0.5 * (uˣ.φˣ[i⁻] + uˣ.φˣ[i⁺])
            @. uˣ.fʸ[i⁻] = 0.0
            @. uʸ.fˣ[i⁻] = 0.0
            @. uʸ.fʸ[i⁻] = 0.5 * (uʸ.φʸ[i⁻] + uʸ.φʸ[i⁺])

            # impose BC
            if f.isBoundary[1]
                @. uˣ.fˣ[i⁻] = 2 * (u.ϕ[i⁻] -  u⁰(𝒢.x[1]))
                @. uʸ.fʸ[i⁻] = 2 * (u.ϕ[i⁻] -  u⁰(𝒢.x[1]))
            end

            # compute jumps in flux
            @. uˣ.Δf[i⁻] = f.nˣ * (uˣ.φˣ[i⁻] - uˣ.fˣ[i⁻]) + f.nʸ * (uˣ.φʸ[i⁻] - uˣ.fʸ[i⁻])
            @. uʸ.Δf[i⁻] = f.nˣ * (uʸ.φˣ[i⁻] - uʸ.fˣ[i⁻]) + f.nʸ * (uʸ.φʸ[i⁻] - uʸ.fʸ[i⁻])

            # compute surface terms
            uˣ.∮f[iⱽ] = Ωᵏ.M⁺ * f.∮ * (f.C .* uˣ.Δf[i⁻])
            uʸ.∮f[iⱽ] = Ωᵏ.M⁺ * f.∮ * (f.C .* uʸ.Δf[i⁻])
            @. uˣ.ϕ[iⱽ] -= uˣ.∮f[iⱽ]
            @. uʸ.ϕ[iⱽ] -= uʸ.∮f[iⱽ]
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
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            # evaluate numerical fluxes
            @. uˣ.fˣ[i⁻] = 0.5 * (uˣ.ϕ[i⁻] + uˣ.ϕ[i⁺])
            @. uʸ.fʸ[i⁻] = 0.5 * (uʸ.ϕ[i⁻] + uʸ.ϕ[i⁺])
            @. u².fˣ[i⁻] = 0.5 * (u².ϕ[i⁻] + u².ϕ[i⁺])
            @. u².fʸ[i⁻] = 0.5 * (u².ϕ[i⁻] + u².ϕ[i⁺])

            # impose BC on uˣ, uʸ, and u²
            if f.isBoundary[1]
                @. uˣ.fˣ[i⁻] = 0.0
                @. uʸ.fʸ[i⁻] = 0.0
                @. u².fˣ[i⁻] = u².ϕ[i⁻] - u⁰(𝒢.x[1])^2
                @. u².fʸ[i⁻] = u².ϕ[i⁻] - u⁰(𝒢.x[1])^2
            end

            # evaluate numerical flux for u
            C = maximum(abs.(u.ϕ[i⁻]))
            @. u.fˣ[i⁻] = 0.5 * α * u².fˣ[i⁻] - sqrt(ε) * uˣ.fˣ[i⁻] + 0.5 * C * f.nˣ * (u.ϕ[i⁻] - u.ϕ[i⁺])
            @. u.fʸ[i⁻] = 0.0 # make non-zero for 2D burgers eqn

            # compute jump in flux
            @. u.Δf[i⁻] = f.nˣ * (u.φˣ[i⁻] - u.fˣ[i⁻]) + f.nʸ * (u.φʸ[i⁻] - u.fʸ[i⁻])

            # compute surface term
            u.∮f[iⱽ] = Ωᵏ.M⁺ * f.∮ * (f.C .* u.Δf[i⁻])
            @. u.ϕ̇[iⱽ] += u.∮f[iⱽ]
        end
    end

    return nothing
end
