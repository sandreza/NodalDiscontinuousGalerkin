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
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            # evaluate numerical fluxes
            v⁻ = @. abs(f.nˣ * vˣ[i⁻] + f.nʸ * vʸ[i⁻])
            v⁺ = @. abs(f.nˣ * vˣ[i⁺] + f.nʸ * vʸ[i⁺])
            C = maximum([v⁻, v⁺])
            @. u.fˣ[i⁻] = 0.5 * (u.φˣ[i⁻] + u.φˣ[i⁺] + C * f.nˣ * (u.ϕ[i⁻] - u.ϕ[i⁺]))
            @. u.fʸ[i⁻] = 0.5 * (u.φʸ[i⁻] + u.φʸ[i⁺] + C * f.nʸ * (u.ϕ[i⁻] - u.ϕ[i⁺]))

            # impose BC
            if f.isBoundary[1]
                @. u.fˣ[i⁻] = u.φˣ[i⁻]
                @. u.fʸ[i⁻] = u.φʸ[i⁻]
            end

            # compute jump in flux
            @. u.Δf[i⁻] = f.nˣ * (u.φˣ[i⁻] - u.fˣ[i⁻]) + f.nʸ * (u.φʸ[i⁻] - u.fʸ[i⁻])

            # compute surface term
            u.∮f[iⱽ] = Ωᵏ.M⁺ * f.∮ * (f.C .* u.Δf[i⁻])
            @. u.ϕ̇[iⱽ] += u.∮f[iⱽ]
        end
    end

    @. U̇ = u.ϕ̇

    return nothing
end
