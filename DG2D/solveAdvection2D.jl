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

        # compute volume contributions
        ∇⨀!(u.∇ϕ, vˣ .* u.ϕ, vʸ .* u.ϕ, Ωᵏ)
        @. u.ϕ̇[iⱽ] = -u.∇ϕ[iⱽ]

        # compute surface contributions
        for f in Ωᵏ.faces
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            # define field differences at faces
            @. u.Δϕ[i⁻] = u.ϕ[i⁻] - u.ϕ[i⁺]

            # impose BC
            if f.isBoundary[1]
                @. u.Δϕ[i⁻] = u.ϕ[i⁻]
            end

            # evaluate flux
            vⁿ = @. f.nˣ * vˣ[f.i⁻] + f.nʸ * vʸ[f.i⁻]
            @. u.fⁿ[i⁻] = 1//2 * (vⁿ - α * abs(vⁿ)) * u.Δϕ[i⁻]

            # compute surface term
            ∮ᶠu = Ωᵏ.M⁺ * f.∮ * (f.C .* u.fⁿ[i⁻])
            @. u.ϕ̇[iⱽ] += ∮ᶠu
        end
    end

    @. U̇ = u.ϕ̇

    return nothing
end
