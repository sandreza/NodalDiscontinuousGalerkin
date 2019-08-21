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
    q  = fields[3]

    # calculate q
    for Ωᵏ in 𝒢.Ω
        # get volume nodes
        iⱽ = Ωᵏ.iⱽ

        # compute volume contribution to q
        ∇!(u.φˣ, u.φʸ, u.ϕ, Ωᵏ)
        @. q.ϕ[iⱽ] = sqrt(ε) * u.φˣ[iⱽ]

        # compute u²
        @. u².ϕ[iⱽ] = u.ϕ[iⱽ]^2

        # compute surface contributions to q
        for f in Ωᵏ.faces
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            # define field differences at faces
            @. u.Δϕ[i⁻] = 1//2 * (u.ϕ[i⁻] - u.ϕ[i⁺])

            # impose Dirichlet BC on u
            if f.isBoundary[1]
                @. u.Δϕ[i⁻] = 2 * (u.ϕ[i⁻] -  u⁰(𝒢.x[1]))
            end

            # compute surface terms
            ∮ˣu = Ωᵏ.M⁺ * f.∮ * (f.C .* f.nˣ .* u.Δϕ[i⁻])

            # combine them
            @. q.ϕ[iⱽ] -= ∮ˣu
        end

        # define physical flux
        @. u.∇ϕ[iⱽ] = 1//2 * α * u².ϕ[iⱽ] - sqrt(ε) * q.ϕ[iⱽ]

        # compute volume contributions to tendency
        ∇!(u.φˣ, u.φʸ, u.∇ϕ, Ωᵏ)
        @. u.ϕ̇[iⱽ] = -u.φˣ[iⱽ]

        # compute surface contributions to tendency
        for f in Ωᵏ.faces
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            # define field differences at faces
            @.  q.Δϕ[i⁻] = 1//2 * ( q.ϕ[i⁻] -  q.ϕ[i⁺])
            @. u².Δϕ[i⁻] = 1//2 * (u².ϕ[i⁻] - u².ϕ[i⁺])

            # impose BC on q and u²
            if f.isBoundary[1]
                @.  q.Δϕ[i⁻] = 0.0
                @. u².Δϕ[i⁻] = u².ϕ[i⁻] - u⁰(𝒢.x[1])^2
            end

            # evaluate numerical flux
            maxu = maximum(abs.(u.ϕ[i⁻]))
            @. u.fⁿ[i⁻] = f.nˣ * (1//2 * α * u².Δϕ[i⁻] - sqrt(ε) * q.Δϕ[i⁻]) - 1//2 * maxu * u.Δϕ[i⁻]

            # compute surface term
            ∮ᶠu = Ωᵏ.M⁺ * f.∮ * (f.C .* u.fⁿ[i⁻])

            # combine terms
            @. u.ϕ̇[iⱽ] += ∮ᶠu
        end
    end

    return nothing
end
