include("field2D.jl")
include("utils2D.jl")

"""
solveMaxwell!(fields, params)

# Description

    numerical solution to 1D maxwell's equation

# Arguments

-   `fields = (Hˣ, Hʸ, Eᶻ)`: fields to compute
-   `params = (𝒢, α)`: parameters needed for computation

"""
function solveMaxwell2D!(fields, params)
    # unpack params
    𝒢 = params[1] # grid parameters
    α = params[2]

    # unpack fields
    Hˣ = fields[1]
    Hʸ = fields[2]
    Eᶻ = fields[3]

    # perform calculations over elements
    for Ωᵏ in 𝒢.Ω
        # get volume nodes
        iⱽ = Ωᵏ.iⱽ

        # compute volume contributions
        ∇!(Hʸ.∇ϕ, Hˣ.∇ϕ, Eᶻ.ϕ, Ωᵏ)
        ∇⨂!(Eᶻ.∇ϕ, Hˣ.ϕ, Hʸ.ϕ, Ωᵏ)

        @. Hˣ.ϕ̇[iⱽ] = -Hˣ.∇ϕ[iⱽ]
        @. Hʸ.ϕ̇[iⱽ] =  Hʸ.∇ϕ[iⱽ]
        @. Eᶻ.ϕ̇[iⱽ] =  Eᶻ.∇ϕ[iⱽ]

        # compute surface contributions
        for f in Ωᵏ.faces
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            # define field differences at faces
            @. Hˣ.Δϕ[i⁻] = Hˣ.ϕ[i⁻] - Hˣ.ϕ[i⁺]
            @. Hʸ.Δϕ[i⁻] = Hʸ.ϕ[i⁻] - Hʸ.ϕ[i⁺]
            @. Eᶻ.Δϕ[i⁻] = Eᶻ.ϕ[i⁻] - Eᶻ.ϕ[i⁺]

            # impose reflective BC
            if f.isBoundary[1]
                @. Hˣ.Δϕ[i⁻] = 0
                @. Hʸ.Δϕ[i⁻] = 0
                @. Eᶻ.Δϕ[i⁻] = 2 * Eᶻ.ϕ[i⁻]
            end

            # evaluate upwind fluxes
            nˣΔH = @. f.nˣ * (f.nˣ * Hˣ.Δϕ[i⁻] + f.nʸ * Hʸ.Δϕ[i⁻])
            nʸΔH = @. f.nʸ * (f.nˣ * Hˣ.Δϕ[i⁻] + f.nʸ * Hʸ.Δϕ[i⁻])

            # minus isn't defined for these fluxes?????
            @. Hˣ.fⁿ[i⁻] =  f.nʸ * Eᶻ.Δϕ[i⁻] + α * (nˣΔH - Hˣ.Δϕ[i⁻])
            @. Hʸ.fⁿ[i⁻] = -f.nˣ * Eᶻ.Δϕ[i⁻] + α * (nʸΔH - Hʸ.Δϕ[i⁻])
            @. Eᶻ.fⁿ[i⁻] = -f.nˣ * Hʸ.Δϕ[i⁻] + f.nʸ * Hˣ.Δϕ[i⁻] - α * Eᶻ.Δϕ[i⁻]

            # compute RHS of PDE's
            ∮Hˣ = 1//2 * Ωᵏ.M⁺ * f.∮ * (f.C .* Hˣ.fⁿ[i⁻])
            ∮Hʸ = 1//2 * Ωᵏ.M⁺ * f.∮ * (f.C .* Hʸ.fⁿ[i⁻])
            ∮Eᶻ = 1//2 * Ωᵏ.M⁺ * f.∮ * (f.C .* Eᶻ.fⁿ[i⁻])

            @. Hˣ.ϕ̇[iⱽ] += ∮Hˣ
            @. Hʸ.ϕ̇[iⱽ] += ∮Hʸ
            @. Eᶻ.ϕ̇[iⱽ] += ∮Eᶻ
        end
    end

    return nothing
end
