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
    𝑓 = params[end]

    @. 𝑓.ϕ = U

    # perform calculations over elements
    for Ωᵏ in 𝒢.Ω
        # get views of volume elements
        u  = view(𝑓.ϕ,  Ωᵏ.iⱽ)
        u̇  = view(𝑓.ϕ̇,  Ωᵏ.iⱽ)
        ∇u = view(𝑓.∇ϕ, Ωᵏ.iⱽ)

        # impose BC
        # @. u = 0.0

        # compute volume contributions
        ∇⨀!(∇u, vˣ[Ωᵏ.iⱽ] .* u, vʸ[Ωᵏ.iⱽ] .* u, Ωᵏ)
        @. u̇ = -∇u

        # compute surface contributions
        for f in Ωᵏ.faces
            # get views of surface elements
            Δu = view(𝑓.Δϕ, f.i⁻)
            fⁿ = view(𝑓.fⁿ, f.i⁻)

            # define field differences at faces
            @. Δu = 𝑓.ϕ[f.i⁻] - 𝑓.ϕ[f.i⁺]

            # evaluate flux
            vⁿ = @. f.nˣ * vˣ[f.i⁻] + f.nʸ * vʸ[f.i⁻]
            @. fⁿ = 1//2 * (vⁿ - α * abs(vⁿ)) * Δu

            # compute surface term
            ∮ᶠu = Ωᵏ.M⁺ * f.∮ * (f.C .* fⁿ)
            @. u̇ += ∮ᶠu
        end
    end

    @. U̇ = 𝑓.ϕ̇

    return nothing
end
