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

    # define field differences at faces
    @. Hˣ.Δϕ = Hˣ.ϕ[𝒢.nodes⁻] - Hˣ.ϕ[𝒢.nodes⁺]
    @. Hʸ.Δϕ = Hʸ.ϕ[𝒢.nodes⁻] - Hʸ.ϕ[𝒢.nodes⁺]
    @. Eᶻ.Δϕ = Eᶻ.ϕ[𝒢.nodes⁻] - Eᶻ.ϕ[𝒢.nodes⁺]

    # impose reflective BC
    @. Hˣ.Δϕ[𝒢.mapᴮ] = 0
    @. Hʸ.Δϕ[𝒢.mapᴮ] = 0
    @. Eᶻ.Δϕ[𝒢.mapᴮ] = 2 * Eᶻ.ϕ[𝒢.nodesᴮ]

    # perform calculations over elements
    let nGL = nBP = 0
        for k in 1:𝒢.ℳ.K
            # get element and number of GL points
            Ωᵏ = 𝒢.Ω[k]
            nGLᵏ = (nGL + 1):(nGL + Ωᵏ.nGL)
            nBPᵏ = (nBP + 1):(nBP + Ωᵏ.nBP)
            nGL += Ωᵏ.nGL
            nBP += Ωᵏ.nBP

            # get views of computation elements
            uHˣ = view(Hˣ.ϕ, nGLᵏ)
            uHʸ = view(Hʸ.ϕ, nGLᵏ)
            uEᶻ = view(Eᶻ.ϕ, nGLᵏ)

            ϕ̇Hˣ = view(Hˣ.ϕ̇, nGLᵏ)
            ϕ̇Hʸ = view(Hʸ.ϕ̇, nGLᵏ)
            ϕ̇Eᶻ = view(Eᶻ.ϕ̇, nGLᵏ)

            ∇Hˣ = view(Hˣ.∇ϕ, nGLᵏ)
            ∇Hʸ = view(Hʸ.∇ϕ, nGLᵏ)
            ∇Eᶻ = view(Eᶻ.∇ϕ, nGLᵏ)

            ΔHˣ = view(Hˣ.Δϕ, nBPᵏ)
            ΔHʸ = view(Hʸ.Δϕ, nBPᵏ)
            ΔEᶻ = view(Eᶻ.Δϕ, nBPᵏ)

            fHˣ = view(Hˣ.fⁿ, nBPᵏ)
            fHʸ = view(Hʸ.fⁿ, nBPᵏ)
            fEᶻ = view(Eᶻ.fⁿ, nBPᵏ)

            # evaluate upwind fluxes
            nˣΔH = @. Ωᵏ.nˣ * (Ωᵏ.nˣ * ΔHˣ + Ωᵏ.nʸ * ΔHʸ)
            nʸΔH = @. Ωᵏ.nʸ * (Ωᵏ.nˣ * ΔHˣ + Ωᵏ.nʸ * ΔHʸ)

            # minus isn't defined for these fluxes?????
            @. fHˣ =      Ωᵏ.nʸ * ΔEᶻ + α * (nˣΔH + (-1 * ΔHˣ))
            @. fHʸ = -1 * Ωᵏ.nˣ * ΔEᶻ + α * (nʸΔH + (-1 * ΔHʸ))
            @. fEᶻ = -1 * Ωᵏ.nˣ * ΔHʸ + Ωᵏ.nʸ * ΔHˣ + (-1 * α * ΔEᶻ)

            # local derivatives of the fields
            ∇!(∇Hʸ, ∇Hˣ, uEᶻ, Ωᵏ)
            ∇⨂!(∇Eᶻ, uHˣ, uHʸ, Ωᵏ)

            # compute RHS of PDE's
            liftHˣ = 1//2 * Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* fHˣ)
            liftHʸ = 1//2 * Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* fHʸ)
            liftEᶻ = 1//2 * Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* fEᶻ)

            @. ϕ̇Hˣ = -∇Hˣ + liftHˣ
            @. ϕ̇Hʸ =  ∇Hʸ + liftHʸ
            @. ϕ̇Eᶻ =  ∇Eᶻ + liftEᶻ
        end
    end

    return nothing
end
