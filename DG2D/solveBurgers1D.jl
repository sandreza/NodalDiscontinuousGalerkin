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
    𝑓ᵘ = fields[1]
    𝑓² = fields[2]
    𝑓ᵖ = fields[3]

    # define field differences at faces
    @. 𝑓ᵘ.Δϕ = 𝑓ᵘ.ϕ[𝒢.nodes⁻] - 𝑓ᵘ.ϕ[𝒢.nodes⁺]
    @. 𝑓².Δϕ = 1//2 * (𝑓ᵘ.ϕ[𝒢.nodes⁻]^2 - 𝑓ᵘ.ϕ[𝒢.nodes⁺]^2)

    # impose Dirichlet BC on u
    # @. 𝑓ᵘ.ϕ[𝒢.mapᴮ] = 2 * (𝑓ᵘ.ϕ[𝒢.nodesᴮ] - u⁰(𝒢.x[1]))
    # @. 𝑓².ϕ[𝒢.mapᴮ] = 𝑓ᵘ.ϕ[𝒢.nodesᴮ]^2 - u⁰(𝒢.x[1])^2

    # calculate max value of u (might need to be a face by face calculation later)
    maxu = maximum(abs.(𝑓ᵘ.ϕ))

    # calculate q
    let nGL = nBP = 0
        for Ωᵏ in 𝒢.Ω
            # get number of GL points
            GLᵏ  = (nGL + 1):(nGL + Ωᵏ.nGL)
            BPᵏ  = (nBP + 1):(nBP + Ωᵏ.nBP)
            nGL += Ωᵏ.nGL
            nBP += Ωᵏ.nBP

            # get views of computation elements
            u  = view(𝑓ᵘ.ϕ,  GLᵏ)
            uˣ = view(𝑓ᵘ.φˣ, GLᵏ)
            uʸ = view(𝑓ᵘ.φʸ, GLᵏ)
            Δu = view(𝑓ᵘ.Δϕ, BPᵏ)

            q  = view(𝑓ᵖ.ϕ,  GLᵏ)

            # interior terms
            ∇!(uˣ, uʸ, u, Ωᵏ)

            # surface terms
            ∮ˣu = 1//2 * Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nˣ .* Δu)

            # combine them
            @. q = sqrt(ε) * uˣ - ∮ˣu
        end
    end

    # define field differences at faces
    @. 𝑓ᵖ.Δϕ = 1//2 * (𝑓ᵖ.ϕ[𝒢.nodes⁻] - 𝑓ᵖ.ϕ[𝒢.nodes⁺])

    # impose Dirichlet BC on q
    @. 𝑓ᵖ.Δϕ[𝒢.mapᴮ] = 0.0

    # perform calculations over elements
    let nGL = nBP = 0
        for Ωᵏ in 𝒢.Ω
            # get number of GL points
            GLᵏ  = (nGL + 1):(nGL + Ωᵏ.nGL)
            BPᵏ  = (nBP + 1):(nBP + Ωᵏ.nBP)
            nGL += Ωᵏ.nGL
            nBP += Ωᵏ.nBP

            # get views of computation elements
            u   = view(𝑓ᵘ.ϕ,  GLᵏ)
            u̇   = view(𝑓ᵘ.ϕ̇,  GLᵏ)
            ∇u  = view(𝑓ᵘ.∇ϕ, GLᵏ)
            Δu  = view(𝑓ᵘ.Δϕ, BPᵏ)
            uˣ  = view(𝑓ᵘ.φˣ, GLᵏ)
            uʸ  = view(𝑓ᵘ.φʸ, GLᵏ)
            fⁿ  = view(𝑓ᵘ.fⁿ, BPᵏ)

            q   = view(𝑓ᵖ.ϕ,  GLᵏ)
            Δq  = view(𝑓ᵖ.Δϕ, BPᵏ)

            Δu² = view(𝑓².Δϕ, BPᵏ)

            # evaluate numerical flux
            @. fⁿ = Ωᵏ.nˣ * (α * Δu²/2 - sqrt(ε) * Δq) - 1//2 * maxu * Δu

            # compute surface term
            ∮u = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* fⁿ)

            # define physical flux in the x direction
            @. ∇u = 1//2 * α * u^2 - sqrt(ε) * q

            # define derivatives of physical flux
            ∇!(uˣ, uʸ, ∇u, Ωᵏ)
            @. ∇u = uˣ

            # combine terms
            @. u̇ = -∇u + ∮u
        end
    end

    return nothing
end
