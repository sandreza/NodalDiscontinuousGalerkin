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
    𝑓 = params[end]

    @. 𝑓.ϕ = U

    # define field differences at faces
    @. 𝑓.Δϕ = 𝑓.ϕ[𝒢.nodes⁻] - 𝑓.ϕ[𝒢.nodes⁺]

    # impose BC
    # @. 𝑓.ϕ[𝒢.nodesᴮ] = 0.0

    # perform calculations over elements
    let nGL = nBP = 0
        for Ωᵏ in 𝒢.Ω
            # get number of GL points
            GLᵏ  = (nGL + 1):(nGL + Ωᵏ.nGL)
            BPᵏ  = (nBP + 1):(nBP + Ωᵏ.nBP)
            nGL += Ωᵏ.nGL
            nBP += Ωᵏ.nBP

            # get views of params
            vˣ = view(params[3], GLᵏ)
            vʸ = view(params[4], GLᵏ)

            # get views of computation elements
            u  = view(𝑓.ϕ,  GLᵏ)
            u̇  = view(𝑓.ϕ̇,  GLᵏ)
            ∇u = view(𝑓.∇ϕ, GLᵏ)
            Δu = view(𝑓.Δϕ, BPᵏ)
            f  = view(𝑓.fⁿ, BPᵏ)

            # local derivatives of the fields
            ∇⨀!(∇u, vˣ .* u, vʸ .* u, Ωᵏ)

            # evaluate flux
            vⁿ = @. Ωᵏ.nˣ * vˣ[Ωᵏ.fmask][:] + Ωᵏ.nʸ * vʸ[Ωᵏ.fmask][:]
            @. f = 1//2 * (vⁿ - α * abs(vⁿ)) * Δu

            # compute surface term
            lift = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* f)

            # compute RHS of PDE's
            @. u̇ = -∇u + lift
        end
    end

    @. U̇ = 𝑓.ϕ̇

    return nothing
end
