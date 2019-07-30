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
    h = params[end]

    @. h.u = U

    # define field differences at faces
    @. h.Δu = h.u[𝒢.nodes⁻] - h.u[𝒢.nodes⁺]

    # impose BC
    # @. h.u[𝒢.nodesᴮ] = 0.0

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
            u  = view(h.u,  GLᵏ)
            u̇  = view(h.u̇,  GLᵏ)
            ∇u = view(h.∇u, GLᵏ)
            Δu = view(h.Δu, BPᵏ)
            f  = view(h.fⁿ, BPᵏ)

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

    @. U̇ = h.u̇

    return nothing
end
