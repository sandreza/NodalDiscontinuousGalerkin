include("field2D.jl")
include("utils2D.jl")

"""
dg_maxwell!(u̇, u, params)

# Description

    numerical solution to 1D maxwell's equation

# Arguments

-   `u̇ = (Eʰ, Hʰ)`: container for numerical solutions to fields
-   `u  = (E , H )`: container for starting field values
-   `params = (𝒢, E, H, ext)`: mesh, E sol, H sol, and material parameters

"""
function dg_advection2D!(U̇, U, params, t)
    # unpack params
    𝒢 = params[1] # grid parameters
    α = params[2]
    h = params[end]

    @. h.u = U
    @. h.u̇ = U̇

    # unpack fields
    # h = fields[1]

    # define field differences at faces
    # println(h.u[𝒢.nodes⁻])
    # println(h.u[𝒢.nodes⁺])
    @. h.Δu = h.u[𝒢.nodes⁻] - h.u[𝒢.nodes⁺]

    # impose BC
    # @. h.u[𝒢.nodesᴮ] = 0.0

    # perform calculations over elements
    let nGL = nBP = 0
        for k in 1:𝒢.ℳ.K
            # get element and number of GL points
            Ωᵏ = 𝒢.Ω[k]
            nGLᵏ = (nGL + 1):(nGL + Ωᵏ.nGL)
            nBPᵏ = (nBP + 1):(nBP + Ωᵏ.nBP)
            nGL += Ωᵏ.nGL
            nBP += Ωᵏ.nBP

            # println("\n$k")

            # get views of params
            vˣ = view(params[3], nGLᵏ)
            vʸ = view(params[4], nGLᵏ)

            # get views of computation elements
            u = view(h.u, nGLᵏ)
            u̇ = view(h.u̇, nGLᵏ)
            ∇u = view(h.∇u, nGLᵏ)
            Δu = view(h.Δu, nBPᵏ)
            f = view(h.f, nBPᵏ)

            # println(u)
            # println(Δu)

            # evaluate fluxes
            mask = Ωᵏ.fmask
            n̂ˣ = Ωᵏ.n̂[:,1]
            n̂ʸ = Ωᵏ.n̂[:,2]
            @. f = 1//2 * (vˣ[mask][:] * Δu * n̂ˣ + vʸ[mask][:] * Δu * n̂ʸ)
            # println(f)

            # println(∇u)

            # local derivatives of the fields
            ∇⨀!(∇u, vˣ .* u, vʸ .* u, Ωᵏ)

            # println(∇u)

            # compute RHS of PDE's
            lift = 1//2 * Ωᵏ.lift * (Ωᵏ.volume .* f)
            @. u̇ = -∇u + lift
        end
    end

    @. U = h.u
    @. U̇ = h.u̇

    return nothing
end
