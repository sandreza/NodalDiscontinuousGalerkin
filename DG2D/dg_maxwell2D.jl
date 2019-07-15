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
function dg_maxwell2D!(fields, params)
    # unpack params
    𝒢 = params[1] # grid parameters
    α = params[2]

    # unpack fields
    Hˣ = fields[1]
    Hʸ = fields[2]
    Eᶻ = fields[3]

    # define field differences at faces
    # need to make Δu same length as other arrays
    # each vmap is half the size of the whole array
    @. Hˣ.Δu = Hˣ.u[𝒢.vmap⁻] - Hˣ.u[𝒢.vmap⁺]
    @. Hʸ.Δu = Hʸ.u[𝒢.vmap⁻] - Hʸ.u[𝒢.vmap⁺]
    @. Eᶻ.Δu = Eᶻ.u[𝒢.vmap⁻] - Eᶻ.u[𝒢.vmap⁺]

    # impose reflective BC
    @. Hˣ.Δu[𝒢.mapᴮ] = 0
    @. Hʸ.Δu[𝒢.mapᴮ] = 0
    @. Eᶻ.Δu[𝒢.mapᴮ] = 2 * Eᶻ.u[𝒢.vmapᴮ]

    # perform calculations over elements
    let nGL = nBP = 0
        for k in 𝒢.ℳ.K
            # get element and number of GL points
            Ωᵏ = 𝒢.Ω[k]
            nGLᵏ = (nGL + 1):(nGL + length(Ωᵏ.x[:,1]))
            nBPᵏ = (nBP + 1):(nBP + Ωᵏ.nBP)
            nGL += length(Ωᵏ.x[:,1])

            # get views of computation elements
            uHˣ = view(Hˣ.u, nGLᵏ)
            uHʸ = view(Hʸ.u, nGLᵏ)
            uEᶻ = view(Eᶻ.u, nGLᵏ)

            u̇Hˣ = view(Hˣ.u̇, nGLᵏ)
            u̇Hʸ = view(Hʸ.u̇, nGLᵏ)
            u̇Eᶻ = view(Eᶻ.u̇, nGLᵏ)

            ∇Hˣ = view(Hˣ.∇u, nGLᵏ)
            ∇Hʸ = view(Hʸ.∇u, nGLᵏ)
            ∇Eᶻ = view(Eᶻ.∇u, nGLᵏ)

            ΔHˣ = Array(view(Hˣ.Δu, nBPᵏ))
            ΔHʸ = Array(view(Hʸ.Δu, nBPᵏ))
            ΔEᶻ = Array(view(Eᶻ.Δu, nBPᵏ))

            fHˣ = view(Hˣ.f, nBPᵏ)
            fHʸ = view(Hʸ.f, nBPᵏ)
            fEᶻ = view(Eᶻ.f, nBPᵏ)

            # evaluate upwind fluxes
            n̂ˣ = Ωᵏ.n̂[:,1]
            n̂ʸ = Ωᵏ.n̂[:,2]
            n̂ˣΔH = @. (n̂ˣ * ΔHˣ + n̂ʸ * ΔHʸ) * n̂ˣ
            n̂ʸΔH = @. (n̂ˣ * ΔHˣ + n̂ʸ * ΔHʸ) * n̂ʸ

            # minus isn't defined for these fluxes?????
            @. fHˣ =      n̂ʸ * ΔEᶻ + α * (n̂ˣΔH + (-1 * ΔHˣ))
            @. fHʸ = -1 * n̂ˣ * ΔEᶻ + α * (n̂ʸΔH + (-1 * ΔHʸ))
            @. fEᶻ = -1 * n̂ˣ * ΔHʸ + n̂ʸ * ΔHˣ + (-1 * α * ΔEᶻ)

            # local derivatives of the fields
            ∇Hʸ,-∇Hˣ = ∇(uEᶻ, Ωᵏ)
            ∇Eᶻ = ∇⨂(uHˣ, uHʸ, Ωᵏ)

            # compute RHS of PDE's
            u̇Hˣ += ∇Hˣ + 1//2 * Ωᵏ.lift * (Ωᵏ.volume .* fHˣ)
            u̇Hʸ += ∇Hʸ + 1//2 * Ωᵏ.lift * (Ωᵏ.volume .* fHʸ)
            u̇Eᶻ += ∇Eᶻ + 1//2 * Ωᵏ.lift * (Ωᵏ.volume .* fEᶻ)
        end
    end

    return nothing
end
