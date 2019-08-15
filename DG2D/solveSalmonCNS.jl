include("field2D.jl")
include("utils2D.jl")

"""
solveSalmonCNS!(fields, params)

# Description

    numerical solution to 1D maxwell's equation

# Arguments

-   `fields = (u, v, p)`: velocity in each dimension and a pressure
-   `params = (𝒢, ν, c²)`: grid struct, viscosity, and speed of sound
-   `BCᵈ = (Dᵘ, Dᵛ, Dᵖ)`: dirichlet boundary conditions for each field
-   `BCⁿ = (Nᵘ, Nᵛ, Nᵖ)`:   neumann boundary conditions for each field

"""
function solveSalmonCNS!(fields, params; BCᵈ = [nothing, nothing, nothing], BCⁿ = [nothing, nothing, nothing])
    # unpack parameters
    𝒢  = params[1]
    ν  = params[2]
    c² = params[3]

    for (𝑓, D) in zip(fields, BCᵈ)
        # define field differences at faces
        @. 𝑓.Δϕ = 𝑓.ϕ[𝒢.nodes⁻] - 𝑓.ϕ[𝒢.nodes⁺]

        # apply dirichlet boundary conditions
        if D != nothing
            dirichlet!(𝑓, D)
        end
    end

    # unpack fields
    𝑓ᵘ = fields[1]
    𝑓ᵛ = fields[2]
    𝑓ᵖ = fields[3]

    # compute pressure fluxes
    # might need to initialize all fluxes to zero first
    @. 𝑓ᵖ.φˣ[𝒢.nodes⁻] = c² * 𝑓ᵘ.Δϕ
    @. 𝑓ᵖ.φʸ[𝒢.nodes⁻] = c² * 𝑓ᵛ.Δϕ

    # start with pressure jump for appropriate velocity fluxes
    @. 𝑓ᵘ.φˣ[𝒢.nodes⁻] = -𝑓ᵖ.Δϕ
    @. 𝑓ᵘ.φʸ[𝒢.nodes⁻] = 0
    @. 𝑓ᵛ.φˣ[𝒢.nodes⁻] = 0
    @. 𝑓ᵛ.φʸ[𝒢.nodes⁻] = -𝑓ᵖ.Δϕ

    # compute velocity fluxes for each element
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

            v  = view(𝑓ᵛ.ϕ,  GLᵏ)
            vˣ = view(𝑓ᵛ.φˣ, GLᵏ)
            vʸ = view(𝑓ᵛ.φʸ, GLᵏ)
            Δv = view(𝑓ᵛ.Δϕ, BPᵏ)

            p  = view(𝑓ᵖ.ϕ,  GLᵏ)

            # compute surface integrals
            ∮ˣu = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nˣ .* Δu)
            ∮ʸu = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nʸ .* Δu)
            ∮ˣv = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nˣ .* Δv)
            ∮ʸv = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nʸ .* Δv)

            # compute gradients for laplacian
            ∇!(uˣ, uʸ, u, Ωᵏ)
            ∇!(vˣ, vʸ, v, Ωᵏ)

            # compute velocity fluxes
            @. uˣ +=  c² / p * v * v + ν * (uˣ + ∮ˣu)
            @. uʸ += -c² / p * v * u + ν * (uʸ + ∮ʸu)
            @. vˣ += -c² / p * u * v + ν * (vˣ + ∮ˣv)
            @. vʸ +=  c² / p * u * u + ν * (vʸ + ∮ʸv)
        end
    end

    for (𝑓, N) in zip(fields, BCⁿ)
        # Form field differences at faces for x and y partial derivatives
        @. 𝑓.fˣ = 𝑓.φˣ[𝒢.nodes⁻] - 1//2 * (𝑓.φˣ[𝒢.nodes⁺] + 𝑓.φˣ[𝒢.nodes⁻])
        @. 𝑓.fʸ = 𝑓.φʸ[𝒢.nodes⁻] - 1//2 * (𝑓.φʸ[𝒢.nodes⁺] + 𝑓.φʸ[𝒢.nodes⁻])

        # enfore boundary conditions for flux (neumann)
        if N != nothing
            neumann!(𝑓, N)
        end
    end

    # compute tendecy for each element
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
            ∇u = view(𝑓ᵘ.∇ϕ, GLᵏ)

            v  = view(𝑓ᵛ.ϕ,  GLᵏ)
            vˣ = view(𝑓ᵛ.φˣ, GLᵏ)
            vʸ = view(𝑓ᵛ.φʸ, GLᵏ)
            ∇v = view(𝑓ᵛ.∇ϕ, GLᵏ)

            p  = view(𝑓ᵖ.ϕ,  GLᵏ)
            pˣ = view(𝑓ᵖ.φˣ, GLᵏ)
            pʸ = view(𝑓ᵖ.φʸ, GLᵏ)
            ∇p = view(𝑓ᵖ.∇ϕ, GLᵏ)

            # compute laplacian
            ∇⨀!(∇u, uˣ, uʸ, Ωᵏ) #### must come before gradient
            ∇⨀!(∇v, vˣ, vʸ, Ωᵏ) #### must come before gradient

            # compute partials for curl
            ∇!(uˣ, uʸ, u, Ωᵏ)   #### gradient overwrites values
            ∇!(vˣ, vʸ, v, Ωᵏ)   #### gradient overwrites values

            # compute full inner derivative
            @. ∇u =  c² / p * v * (vˣ - uʸ) - pˣ + ν * ∇u
            @. ∇v = -c² / p * u * (vˣ - uʸ) - pʸ + ν * ∇v

            # compute pressure derivative
            ∇⨀!(∇p, u, v, Ωᵏ)
            @. ∇p *= -c²

            for 𝑓 in fields
                # compute field differences at face points
                @. 𝑓.fⁿ[BPᵏ] = Ωᵏ.nˣ * 𝑓.fˣ[BPᵏ] + Ωᵏ.nʸ * 𝑓.fʸ[BPᵏ]

                # compute surface term
                ∮𝑓 = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* 𝑓.fⁿ[BPᵏ])

                # compute RHS of PDE's
                @. 𝑓.ϕ̇[GLᵏ] = 𝑓.∇ϕ[GLᵏ] + ∮𝑓
            end
        end
    end
end
