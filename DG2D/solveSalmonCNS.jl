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

    for (ϕ, D) in zip(fields, BCᵈ)
        # define field differences at faces
        @. ϕ.Δu = ϕ.u[𝒢.nodes⁻] - ϕ.u[𝒢.nodes⁺]

        # apply dirichlet boundary conditions
        if D != nothing
            dirichlet!(ϕ, D)
        end
    end

    # unpack fields
    ϕᵘ = fields[1]
    ϕᵛ = fields[2]
    ϕᵖ = fields[3]

    # compute pressure fluxes
    # might need to initialize all fluxes to zero first
    @. ϕᵖ.φˣ[𝒢.nodes⁻] = c² * ϕᵘ.Δu
    @. ϕᵖ.φʸ[𝒢.nodes⁻] = c² * ϕᵛ.Δu

    # start with pressure jump for appropriate velocity fluxes
    @. ϕᵘ.φˣ[𝒢.nodes⁻] = -ϕᵖ.Δu
    @. ϕᵛ.φʸ[𝒢.nodes⁻] = -ϕᵖ.Δu

    # compute velocity fluxes for each element
    let nGL = nBP = 0
        for Ωᵏ in 𝒢.Ω
            # get number of GL points
            GLᵏ  = (nGL + 1):(nGL + Ωᵏ.nGL)
            BPᵏ  = (nBP + 1):(nBP + Ωᵏ.nBP)
            nGL += Ωᵏ.nGL
            nBP += Ωᵏ.nBP

            # get views of computation elements
            u  = view(ϕᵘ.u,  GLᵏ)
            uˣ = view(ϕᵘ.φˣ, GLᵏ)
            uʸ = view(ϕᵘ.φʸ, GLᵏ)
            Δu = view(ϕᵘ.Δu, BPᵏ)

            v  = view(ϕᵛ.u,  GLᵏ)
            vˣ = view(ϕᵛ.φˣ, GLᵏ)
            vʸ = view(ϕᵛ.φʸ, GLᵏ)
            Δv = view(ϕᵛ.Δu, BPᵏ)

            p  = view(ϕᵖ.u,  GLᵏ)

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
            @. uʸ  = -c² / p * v * u + ν * (uʸ + ∮ʸu)
            @. vˣ  = -c² / p * u * v + ν * (vˣ + ∮ˣv)
            @. vʸ +=  c² / p * u * u + ν * (vʸ + ∮ʸv)
        end
    end

    for (ϕ, N) in zip(fields, BCⁿ)
        # Form field differences at faces for x and y partial derivatives
        @. ϕ.fˣ = ϕ.φˣ[𝒢.nodes⁻] - 1//2 * (ϕ.φˣ[𝒢.nodes⁺] + ϕ.φˣ[𝒢.nodes⁻])
        @. ϕ.fʸ = ϕ.φʸ[𝒢.nodes⁻] - 1//2 * (ϕ.φʸ[𝒢.nodes⁺] + ϕ.φʸ[𝒢.nodes⁻])

        # enfore boundary conditions for flux (neumann)
        if N != nothing
            neumann!(ϕ, N)
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
            u  = view(ϕᵘ.u,  GLᵏ)
            uˣ = view(ϕᵘ.φˣ, GLᵏ)
            uʸ = view(ϕᵘ.φʸ, GLᵏ)
            ∇u = view(ϕᵘ.∇u, GLᵏ)

            v  = view(ϕᵛ.u,  GLᵏ)
            vˣ = view(ϕᵛ.φˣ, GLᵏ)
            vʸ = view(ϕᵛ.φʸ, GLᵏ)
            ∇v = view(ϕᵛ.∇u, GLᵏ)

            p  = view(ϕᵖ.u,  GLᵏ)
            pˣ = view(ϕᵖ.φˣ, GLᵏ)
            pʸ = view(ϕᵖ.φʸ, GLᵏ)
            ∇p = view(ϕᵖ.∇u, GLᵏ)

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

            for ϕ in fields
                # compute field differences at face points
                @. ϕ.fⁿ[BPᵏ] = Ωᵏ.nˣ * ϕ.fˣ[BPᵏ] + Ωᵏ.nʸ * ϕ.fʸ[BPᵏ]

                # compute surface term
                ∮ϕ = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* ϕ.fⁿ[BPᵏ])

                # compute RHS of PDE's
                @. ϕ.u̇[GLᵏ] = ϕ.∇u[GLᵏ] + ∮ϕ
            end
        end
    end
end
