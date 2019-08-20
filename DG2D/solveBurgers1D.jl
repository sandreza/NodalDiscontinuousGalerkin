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
    𝑓ᵛ = fields[2]
    𝑓ᵖ = fields[3]

    # calculate q
    for Ωᵏ in 𝒢.Ω
        # get view of volume elements
        u  = view(𝑓ᵘ.ϕ,  Ωᵏ.iⱽ)
        u̇  = view(𝑓ᵘ.ϕ̇,  Ωᵏ.iⱽ)
        ∇u = view(𝑓ᵘ.∇ϕ, Ωᵏ.iⱽ)
        uˣ = view(𝑓ᵘ.φˣ, Ωᵏ.iⱽ)
        uʸ = view(𝑓ᵘ.φʸ, Ωᵏ.iⱽ)

        v  = view(𝑓ᵛ.ϕ,  Ωᵏ.iⱽ)
        q  = view(𝑓ᵖ.ϕ,  Ωᵏ.iⱽ)

        # compute volume contribution to q
        ∇!(uˣ, uʸ, u, Ωᵏ)
        @. q = sqrt(ε) * uˣ

        # compute u²
        @. v = u^2

        # compute surface contributions to q
        for f in Ωᵏ.faces
            # get views of surface elements
            u⁻ = view(𝑓ᵘ.ϕ , f.i⁻)
            u⁺ = view(𝑓ᵘ.ϕ , f.i⁺)
            Δu = view(𝑓ᵘ.Δϕ, f.i⁻)

            # define field differences at faces
            @. Δu = 1//2 * (u⁻ - u⁺)

            # impose Dirichlet BC on u
            if f.isBoundary[1]
                @. Δu = u⁰(𝒢.x[1]) - u⁻
            end

            # compute surface terms
            ∮ˣu = Ωᵏ.M⁺ * f.∮ * (f.C .* f.nˣ .* Δu)
            # combine them
            @. q -= ∮ˣu
        end

        # define physical flux
        @. ∇u = 1//2 * α * u^2 - sqrt(ε) * q

        # compute volume contributions to tendency
        ∇!(uˣ, uʸ, ∇u, Ωᵏ)
        @. u̇ = -uˣ

        # compute surface contributions to tendency
        for f in Ωᵏ.faces
            # get views of surface elements
            u⁻ = view(𝑓ᵘ.ϕ , f.i⁻)
            Δu = view(𝑓ᵘ.Δϕ, f.i⁻)
            fⁿ = view(𝑓ᵘ.fⁿ, f.i⁻)

            v⁻ = view(𝑓ᵛ.ϕ , f.i⁻)
            v⁺ = view(𝑓ᵛ.ϕ , f.i⁺)
            Δv = view(𝑓ᵛ.Δϕ, f.i⁻)

            q⁻ = view(𝑓ᵖ.ϕ , f.i⁻)
            q⁺ = view(𝑓ᵖ.ϕ , f.i⁺)
            Δq = view(𝑓ᵖ.Δϕ, f.i⁻)

            # define field differences at faces
            @. Δq = 1//2 * (q⁻ - q⁺)
            @. Δv = 1//2 * (v⁻ - v⁺)

            # impose BC on q and u²
            if f.isBoundary[1]
                @. Δq  = 0.0
                @. Δv = u⁰(𝒢.x[1])^2 - v⁻
            end

            # evaluate numerical flux
            maxu = maximum(abs.(u⁻))
            @. fⁿ = f.nˣ * (α * Δv/2 - sqrt(ε) * Δq) - 1//2 * maxu * Δu

            # compute surface term
            ∮ᶠu = Ωᵏ.M⁺ * f.∮ * (f.C .* fⁿ)

            # combine terms
            @. u̇ += ∮ᶠu
        end
    end

    return nothing
end
