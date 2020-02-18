include("field2D.jl")
include("utils2D.jl")
include("boundaryConditions2D.jl")

"""
computeTau(𝒢::Grid2D)

# Description

-   Computes the tau parameter in NDG
-   requires adding a face structure./// so not having yet

# Arguments

-   `𝒢`: Grid object

# Returns

-   `τ`: the value of τ at every grid point. (in the code could be either)
"""
function computeTau(𝒢::Grid2D)
    mat⁻ = @. 𝒢.J[𝒢.nodes⁻] / 𝒢.sJ[:]
    mat⁺ = @. 𝒢.J[𝒢.nodes⁺] / 𝒢.sJ[:]

    mat⁻ = @. Ωᵏ.volume[𝒢.nodes⁻[BPᵏ]]



    for i in 1:length(mat⁺)
        if mat⁺[i] < mat⁻[i]
            τ[i] = 2 * mat⁺[i]
        else
            τ[i] = 2 * mat⁻[i]
        end
    end

    τ = zeros(𝒢.nBP)
    let nBP = 0
        for Ωᵏ in 𝒢.Ω
            BPᵏ  = (nBP + 1):(nBP + Ωᵏ.nBP)
            nBP += Ωᵏ.nBP

            @. τ[BPᵏ] = Ωᵏ.nGL / τ[BPᵏ]
        end
    end

    return τ
end

# builds the affine operator (one column at a time) (sparse matrix)
# here Δ[u] = L[u] + b (b is where the boundary conditions go as a forcing term)
function constructHelmholtzOperator(𝑓::Field2D, 𝒢::Grid2D, params; BCᵈ::Union{DirichletBC, Nothing} = nothing, BCⁿ::Union{NeumannBC2D, Nothing} = nothing)
    L = spzeros(𝒢.nGL, 𝒢.nGL)

    @. 𝑓.ϕ = 0.0
    Δq = copy(𝑓.ϕ)
    q  = copy(𝑓.ϕ)
    b  = copy(𝑓.ϕ)

    # affine part of operator
    solveHelmholtz!(b, q, 𝑓, 𝒢, params, BCᵈ = BCᵈ, BCⁿ = BCⁿ)
    @. q = 0.0

    for i in 1:𝒢.nGL
        q[i] = 1.0
        solveHelmholtz!(Δq, q, 𝑓, 𝒢, params, BCᵈ = BCᵈ, BCⁿ = BCⁿ)
        @. L[:,i] = Δq[:] - b[:]
        q[i] = 0.0
    end

    dropϵzeros!(L)

    return L, b
end

"""
solveHelmholtz!(Δu, u, 𝑓::Element2D, 𝒢::Field2D, params, BCᵈ::DirichletBC, BCⁿ::NeumannBC2D)


# Description

- Evaluate the helmholtz operator

# Arguments

- `Δu` : helholtz operator acting on u
- `u` :  the thing we want to take laplacian of
- `𝑓` : the field
- `𝒢` : the grid
- `params`: any penalty parameters that we would like to include
- `bc` : boundary condition tuple with indices
- `dbc` : boundary condition tuple with indices

"""
function solveHelmholtz!(ΔU, U, 𝑓::Field2D, 𝒢::Grid2D, params; BCᵈ::Union{DirichletBC, Nothing} = nothing, BCⁿ::Union{NeumannBC2D, Nothing} = nothing)
    # unpack parameters
    τ = params[1]
    γ = params[2]

    @. 𝑓.ϕ = U

    # Form q-flux differences at faces
    @. 𝑓.Δϕ = 𝑓.ϕ[𝒢.nodes⁻] - 1//2 * (𝑓.ϕ[𝒢.nodes⁻] + 𝑓.ϕ[𝒢.nodes⁺])

    # Choose boundary condition type, dirichlet
    if BCᵈ != nothing
        dirichlet!(𝑓, BCᵈ)
    end

    # compute fluxes for each element
    let nGL = nBP = 0
        for Ωᵏ in 𝒢.Ω
            # get number of GL points
            GLᵏ  = (nGL + 1):(nGL + Ωᵏ.nGL)
            BPᵏ  = (nBP + 1):(nBP + Ωᵏ.nBP)
            nGL += Ωᵏ.nGL
            nBP += Ωᵏ.nBP

            # get views of computation elements
            u  = view(𝑓.ϕ,  GLᵏ)
            φˣ = view(𝑓.φˣ, GLᵏ)
            φʸ = view(𝑓.φʸ, GLᵏ)
            Δu = view(𝑓.Δϕ, BPᵏ)

            ∮ˣu = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nˣ .* Δu)
            ∮ʸu = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nʸ .* Δu)

            # lhs of the semi-discerte PDE, ∇⋅(q) = f , q  = ∇u, qˣ = ∂ˣu, qʸ = ∂ʸu
            # first get ∇q + flux terms
            ∇!(φˣ, φʸ, u, Ωᵏ)
            @. φˣ -= ∮ˣu
            @. φʸ -= ∮ʸu
        end
    end

    # Form field differences at faces for x and y partial derivatives
    @. 𝑓.fˣ = 𝑓.φˣ[𝒢.nodes⁻] - 1//2 * (𝑓.φˣ[𝒢.nodes⁺] + 𝑓.φˣ[𝒢.nodes⁻])
    @. 𝑓.fʸ = 𝑓.φʸ[𝒢.nodes⁻] - 1//2 * (𝑓.φʸ[𝒢.nodes⁺] + 𝑓.φʸ[𝒢.nodes⁻])

    # enfore boundary conditions for flux (neumann)
    if BCⁿ != nothing
        neumann!(𝑓, BCⁿ)
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
            u  = view(𝑓.ϕ,  GLᵏ)
            u̇  = view(𝑓.ϕ̇,  GLᵏ)
            ∇u = view(𝑓.∇ϕ, GLᵏ)
            φˣ = view(𝑓.φˣ, GLᵏ)
            φʸ = view(𝑓.φʸ, GLᵏ)

            Δu = view(𝑓.Δϕ, BPᵏ)
            fˣ = view(𝑓.fˣ, BPᵏ)
            fʸ = view(𝑓.fʸ, BPᵏ)
            fⁿ = view(𝑓.fⁿ, BPᵏ)
            r  = view(𝑓.r,  GLᵏ)

            # compute divergence of flux, volume term
            ∇⨀!(∇u, φˣ, φʸ, Ωᵏ)

            # modify with τ, remember fⁿ is field differences at face points
            @. fⁿ = Ωᵏ.nˣ * fˣ + Ωᵏ.nʸ * fʸ + τ * Δu

            # compute surface term
            lift = Ωᵏ.M⁺ * Ωᵏ.∮ * (Ωᵏ.volume .* fⁿ)

            # combine the terms
            @. r = ∇u - lift - γ * u

            # multiply by J * M for making it symmetric
            mul!(u̇, Ωᵏ.M, r)
            @. u̇ *= Ωᵏ.J
        end
    end

    @. ΔU = 𝑓.ϕ̇

    return nothing
end
