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
function constructHelmholtzOperator(ϕ::Field2D, 𝒢::Grid2D, params; BCᵈ::Union{DirichletBC, Nothing} = nothing, BCⁿ::Union{NeumannBC2D, Nothing} = nothing)
    L = spzeros(𝒢.nGL, 𝒢.nGL)

    @. ϕ.u = 0.0
    Δq = copy(ϕ.u)
    q  = copy(ϕ.u)
    b  = copy(ϕ.u)

    # affine part of operator
    solveHelmholtz!(b, q, ϕ, 𝒢, params, BCᵈ = BCᵈ, BCⁿ = BCⁿ)
    @. q = 0.0

    for i in 1:𝒢.nGL
        q[i] = 1.0
        solveHelmholtz!(Δq, q, ϕ, 𝒢, params, BCᵈ = BCᵈ, BCⁿ = BCⁿ)
        @. L[:,i] = Δq[:] - b[:]
        q[i] = 0.0
    end

    dropϵzeros!(L)

    return L, b
end

"""
solveHelmholtz!(Δu, u, ϕ::Element2D, 𝒢::Field2D, params, BCᵈ::DirichletBC, BCⁿ::NeumannBC2D)


# Description

- Evaluate the helmholtz operator

# Arguments

- `Δu` : helholtz operator acting on u
- `u` :  the thing we want to take laplacian of
- `ϕ` : the field
- `𝒢` : the grid
- `params`: any penalty parameters that we would like to include
- `bc` : boundary condition tuple with indices
- `dbc` : boundary condition tuple with indices

"""
function solveHelmholtz!(ΔU, U, ϕ::Field2D, 𝒢::Grid2D, params; BCᵈ::Union{DirichletBC, Nothing} = nothing, BCⁿ::Union{NeumannBC2D, Nothing} = nothing)
    # unpack parameters
    τ = params[1]
    γ = params[2]

    @. ϕ.u = U

    # Form q-flux differences at faces
    @. ϕ.Δu = ϕ.u[𝒢.nodes⁻] - 1//2 * (ϕ.u[𝒢.nodes⁻] + ϕ.u[𝒢.nodes⁺])

    # Choose boundary condition type, dirichlet
    if BCᵈ != nothing
        dirichlet!(ϕ, BCᵈ)
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
            u  = view(ϕ.u,  GLᵏ)
            φˣ = view(ϕ.φˣ, GLᵏ)
            φʸ = view(ϕ.φʸ, GLᵏ)
            Δu = view(ϕ.Δu, BPᵏ)

            liftˣ = inv(Ωᵏ.M) * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nˣ .* Δu)
            liftʸ = inv(Ωᵏ.M) * Ωᵏ.∮ * (Ωᵏ.volume .* Ωᵏ.nʸ .* Δu)

            # lhs of the semi-discerte PDE, ∇⋅(q) = f , q  = ∇u, qˣ = ∂ˣu, qʸ = ∂ʸu
            # first get ∇q + flux terms
            ∇!(φˣ, φʸ, u, Ωᵏ)
            @. φˣ -= liftˣ
            @. φʸ -= liftʸ
        end
    end

    # Form field differences at faces for x and y partial derivatives
    @. ϕ.fˣ = ϕ.φˣ[𝒢.nodes⁻] - 1//2 * (ϕ.φˣ[𝒢.nodes⁺] + ϕ.φˣ[𝒢.nodes⁻])
    @. ϕ.fʸ = ϕ.φʸ[𝒢.nodes⁻] - 1//2 * (ϕ.φʸ[𝒢.nodes⁺] + ϕ.φʸ[𝒢.nodes⁻])

    # enfore boundary conditions for flux (neumann)
    if BCⁿ != nothing
        neumann!(ϕ, BCⁿ)
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
            u  = view(ϕ.u,  GLᵏ)
            u̇  = view(ϕ.u̇,  GLᵏ)
            ∇u = view(ϕ.∇u, GLᵏ)
            φˣ = view(ϕ.φˣ, GLᵏ)
            φʸ = view(ϕ.φʸ, GLᵏ)

            Δu = view(ϕ.Δu, BPᵏ)
            # the bug is here
            fˣ = view(ϕ.fˣ, BPᵏ)
            fʸ = view(ϕ.fʸ, BPᵏ)
            fⁿ = view(ϕ.fⁿ, BPᵏ)
            r  = view(ϕ.r,  GLᵏ)

            # modify with τ, remember fⁿ is field differences at face points
            @. fⁿ = Ωᵏ.nˣ * fˣ + Ωᵏ.nʸ * fʸ + τ * Δu

            # compute divergence of flux, volume term
            ∇⨀!(∇u, φˣ, φʸ, Ωᵏ)

            # compute surface term
            lift = inv(Ωᵏ.M) * Ωᵏ.∮ * (Ωᵏ.volume .* fⁿ)

            # combine the terms
            @. r = ∇u - lift - γ * u

            # multiply by J * M for making it symmetric
            mul!(u̇, Ωᵏ.M, r)
            @. u̇ *= Ωᵏ.J
        end
    end

    @. ΔU = ϕ.u̇

    return nothing
end
