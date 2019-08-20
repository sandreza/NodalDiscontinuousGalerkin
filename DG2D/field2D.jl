include("grid2D.jl")

abstract type AbstractField2D end

"""
Field2D(𝒢::Grid2D)

# Description

    Contains all the computational elements necessary to evolve a field in time

# Arguments

-   `mesh`: a mesh to compute on

# Return Values:

-   `u` : the field to be computed
-   `u̇`: numerical solutions for the field
-   `flux`: the numerical flux for the computation

"""
struct Field2D{T} <: AbstractField2D
    # volume terms
    ϕ::T
    ϕ̇::T
    ∇ϕ::T
    φˣ::T
    φʸ::T

    # surface terms
    ϕ⁺::T
    Δϕ::T
    fˣ::T
    fʸ::T
    fⁿ::T

    # residual
    r::T

    function Field2D(𝒢::Grid2D)
        # set up the solution
        ϕ  = zeros(𝒢.nGL)
        ϕ̇  = zeros(𝒢.nGL)
        ∇ϕ = zeros(𝒢.nGL)
        φˣ = zeros(𝒢.nGL)
        φʸ = zeros(𝒢.nGL)

        ϕ⁺ = zeros(𝒢.nGL)
        Δϕ = zeros(𝒢.nGL)
        fˣ = zeros(𝒢.nGL)
        fʸ = zeros(𝒢.nGL)
        fⁿ = zeros(𝒢.nGL)

        r  = zeros(𝒢.nGL)

    return new{typeof(ϕ)}(ϕ,ϕ̇,∇ϕ,φˣ,φʸ, ϕ⁺,Δϕ,fˣ,fʸ,fⁿ, r)
    end
end
