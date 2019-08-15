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
    ϕ::T
    ϕ̇::T
    ∇ϕ::T
    φˣ::T
    φʸ::T

    Δϕ::T
    fˣ::T
    fʸ::T
    fⁿ::T

    r::T

    function Field2D(𝒢::Grid2D)
        # set up the solution
        ϕ  = zeros(𝒢.nGL)
        ϕ̇  = zeros(𝒢.nGL)
        ∇ϕ = zeros(𝒢.nGL)
        φˣ = zeros(𝒢.nGL)
        φʸ = zeros(𝒢.nGL)

        Δϕ = zeros(𝒢.nBP)
        fˣ = zeros(𝒢.nBP)
        fʸ = zeros(𝒢.nBP)
        fⁿ = zeros(𝒢.nBP)

        r  = zeros(𝒢.nGL)

        return new{typeof(ϕ)}(ϕ,ϕ̇,∇ϕ,φˣ,φʸ, Δϕ,fˣ,fʸ,fⁿ, r)
    end
end
