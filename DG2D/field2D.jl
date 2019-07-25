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
    u::T
    u̇::T
    ∇u::T
    φˣ::T
    φʸ::T

    Δu::T
    fˣ::T
    fʸ::T
    fⁿ::T

    r::T

    function Field2D(𝒢::Grid2D)
        # set up the solution
        u  = zeros(𝒢.nGL)
        u̇  = zeros(𝒢.nGL)
        ∇u = zeros(𝒢.nGL)
        φˣ = zeros(𝒢.nGL)
        φʸ = zeros(𝒢.nGL)

        Δu = zeros(𝒢.nBP)
        fˣ = zeros(𝒢.nBP)
        fʸ = zeros(𝒢.nBP)
        fⁿ = zeros(𝒢.nBP)

        r  = zeros(𝒢.nGL)

        return new{typeof(u)}(u,u̇,∇u,φˣ,φʸ, Δu,fˣ,fʸ,fⁿ, r)
    end
end
