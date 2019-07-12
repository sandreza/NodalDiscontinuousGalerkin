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
    Δu::T
    ∇u::T
    f::T
    r::T

    function Field2D(𝒢::Grid2D)
        # set up the solution
        u  = zeros(𝒢.nGL)
        u̇  = zeros(𝒢.nGL)
        Δu = zeros(𝒢.nGL)
        ∇u = zeros(𝒢.nGL)
        f  = zeros(𝒢.nGL)
        r  = zeros(𝒢.nGL)

        return new{typeof(u)}(u, u̇, Δu, ∇u, f)
    end
end
