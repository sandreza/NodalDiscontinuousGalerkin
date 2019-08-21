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
    # field value and tendency
    ϕ::T
    ϕ̇::T
    
    # volume contributions to tendency
    𝚽::T

    # physical fluxes
    φˣ::T
    φʸ::T

    # numerical fluxes
    fˣ::T
    fʸ::T

    # jump in the flux
    Δf::T

    # surface contributions to the tendency
    ∮f::T

    # residual for RK4 methods
    r::T

    function Field2D(𝒢::Grid2D)

        ϕ  = zeros(𝒢.nGL)
        ϕ̇  = zeros(𝒢.nGL)
        𝚽 = zeros(𝒢.nGL)
        φˣ = zeros(𝒢.nGL)
        φʸ = zeros(𝒢.nGL)

        fˣ = zeros(𝒢.nGL)
        fʸ = zeros(𝒢.nGL)
        Δf = zeros(𝒢.nGL)
        ∮f = zeros(𝒢.nGL)

        r  = zeros(𝒢.nGL)

    return new{typeof(ϕ)}(ϕ,ϕ̇,𝚽, φˣ,φʸ, fˣ,fʸ, Δf,∮f, r)

    end
end
