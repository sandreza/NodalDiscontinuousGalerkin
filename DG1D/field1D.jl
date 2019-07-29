include("mesh.jl")

using Plots
using BenchmarkTools


"""
Field1D(mesh)

# Description

    initialize Field1D struct

# Arguments

-   `mesh`: a mesh to compute on

# Return Values:

-   `u` : the field to be computed
-   `u̇`: numerical solutions for the field
-   `flux`: the numerical flux for the computation

"""
struct Field1D{T}
    u::T
    u̇::T
    flux::T
    r::T

    function Field1D(mesh)
        # set up the solution
        u    = copy(mesh.x)
        u̇    = copy(mesh.x)
        flux = zeros(mesh.nFP * mesh.nFaces, mesh.K)
        r    = copy(mesh.x)

        return new{typeof(u)}(u, u̇, flux, r)
    end
end
