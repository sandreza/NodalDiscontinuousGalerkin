
abstract type AbstractBoundaryCondition end
abstract type AbstractDirichletBC <: AbstractBoundaryCondition end
abstract type AbstractNeumannBC2D <: AbstractBoundaryCondition end

"""
DirichletBC(nodes, map, f)

# Description

    initialize a DirichletBC struct containing the boundary nodes, the boundary map, and a function to apply at the boundary

# Arguments

-   `nodes`: the boundary nodes
-   `map`: the boundary map
-   `f`: function to apply at the boundary

# Return Values:

    return a properly initiliazed DirichletBC object

"""
struct DirichletBC{S,T} <:AbstractDirichletBC
    nodes::S
    map::S
    f::T

    function DirichletBC(nodes, map, f)
        return new{typeof(nodes),typeof(f),}(nodes, map, f)
    end
end

"""
dirichlet!(ϕ::Field2D, B::DirichletBC)

# Description

    Impose dirichlet boundary conditions on a field

# Arguments

-   `ϕ`: a field
-   `B`: a DirichletBC struct

# Return Values

"""
function dirichlet!(ϕ::Field2D, B::DirichletBC)
    @. ϕ.Δu[B.map] = ϕ.u[B.nodes] - B.f

    return nothing
end

"""
NeumannBC2D(nodes, map, fˣ, fʸ)

# Description

    initialize a DirichletBC struct containing the boundary nodes, the boundary map, and functions to apply at the boundary

# Arguments

-   `nodes`: the boundary nodes
-   `map`: the boundary map
-   `fˣ`: function to apply at the boundary to the x derivative
-   `fʸ`: function to apply at the boundary to the y derivative

# Return Values:

    return a properly initiliazed NeumannBC2D object

"""
struct NeumannBC2D{S,T} <:AbstractNeumannBC2D
    nodes::S
    map::S
    fˣ::T
    fʸ::T

    function NeumannBC2D(nodes, map, fˣ, fʸ)
        return new{typeof(nodes),typeof(fˣ),}(nodes, map, fˣ, fʸ)
    end
end

"""
neumann!(ϕ::Field2D, B:NeumannBC2D)

# Description

    Impose neumann boundary conditions on a field

# Arguments

-   `ϕ`: a field
-   `B`: a NeumannBC2D struct

# Return Values

"""
function neumann!(ϕ::Field2D, B::NeumannBC2D)
    @. ϕ.fˣ[B.map] = ϕ.φˣ[B.nodes] - B.fˣ
    @. ϕ.fʸ[B.map] = ϕ.φʸ[B.nodes] - B.fʸ

    return nothing
end
