include("grid2D.jl")
include("dg_advection2D.jl")
include("dg_helmholtz2D.jl")
include("../src/CuthillMckee.jl")

using BandedMatrices
using LinearAlgebra
using Plots

# make mesh
K = 2
L = 2
xmin = ymin = -1.0
xmax = ymax = 1.0
ℳ = rectmesh2D(xmin, xmax, ymin, ymax, K, L)

filename = "Maxwell05.neu"
filepath = "./DG2D/grids/"
filename = filepath * filename
# ℳ = meshreader_gambit2D(filename)

# set number of DG elements and poly order
N = 2

# make grid
𝒢 = Grid2D(ℳ, N, periodic=true)
x̃ = 𝒢.x[:,1]
ỹ = 𝒢.x[:,2]
dof = 𝒢.nGL
println("The degrees of freedom are $dof")
# plotgrid2D(𝒢)

# make field objects
ϕ = Field2D(𝒢)

# Boundary conditions
BCᵈ = DirichletBC(𝒢.nodesᴮ, 𝒢.mapᴮ, 0.0)
# BCᵈ = nothing
# BCⁿ = NeumannBC2D(𝒢.nodesᴮ, 𝒢.mapᴮ, 0.0, 0.0)
BCⁿ = nothing

#compute tau and define γ
γ = 10.0
τ = 1
params = [τ, γ]

# for the first helmholtz equation
# may take a while for larger matrices
∇², b = helmholtz_setup(ϕ, 𝒢, params, BCᵈ = BCᵈ, BCⁿ = BCⁿ)

# make sure its numericall symmetric
symmetric_check = sum(abs.(∇² .- (∇² + ∇²')./2)) / length(∇²) / maximum(abs.(∇²))
if symmetric_check > eps(1.0)
    println("warning the matrix is not numerically symmetric")
    ∇² = (∇² + ∇²')/2
else
    ∇² = (∇² + ∇²')/2
end
