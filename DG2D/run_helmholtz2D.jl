include("grid2D.jl")
include("dg_advection2D.jl")
include("dg_helmholtz2D.jl")
include("../src/CuthillMckee.jl")

using LinearAlgebra
using Plots

# make mesh
K = 2
L = 2
xmin = ymin = -1.0
xmax = ymax = 1.0
# ℳ = rectmesh2D(xmin, xmax, ymin, ymax, K, L)

filename = "Maxwell2.neu"
filepath = "./DG2D/grids/"
filename = filepath * filename
ℳ = meshreader_gambit2D(filename)

# set number of DG elements and poly order
N = 3

# make grid
𝒢 = Grid2D(ℳ, N, periodic=false)
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

display(Array(∇²))

# make sure its numericall symmetric
symmetric_check = sum(abs.(∇² .- (∇² + ∇²')./2)) / length(∇²) / maximum(abs.(∇²))
if symmetric_check > eps(1.0)
    println("warning the matrix is not numerically symmetric")
    ∇² = (∇² + ∇²')/2
else
    ∇² = (∇² + ∇²')/2
end

# output some matrix properties
println("The size of the matrix is $(size(∇²))")
i,j = findnz(∇²)
println("The bandwidth of the matrix is $(maximum(i-j)+1)")
println("The sparsity is $(length(nonzeros(∇²)) / length(∇²))")

# first create an exact solution
exact(x,y,α,β) = cos(π/2 * x * α) * cos(π/2 * y * β)

# then create a forcing function
forcing(x,y,α,β) = -((α * π/2)^2 + (β * π/2)^2 + γ) * cos(π/2 * x * α) * cos(π/2 * y * β)

# evaluate at grid points with given values for α and β
# odd for dirichlet, even for neumann
α = β = 1
frhs = [ forcing(x̃[i], ỹ[i], α, β) for i in 1:𝒢.nGL]
fsol = [   exact(x̃[i], ỹ[i], α, β) for i in 1:𝒢.nGL]

# adjust for J * mass matrix component
let nGL = 0
    for Ωᵏ in 𝒢.Ω
        GLᵏ  = (nGL + 1):(nGL + Ωᵏ.nGL)
        nGL += Ωᵏ.nGL

        frhs[GLᵏ] = Ωᵏ.J .* (Ωᵏ.M * frhs[GLᵏ])
    end
end

# subtract affine part
Δu = -(frhs - b)

# now to compute the solution
∇² = cholesky(-∇²)
u = ∇² \ Δu

# modify for neumann
u = u .- sum(u)/length(u) .+ sum(fsol)/length(fsol)

# check error
w2inf = rel_error(u, fsol)
println("The relative error in computing the solution is $(w2inf)")
println("----------------")
