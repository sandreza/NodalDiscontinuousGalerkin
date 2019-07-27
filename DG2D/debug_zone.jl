
# plot the boundary nodes
# scatter(mesh.x[mesh.nodesᴮ,1], mesh.x[mesh.nodesᴮ,2] , legend = false)

#=
for i in 1:length(mesh.nodesᴮ)
    ind1 = mesh.nodesᴮ[i]
    ind2 = mesh.mapᴮ[i]
    println("----------")
    println("the normal at point $(mesh.Ω[1].x[ind1,:])")
    println("is $(mesh.Ω[1].n̂[ind2,:])")
    println("----------")
end
=#

#=
for i in 1:length(mesh.x[:,1])
    println("----------")
    println("the point $(mesh.Ω[1].x[i,:]) is $(i)")
    println("----------")
end
=#

# set number of DG elements and poly order
N = 2
const debug = false
# make grid
𝒢 = Grid2D(ℳ, N, periodic=false)
mesh = 𝒢
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
γ = 00.0
τ = 0000.0
params = [τ, γ]

# for the first helmholtz equation
# may take a while for larger matrices
#@. 𝒢.Ω[1].ℰ = 0.0
m,n  = size(𝒢.Ω[1].ℰ )
e1 = zeros(n)
e1[1] = 1.0
e1[12] = 1.0
e3 = zeros(n)
e3[3] = 1.0
e3[3] = 1.0
∇², b = helmholtz_setup(ϕ, 𝒢, params, BCᵈ = BCᵈ, BCⁿ = BCⁿ);
interior = setdiff(collect(1:length(mesh.x[:,1])), mesh.nodesᴮ);
check = ∇²[interior, interior] - (∇²[interior, interior] + ∇²[interior, interior]') ./ 2;
println("check symmetry of interior nodes")
display(Array(∇²[interior, interior]))
println("check full")
display(Array(∇²))

#scale factors
rˣ,sˣ,rʸ,sʸ = partials(mesh.Ω[1].rˣ)
#manually constructed laplacian
md1 = mesh.Ω[1].M * ( mesh.Ω[1].D[1] * mesh.Ω[1].D[1] + mesh.Ω[1].D[2] * mesh.Ω[1].D[2] )
println("constructed by hand (only for lift = 0) ")
display(md1)
tmp = inv(mesh.Ω[1].M ) * ∇²
tmp = sparse(tmp)
dropϵzeros!(tmp)

display(rel_error(md1,∇²) )
asym =  ∇² - ∇²'
dropϵzeros!(asym)
println("The asymmetry is ")
display(Array(asym))
###
# load the 1D operator for checking

include("../DG1D/dg1D.jl")
include("../DG1D/dg_poisson.jl")
include("../DG1D/dg_heat.jl")
include("../DG1D/dg_advection.jl")

using Plots
using BenchmarkTools
using BandedMatrices

# choose eqn type
periodic = false #need to keep as false
timings = true   #to see how different linear solvers perform

# set number of DG elements and polynomial order
K = 2^0 #number of elements
n = N #polynomial order,

# for 64 total dof, K = 2^3, n = 2^3 -1 is the break even point b/w sparse and full
# for K = 2^4, n = 2^2 - 1 sparse does better
# for K = 2^2, n = 2^4 - 1 full does better

println("The degrees of freedom are ")
println((n+1) * K)

# set domain parameters
L    = 2
xmin = 0.0
xmax = L

# generate mesh variables
𝒢1 = Mesh(K, n, xmin, xmax)
mesh1d = Mesh(K, n, xmin, xmax)
# generate internal variables
ι = dg(𝒢1)

# set external parameters
ϰ = 1.0   #
α = 1.0   # parameter for solution, 1.0 is the example in the book
τ = 0.0  # penalty parameter
ε = (ϰ, α)

# easy access
x  = 𝒢1.x
u  = ι.u
u̇ = ι.u̇
q = copy(u)
dq = copy(ι.flux)


params = (𝒢1, ι, ε, periodic, q, dq, τ)

d1∇² = poisson_setup(𝒢1, periodic, τ)

# construct identity matrices
Iⁿ = Matrix(I, n+1, n+1)
Iᵐ = Matrix(I, n+1, n+1)

rel_error(Δ1D , tmp)
Δ1D = kron(Iᵐ, mesh1d.D * mesh1d.D) + kron(mesh1d.D * mesh1d.D, Iᵐ)
# px = kron(mesh1d.M, Iⁿ) *  kron(Iᵐ, mesh1d.M)  * mesh.Ω[1].D[1] * mesh.Ω[1].D[1]
###


###
# checking lift operator

𝒢.Ω[1].ℰ * e1

###
