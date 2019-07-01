include("dg1D.jl")
include("dg_poisson.jl")
include("dg_heat.jl")
include("dg_advection.jl")

using Plots
using BenchmarkTools
using DifferentialEquations
using BandedMatrices

# choose eqn type
periodic = false #need to keep as false
timings = true   #to see how different linear solvers perform

# set number of DG elements and polynomial order
K = 2^4 #number of elements
n = 2^3 - 1 #polynomial order,

# for 64 total dof, K = 2^3, n = 2^3 -1 is the break even point b/w sparse and full
# for K = 2^4, n = 2^2 - 1 sparse does better
# for K = 2^2, n = 2^4 - 1 full does better

println("The degrees of freedom are ")
println((n+1) * K)

# set domain parameters
L    = 2π
xmin = 0.0
xmax = L

# generate mesh variables
𝒢 = mesh(K, n, xmin, xmax)

# generate internal variables
ι = dg(𝒢)

# set external parameters
ϰ = 1.0   #
α = 1.0   # parameter for solution, 1.0 is the example in the book
τ = 1.0
ε = external_params(ϰ, α)

# easy access
x  = 𝒢.x
u  = ι.u
uʰ = ι.uʰ
q = copy(u)
dq = copy(ι.flux)

if periodic
    make_periodic1D!(𝒢.vmapP, ι.u)
end
f = 𝒢.M * sin.(α .* x) .* α^2
@. f *= 1 / 𝒢.rx
sol = -sin.(α * x)

params = (𝒢, ι, ε, periodic, q, dq, τ)

∇² = poisson_setup(𝒢, periodic, τ)

∇² = Symmetric(∇²)
display(∇²)



s∇²  = sparse(∇²)
bands = sum(∇²[:,1] .!= 0.)-1
b∇² = BandedMatrix(zeros((n+1) * K,(n+1) * K), (bands,bands))
@. b∇² = ∇²
tmp = f[:]
comp_sol = ∇² \ tmp
@. f[:] = comp_sol
wrongness = norm(sol - f) ./ norm(sol)
println("The relative error is ")
println(wrongness)
eig_val, eig_vec =  eigen(∇²)
println("The first 10 eigenvalues are ")
println(sort(eig_val,rev=true)[1:10])
p1 = spy(s∇²)
p2 = plot(x, f, legend = false)
display(plot(p1,p2))
println("The sparsity is # nonzero / # entries")
sparsity = length(s∇².rowval) / length(s∇²)
println(sparsity)
#check to see how long it takes to solve the system
if timings == true
    println("Full solve")
    @btime comp_sol = ∇² \ tmp
    println("sparse solve")
    @btime comp_sol = s∇² \ tmp
    println("banded solve")
    @btime comp_sol = b∇² \ tmp
end
