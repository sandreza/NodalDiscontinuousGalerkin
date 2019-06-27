include("dg1D.jl")
include("dg_maxwell.jl")

using Plots
using BenchmarkTools
using DifferentialEquations

# set number of DG elements and polynomial order
K = 2^2 # number of elements
n = 2^2-1 # polynomial order,
# (for 2^8, 2^4 with 2^4-1 seems best)
println("The degrees of freedom are ")
println((n+1) * K)

# set domain parameters
L    = 4
xmin = -2.0
xmax = xmin + L
𝒢 = mesh(K, n, xmin, xmax)
x = 𝒢.x

# determine timestep
Δx  = minimum(x[2,:] - x[1,:])
CFL = 0.75
dt  = CFL * Δx / v
dt *= 0.5 / 1

# set material parameters
ϵ(x) = x > 0 ? 2 : 1
μ(x) = x > 0 ? 1 : 1
ext  = material_params(ϵ.(x), μ.(x))

# initial conditions
E = dg(𝒢)
H = dg(𝒢)
@. E.u = sin(π*x) * (x < 0)
@. H.u = 0

# solve equations
tspan  = (0.0, 10.0)
params = (𝒢, E, H, ext)
rhs! = dg_maxwell!

u  = [E.u , H.u ]
uʰ = [E.uʰ, H.uʰ]

dg_maxwell!( uʰ, u, params, 0)

prob = ODEProblem(rhs!, u, tspan, params);
sol  = solve(prob, Tsit5(), dt=dt, adaptive = false);
