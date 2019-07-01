include("dg1D.jl")
include("dg_maxwell.jl")

using Plots
using BenchmarkTools
using DifferentialEquations

# set number of DG elements and polynomial order
K = 2^4 # number of elements
n = 2^3-1 # polynomial order,

println("The degrees of freedom are ")
println((n+1) * K)

# set domain parameters
L    = 2
xmin = -1.0
xmax = xmin + L
𝒢 = mesh(K, n, xmin, xmax)
x = 𝒢.x

# set material parameters
ϵ = ones(n+1,1) * hcat(ones(1, Int(K/2)), 1.5 .* ones(1, Int(K/2)))
μ = ones(size(ϵ))
ext  = material_params(ϵ, μ)
v = @. sqrt( 1 / (ϵ * μ))
vmax = maximum(v)

# determine timestep
Δx  = minimum(x[2,:] - x[1,:])
CFL = 1.0
dt  = CFL * Δx / vmax

# initial conditions
E = dg(𝒢)
H = dg(𝒢)
@. E.u = sin(π*x) * (x < 0)
@. H.u = 0*x

# solve equations
tspan  = (0.0, 10)
params = (𝒢, E, H, ext)
rhs! = dg_maxwell!

u = [E.u, H.u]
u̇ = [E.u̇, H.u̇]

# dg_maxwell!(u̇, u, params, 0)

# prob = ODEProblem(rhs!, u, tspan, params);
# sol  = solve(prob, Tsit5(), dt=dt, adaptive = false);
sol = rk_solver!(dg_maxwell!, u̇, u, params, tspan, dt)

nt = length(sol)
num = 100
step = floor(Int, nt/num)
num = floor(Int, nt/step)
indices = step * collect(1:num)
pushfirst!(indices, 1)
push!(indices, nt)

for i in indices
   plt = plot(x, sol[i][1], xlims=(xmin,xmax), ylims = (-1.1,1.1), color = "yellow",  leg = false)
   plot!(     x, sol[i][2], xlims=(xmin,xmax), ylims = (-1.1,1.1), color = "blue", leg = false)
   display(plt)
   sleep(0.05)
end
