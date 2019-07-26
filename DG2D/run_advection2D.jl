include("grid2D.jl")
include("dg_advection2D.jl")

using Plots
using OrdinaryDiffEq
using ReverseDiff

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
# plotgrid2D(𝒢)

# display(𝒢.Ω[1].rˣ[1, :, :])
# println(𝒢.Ω[1].volume)
# display(𝒢.Ω[1].n̂)
# display(𝒢.Ω[1].lift)

dof = 𝒢.nGL
println("The degrees of freedom are $dof")

# determine timestep
vmax = 10 # no material here
Δx = minspacing2D(𝒢)
CFL = 0.75
dt  = CFL * Δx / vmax
println("Time step is $dt")

# make field objects
u = Field2D(𝒢)

# initialize conditions
σ = 1000.0
x⁰ = 0.0
y⁰ = 0.0
u⁰(x, y, σ) = 10 * exp(-σ * ((x - x⁰)^2 + (y - y⁰)^2)) # * cos(π/2 * x) * cos(π/2 * y)
# u⁰(x, y) = 10*(y-y⁰)^2 # 10*(x-x⁰)^2
# ∇u(x, y) = 20*(x-x⁰)   # - 20*(y-y⁰)
@. u.u = [u⁰(x̃[i], ỹ[i], σ) for i in 1:𝒢.nGL]

# parameters
α  = 1. # determine upwind or central flux
vˣ = zeros(𝒢.nGL)
vʸ = zeros(𝒢.nGL)
@. vˣ = 1.
@. vʸ = 0.

# solve equations
stoptime = 2.
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

fields = [u]
params = (𝒢, α, vˣ, vʸ, u)
tspan = (0.0, stoptime)

"""
dg_advection2D!(u.u̇, u.u, params, dt)
Δ∇u = [∇u(x̃[i], ỹ[i]) - u.∇u[i] for i in 1:𝒢.nGL]
println(Δ∇u)
plot(surface(x̃[:], ỹ[:], Δ∇u, zlims = (0.0, 1.0), camera = (0, 90)))
"""

# solutions = rk_solver!(dg_advection2D!, fields, params, dt, Nsteps)
problem = ODEProblem(dg_advection2D!, u.u, tspan, params);
solutions = solve(problem, Euler(), dt=dt, adaptive = false); # AB3(), RK4(), Tsit5()

Nsteps = floor(Int, length(solutions.u))
step = maximum([floor(Int, Nsteps / 50), 1])
times = 1:step:Nsteps
# times = 1:100
plotfield2D(times, [solutions.u], x̃, ỹ)
