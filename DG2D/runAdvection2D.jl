include("grid2D.jl")
include("solveAdvection2D.jl")

using Plots
using OrdinaryDiffEq

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
N = 2^4

# make grid
𝒢 = Grid2D(ℳ, N, periodic=true)
x̃ = 𝒢.x[:,1]
ỹ = 𝒢.x[:,2]
plotgrid2D(𝒢)

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
σ = 10.0
x⁰ = 0.0
y⁰ = 0.0
u⁰(x, y, σ) = 10 * exp(-σ * ((x - x⁰)^2 + (y - y⁰)^2)) * cos(π/2 * x) * cos(π/2 * y)
# u⁰(x, y) = 10*(y-y⁰)^2 # 10*(x-x⁰)^2
# ∇u(x, y) = 20*(x-x⁰)   # - 20*(y-y⁰)
@. u.u = [u⁰(x̃[i], ỹ[i], σ) for i in 1:𝒢.nGL]

# parameters
α  = 1. # determine upwind or central flux
vˣ = zeros(𝒢.nGL)
vʸ = zeros(𝒢.nGL)
@. vˣ = 1.0
@. vʸ = 1.0

# solve equations
stoptime = 2.
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

fields = [u]
params = (𝒢, α, vˣ, vʸ, u)
tspan = (0.0, stoptime)

# solutions = rk_solver!(solveAdvection2D!, fields, params, dt, Nsteps)
problem = ODEProblem(solveAdvection2D!, u.u, tspan, params);
forward = solve(problem, RK4(), dt=dt, adaptive = false); # AB3(), RK4(), Tsit5()

@. vˣ = -vˣ
@. vʸ = -vʸ

problem = ODEProblem(solveAdvection2D!, u.u, tspan, params);
backward = solve(problem, RK4(), dt=dt, adaptive = false); # AB3(), RK4(), Tsit5()

solutions = [forward.u; backward.u]

Nsteps = floor(Int, length(solutions))
step = maximum([floor(Int, Nsteps / 50), 1])
times = 1:step:Nsteps
# times = 1:100
plotfield2D(times, [solutions], x̃, ỹ)
wrong = rel_error(solutions[1], solutions[end])
println("The relative error of the solution is $wrong")
