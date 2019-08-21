include("grid2D.jl")
include("solveBurgers1D.jl")

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
N = 2^2

# make grid
𝒢 = Grid2D(ℳ, N, periodic=true)
x̃ = 𝒢.x[:,1]
ỹ = 𝒢.x[:,2]
# plotgrid2D(𝒢)

dof = 𝒢.nGL
println("The degrees of freedom are $dof")

# determine timestep
vmax = 10 # no material here
Δx = minspacing2D(𝒢)
CFL = 0.75
dt  = CFL * Δx / vmax
println("Time step is $dt")

# make field objects
u  = Field2D(𝒢)
u² = Field2D(𝒢)
uˣ = Field2D(𝒢)
uʸ = Field2D(𝒢)

# initialize conditions
ε = 0.1;
u⁰(x) = -tanh(( x + 0.5) / (2 * ε)) + 1.0
@. u.ϕ = [u⁰(x̃[i]) for i in 1:𝒢.nGL]

# solve equations
stoptime = 2.
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

# turn non linear turns on/off
α = 1

fields = [u, u², uˣ, uʸ]
params = (𝒢, ε, α)
tspan = (0.0, stoptime)

solutions = rk_solver!(solveBurgers1D!, fields, params, dt, Nsteps)
solutions = solutions[1]

Nsteps = floor(Int, length(solutions))
step = maximum([floor(Int, Nsteps / 50), 1])
times = 1:step:Nsteps
# times = 1:100
plotfield2D(times, [solutions], x̃, ỹ)
