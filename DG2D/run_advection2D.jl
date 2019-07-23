include("grid2D.jl")
include("dg_advection2D.jl")

using Plots
using OrdinaryDiffEq

# make mesh
K = 1
L = 1
xmin = ymin = -1.0
xmax = ymax = 1.0
ℳ = rectmesh2D(xmin, xmax, ymin, ymax, K, L)

filename = "Maxwell05.neu"
filepath = "./DG2D/grids/"
filename = filepath * filename
# ℳ = meshreader_gambit2D(filename)

# set number of DG elements and poly order
N = 16

# make grid
𝒢 = Grid2D(ℳ, N, periodic=true)
x = 𝒢.x[:,1]
y = 𝒢.x[:,2]
# plotgrid2D(𝒢)

display(𝒢.Ω[1].rˣ[1, :, :])
println(𝒢.Ω[1].volume)
display(𝒢.Ω[1].n̂)
display(𝒢.Ω[1].lift)

dof = 𝒢.nGL
println("The degrees of freedom are $dof")

# determine timestep
vmax = 10 # no material here

δx = minimum(setdiff!(union!([abs(x[i+1] - x[i]) for i in 1:length(x)-1]), [0.0]))
δy = minimum(setdiff!(union!([abs(y[i+1] - y[i]) for i in 1:length(y)-1]), [0.0]))
Δx = minimum([δx, δy])
CFL = 0.75
dt  = CFL * Δx / vmax

# make field objects
u = Field2D(𝒢)

# initialize conditions
σ = 100.0
u⁰(x, y, σ) = exp(-σ * x^2 - σ * y^2) * cos(π/2 * x) * cos(π/2 * y)
@. u.u = [u⁰(x[i],y[i],σ) for i in 1:𝒢.nGL]
# @. u.u = 50

# parameters
α = 1 # determine upwind or central flux
vˣ = zeros(𝒢.nGL)
vʸ = zeros(𝒢.nGL)
@. vˣ = 0.0
@. vʸ = 1.0

# solve equations
stoptime = 4.0
Nsteps = ceil(Int, stoptime / dt)
# Nsteps = 10
fields = [u]
params = (𝒢, α, vˣ, vʸ, u)
tspan = (0.0, stoptime)

# solutions = rk_solver!(dg_advection2D!, fields, params, dt, Nsteps)
problem = ODEProblem(dg_advection2D!, u.u, tspan, params);
solutions = solve(problem, RK4(), dt=dt, adaptive = false); # AB3(), RK4(), Tsit5()


Nsteps = length(solutions.u)
step = maximum([floor(Int, Nsteps / 50), 1])
times = 1:step:Nsteps
# times = 1:100
plotfield2D(times, [solutions.u], x, y)
