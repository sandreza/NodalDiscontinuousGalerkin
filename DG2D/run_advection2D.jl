include("grid2D.jl")
include("dg_advection2D.jl")

using Plots

# make mesh
K = 2^3
L = 2^3
xmin = ymin = -1.0
xmax = ymax = 1.0
ℳ = rectmesh2D(xmin, xmax, ymin, ymax, K, L)

filename = "Maxwell1.neu"
filepath = "./DG2D/grids/"
filename = filepath * filename
# ℳ = meshreader_gambit2D(filename)

# set number of DG elements and poly order
N = 2

# make grid
𝒢 = Grid2D(ℳ, N)
x = 𝒢.x[:,1]
y = 𝒢.x[:,2]
# plotgrid2D(𝒢)

dof = 𝒢.nGL
println("The degrees of freedom are $dof")

# determine timestep
vmax = 10 # no material here
Δx  = 𝒢.x[2,2] - 𝒢.x[1,1]
CFL = 0.75
dt  = CFL * Δx / vmax

# make field objects
u = Field2D(𝒢)

# initialize conditions
σ = 10.0
u⁰(x, y, σ) = exp(-σ * x^2 - σ * (y)^2) * cos(π/2 * x) * cos(π/2 * y)
@. u.u = [u⁰(x[i],y[i],σ) for i in 1:𝒢.nGL]
# @. u.u = 50

# parameters
α = 0 # determine upwind or central flux
vˣ = zeros(𝒢.nGL)
vʸ = zeros(𝒢.nGL)
@. vˣ = 0
@. vʸ = 0.5

# solve equations
stoptime = 10
Nsteps = ceil(Int, stoptime / dt)
# Nsteps = 3
fields = [u]
params = (𝒢, α, vˣ, vʸ)

solutions = rk_solver!(dg_advection2D!, fields, params, dt, Nsteps)

step = floor(Int, Nsteps / 50) + 1
times = 1:step:Nsteps+1
plotfield2D(times, solutions, x, y)