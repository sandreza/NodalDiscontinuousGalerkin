include("grid2D.jl")
include("flux2D.jl")
include("solveAdvection2D.jl")

using Plots

# make mesh
K = 3
L = 3
xmin = ymin = -1.0
xmax = ymax = 1.0
ℳ = rectmesh2D(xmin, xmax, ymin, ymax, K, L)

filename = "Maxwell05.neu"
filepath = "./DG2D/grids/"
filename = filepath * filename
# ℳ = meshreader_gambit2D(filename)

# set number of DG elements and poly order
N = 2^3

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
θˣ = Field2D(𝒢)
θʸ = Field2D(𝒢)

# initialize conditions
σ = 10.0
x⁰ = 0.0
y⁰ = 0.0
u⁰(x, y, σ) = 10 * exp(-σ * ((x - x⁰)^2 + (y - y⁰)^2)) * cos(π/2 * x) * cos(π/2 * y)
@. u.ϕ = [u⁰(x̃[i], ỹ[i], σ) for i in 1:𝒢.nGL]

# fluxes
φˣ = Flux2D([θˣ], [-1])
φʸ = Flux2D([θʸ], [-1])

# parameters
vˣ = zeros(𝒢.nGL)
vʸ = zeros(𝒢.nGL)
@. vˣ = 1.0
@. vʸ = 1.0

# solve equations
stoptime = 2.
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

fields = [u]
auxils = [θˣ, θʸ]
fluxes = [φˣ, φʸ]
params = (𝒢, vˣ, vʸ)

forward = rk_solver!(solveAdvection2D!, fields, fluxes, params, dt, Nsteps; auxils = auxils)

@. vˣ = -vˣ
@. vʸ = -vʸ

backward = rk_solver!(solveAdvection2D!, fields, fluxes, params, dt, Nsteps; auxils = auxils)

solutions = [forward[1]; backward[1]]

Nsteps = floor(Int, length(solutions))
step = maximum([floor(Int, Nsteps / 50), 1])
times = 1:step:Nsteps
# times = 1:100

plotfield2D(times, [solutions], x̃, ỹ)
wrong = rel_error(solutions[1], solutions[end])
println("The relative error of the solution is $wrong")
