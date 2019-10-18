include("grid2D.jl")
include("flux2D.jl")
include("solveAdvection2D.jl")

using Plots

# make mesh
K = 10
L = 10

Lˣ = 1e6
H  = 400

xmin = 0
xmax = Lˣ
zmin = -H
zmax = 0
ℳ = rectmesh2D(xmin, xmax, zmin, zmax, K, L)

filename = "Maxwell05.neu"
filepath = "./DG2D/grids/"
filename = filepath * filename
# ℳ = meshreader_gambit2D(filename)

# set number of DG elements and poly order
N = 4

# make grid
𝒢 = Grid2D(ℳ, N, periodic=false)
x̃ = 𝒢.x[:,1]
z̃ = 𝒢.x[:,2]
plotgrid2D(𝒢)

dof = 𝒢.nGL
println("The degrees of freedom are $dof")

# determine timestep
vmax = 10 # no material here
Δx = minspacing2D(𝒢)
CFL = 0.75
dt  = CFL * Δx / vmax
dt = 60
println("Time step is $dt")

# make field objects
θ  = Field2D(𝒢)
θˣ = Field2D(𝒢)
θᶻ = Field2D(𝒢)

# initialize conditions
σ = 1.0
x⁰ = 3//4 * Lˣ
z⁰ = -H/2
θ⁰(x, z, σ) = 10 * exp(-σ * ((x - x⁰)^2 + (z - z⁰)^2))
# @. θ.ϕ = [θ⁰(x̃[i], z̃[i], σ) for i in 1:𝒢.nGL]


θ⁰(z) = 9 + 8z/H
@. θ.ϕ = [θ⁰(z̃[i]) for i in 1:𝒢.nGL]

# fluxes
φˣ = Flux2D([θˣ], [-1])
φᶻ = Flux2D([θᶻ], [-1])

# parameters
u = zeros(𝒢.nGL)
w = zeros(𝒢.nGL)

# stream function
# Ψ(x,z) = cos(π//Lˣ * (x - Lˣ//2)) * cos(π//H * (z + H/2))
ũ(x,z) = -π/Lˣ * cos(π/Lˣ * (x - Lˣ/2)) * sin(π/H * (z + H/2))
w̃(x,z) =  π/H  * sin(π/Lˣ * (x - Lˣ/2)) * cos(π/H * (z + H/2))
@. u = [ũ(x̃[i],z̃[i]) for i in 1:𝒢.nGL]
@. w = [w̃(x̃[i],z̃[i]) for i in 1:𝒢.nGL]

# solve equations
stoptime = 86400.
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

fields = [θ]
auxils = [θˣ, θᶻ]
fluxes = [φˣ, φᶻ]
params = (𝒢, u, w)

forward = rk_solver!(solveAdvection2D!, fields, fluxes, params, dt, Nsteps; auxils = auxils)

@. u = -u
@. w = -w

backward = rk_solver!(solveAdvection2D!, fields, fluxes, params, dt, Nsteps; auxils = auxils)

solutions = [forward[1]; backward[1]]

Nsteps = floor(Int, length(solutions))
step = maximum([floor(Int, Nsteps / 60), 1])
times = 1:step:Nsteps
# times = 1:100

plotfield2D(times, [solutions], x̃, z̃)
wrong = rel_error(solutions[1], solutions[end])
println("The relative error of the solution is $wrong")
