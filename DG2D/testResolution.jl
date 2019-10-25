include("grid2D.jl")
include("flux2D.jl")
include("solveAdvection2D.jl")

using Plots

# make mesh
K = 65

L = 1e6
H = 400
τ = 86400

xmin = 0
xmax = L
zmin = -H
zmax = 0
ℳ = rectmesh2D(xmin, xmax, zmin, zmax, K, K)

# set number of DG elements and poly order
N = 1

# make grid
𝒢 = Grid2D(ℳ, N, periodic=false)
x̃ = 𝒢.x[:,1]
ỹ = 𝒢.x[:,2]
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
x⁰ = 3//4 * L
y⁰ = -H/2
θ⁰(x, y, σ) = 10 * exp(-σ * ((x - x⁰)^2 + (y - y⁰)^2))
# @. θ.ϕ = [θ⁰(x̃[i], ỹ[i], σ) for i in 1:𝒢.nGL]


θ⁰(y) = 9 + 8y/H
@. θ.ϕ = [θ⁰(ỹ[i]) for i in 1:𝒢.nGL]

# fluxes
φˣ = Flux2D([θˣ], [-1])
φᶻ = Flux2D([θᶻ], [-1])

# parameters
u = zeros(𝒢.nGL)
v = zeros(𝒢.nGL)

# stream function
# Ψ(x,y) = L*H/τ * cos(π * (x/L - 1/2)) * cos(π * (y/H + 1/2))
ũ(x,y) = -π*L/τ * cos(π * (x/L - 1/2)) * sin(π * (y/H + 1/2))
ṽ(x,y) =  π*H/τ * sin(π * (x/L - 1/2)) * cos(π * (y/H + 1/2))
@. u = [ũ(x̃[i],ỹ[i]) for i in 1:𝒢.nGL]
@. v = [ṽ(x̃[i],ỹ[i]) for i in 1:𝒢.nGL]

# solve equations
stoptime = 86400.
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

fields = [θ]
auxils = [θˣ, θᶻ]
fluxes = [φˣ, φᶻ]
params = (𝒢, u, v)

forward = rk_solver!(solveAdvection2D!, fields, fluxes, params, dt, Nsteps; auxils = auxils)

solutions = forward[1]

@. u = -u
@. v = -v

# backward = rk_solver!(solveAdvection2D!, fields, fluxes, params, dt, Nsteps; auxils = auxils)

# solutions = [forward[1]; backward[1]]

Nsteps = floor(Int, length(solutions))
step = maximum([floor(Int, Nsteps / 60), 1])
times = 1:step:Nsteps
# times = 1:100

plotfield2D(times, [solutions], x̃, ỹ)
# wrong = rel_error(solutions[1], solutions[end])
# println("The relative error of the solution is $wrong")

max = maximum(solutions[end])
println("The max temperature is $max")
