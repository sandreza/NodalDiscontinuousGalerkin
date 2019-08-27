include("grid2D.jl")
include("solveBurgers1D.jl")

using Plots

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
𝒢 = Grid2D(ℳ, N, periodic=false)
x̃ = 𝒢.x[:,1]
ỹ = 𝒢.x[:,2]
# plotgrid2D(𝒢)

dof = 𝒢.nGL
println("The degrees of freedom are $dof")

# make field objects
u  = Field2D(𝒢)
u² = Field2D(𝒢)
uˣ = Field2D(𝒢)
uʸ = Field2D(𝒢)

# initialize conditions
ν = 0.1;
t⁰ = 0
u⁰(x,t) = -tanh(( x + 0.5 - t) / (2 * ν)) + 1.0
@. u.ϕ = [u⁰(x̃[i],t⁰) for i in 1:𝒢.nGL]

# determine timestep
umax = maximum(abs.(u.ϕ))
Δx = minspacing2D(𝒢)
CFL = 0.25
dt  = CFL * minimum([Δx/umax, Δx^2/ν])
println("Time step is $dt")

# solve equations
stoptime = 2.
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

# turn non linear turns on/off
α = 1

# turn 2D on/off
β = 0

# fluxes
φˣᵤ = Flux2D([u², uˣ], [-0.5 * α, sqrt(ν)])
φʸᵤ = Flux2D([u², uʸ], [-0.5 * α * β, β * sqrt(ν)])
φᵘ  = Flux2D([u], sqrt(ν))

fields = [u]
fluxes = [φᵘ, φˣᵤ, φʸᵤ]
auxils = [u², uˣ, uʸ]
params = (𝒢, ν, α, β)
tspan = (0.0, stoptime)

solutions = rk_solver!(solveBurgers1D!, fields, fluxes, params, dt, Nsteps; auxils = auxils)
solutions = solutions[1]

Nsteps = floor(Int, length(solutions))
step = maximum([floor(Int, Nsteps / 50), 1])
times = 1:step:Nsteps

exacts = []
for time in times
    t = dt * time
    uᵗ = @. [u⁰(x̃[i],t) for i in 1:𝒢.nGL]
    push!(exacts, uᵗ)
end

diffs = []
for (sol, exact) in zip(solutions, exacts)
    diff = @. sol - exact
    push!(diffs, diff)
end


# times = 1:100
plotfield2D(times, [solutions, exacts, diffs], x̃, ỹ)
wrong = rel_error(solutions[end], exacts[end])
println("The relative error of the solution is $wrong")
