include("grid2D.jl")
include("solveMaxwell2D.jl")

using Plots
using JLD2, FileIO

# make mesh
K = 2^2
L = 2^2
xmin = ymin = -1.0
xmax = ymax = 1.0
ℳ = rectmesh2D(xmin, xmax, ymin, ymax, K, L)

name = "Maxwell1.neu"
path = "./DG2D/grids/"
name = path * name
# ℳ = meshreader_gambit2D(name)

# set number of DG elements and poly order
N = 2^3

# make grid
𝒢 = Grid2D(ℳ, N, periodic=false)
x = 𝒢.x[:,1]
y = 𝒢.x[:,2]
plotgrid2D(𝒢)

dof = 𝒢.nGL
println("The degrees of freedom are $dof")

# determine timestep
vmax = 10 # no material here
Δx = minspacing2D(𝒢)
CFL = 0.75
dt  = CFL * Δx / vmax

# make field objects
Eᶻ = Field2D(𝒢)
Hˣ = Field2D(𝒢)
Hʸ = Field2D(𝒢)

# initialize conditions
n = m = 1
@. Eᶻ.u = sin(m*π*x) * sin(n*π*y)
@. Hˣ.u = 0.0
@. Hʸ.u = 0.0

# solve equations
stoptime = 6.0
Nsteps = ceil(Int, stoptime / dt)
fields = (Hˣ, Hʸ, Eᶻ)
α = 1 # determine upwind or central flux
params = (𝒢, α)
rhs! = solveMaxwell2D!

# exact solutions
ω = π/sqrt(m^2 + n^2)
tmp = collect(1:Nsteps+1)
times = @. dt * (tmp - 1)
H̃ˣ = @. -π*n/ω * sin(m*π*x) * cos(n*π*y)
H̃ʸ = @. -π*m/ω * cos(m*π*x) * sin(n*π*y)
Ẽᶻ = @. sin(m*π*x) * sin(n*π*y)

exacts = [[], [], []]
for t in times
    tH̃ˣ = @. H̃ˣ * sin(ω*t)
    tH̃ʸ = @. H̃ʸ * sin(ω*t)
    tẼᶻ = @. Ẽᶻ * cos(ω*t)

    push!(exacts[1], tH̃ˣ)
    push!(exacts[2], tH̃ʸ)
    push!(exacts[3], tẼᶻ)
end

solutions = rk_solver!(solveMaxwell2D!, fields, params, dt, Nsteps)

gr()
theme(:default)
step = floor(Int, Nsteps / 50)

fieldNames = [ "H^{x}", "H^{y}", "E^{z}"]

@animate for t in 1:step:Nsteps
    plots = []
    for (i, sol) in enumerate(solutions)
        ploti = surface(x[:],y[:],sol[t][:], title = fieldNames[i], camera = (30,45))
        push!(plots, ploti)
    end
    display(plot(plots...))
end

@save "maxwell2D_2x2x16.jld2" solutions
