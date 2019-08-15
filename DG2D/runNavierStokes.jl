include("grid2D.jl")
include("solveSalmonCNS.jl")

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
𝒢 = Grid2D(ℳ, N, periodic=true)
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
println("Time step is $dt")

# make field objects
u = Field2D(𝒢)
v = Field2D(𝒢)
p = Field2D(𝒢)

# initialize conditions
@. u.ϕ = 1.0
@. v.ϕ = 0.0
@. p.ϕ = 1.0

# parameters
stoptime = 2.
ν  = 1.0e-1
c² = 1.0

# solve equations
fields = (u, v, p)
params = (𝒢, ν, c²)
rhs!   = solveSalmonCNS!
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

solutions = rk_solver!(rhs!, fields, params, dt, Nsteps)

gr()
theme(:default)
step = floor(Int, Nsteps / 50)
step = 1

fieldNames = [ "u", "v", "p"]

@animate for t in 1:step:Nsteps
    plots = []
    for (i, sol) in enumerate(solutions)
        ploti = surface(x[:],y[:],sol[t][:], title = fieldNames[i], camera = (30,45))
        push!(plots, ploti)
    end
    display(plot(plots...))
end
