include("grid2D.jl")
include("dg_maxwell2D.jl")

using Plots

# make mesh
K = 2^2
L = 2^2
xmin = ymin = -1.0
xmax = ymax = 1.0
# ℳ = rectmesh2D(xmin, xmax, ymin, ymax, K, L)

filename = "Maxwell1.neu"
filepath = "./DG2D/grids/"
filename = filepath * filename
ℳ = meshreader_gambit2D(filename)

# set number of DG elements and poly order
N = 2

# make grid
𝒢 = Grid2D(ℳ, N)
x = 𝒢.x[:,1]
y = 𝒢.x[:,2]
plotgrid2D(𝒢)

dof = 𝒢.nGL
println("The degrees of freedom are $dof")

# determine timestep
vmax = 10 # no material here
Δx  = 𝒢.x[2,2] - 𝒢.x[1,1]
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
stoptime = 10
Nsteps = ceil(Int, stoptime / dt)
fields = (Hˣ, Hʸ, Eᶻ)
α = 0 # determine upwind or central flux
params = (𝒢, α)
rhs! = dg_maxwell2D!

# g_maxwell2D!(fields, params)
# display(Hˣ.u̇)
# display(Hʸ.u̇)
# display(Eᶻ.u̇)

solutions = rk_solver!(dg_maxwell2D!, fields, params, dt, Nsteps)

gr()
step = floor(Int, Nsteps / 50)

fieldNames = [ "H^{x}", "H^{y}", "E^{z}"]

@animate for t in 1:step:Nsteps
    plots = []
    for (i, sol) in enumerate(solutions)
        ploti = surface(x[:],y[:],sol[t][:], title = fieldNames[i], camera = (15,60))
        push!(plots, ploti)
    end
    display(plot(plots...))
end
