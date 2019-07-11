include("mesh2D.jl")
include("dg_maxwell2D.jl")

using Plots

# set number of DG elements and poly order
K = 2^2
L = 2^3
N = 2^3-1
dof = (N+1) * K * L

println("The degrees of freedom are $dof")

# set domain parameters
xmin = ymin = -1.0
xmax = ymax = 1.0

# make grid
ℳ = rectmesh2D(xmin, xmax, ymin, ymax, K, L)
𝒢 = Grid2D(ℳ, N)
x = 𝒢.x[:,1]
y = 𝒢.x[:,2]

# determine timestep
vmax = 1 # no material here
Δx  = 𝒢.x[2,2] - 𝒢.x[1,1]
CFL = 1.0
dt  = CFL * Δx / vmax

# initialize conditions
n = m = 1
Eᶻ = @. sin(m*π*x) * sin(n*π*y)
Hˣ = zeros(length(x))
Hʸ = zeros(length(x))

# solve equations
tspan = (0.0, 10)
params = (𝒢, Hˣ, Hʸ, Eᶻ)
rhs! = dg_maxwell2D!

sol = rk_solver!(dg_maxwell2D!, u̇, u, params, tspan, dt)
