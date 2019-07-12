include("grid2D.jl")
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
fields = (Hˣ, Hʸ, Eᶻ)
α = 1 # determine upwind or central flux
params = (𝒢, α)
rhs! = dg_maxwell2D!

sol = rk_solver!(dg_maxwell2D!, fields, params, stoptime, dt)
