include("grid2D.jl")
include("solveSalmonCNS.jl")
include("solveChorinNS.jl")

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

# make field objects
u = Field2D(𝒢)
v = Field2D(𝒢)

# auxiliary fields
uˣ = AuxiliaryField2D(𝒢)
uʸ = AuxiliaryField2D(𝒢)
vˣ = AuxiliaryField2D(𝒢)
vʸ = AuxiliaryField2D(𝒢)
uu = AuxiliaryField2D(𝒢)
uv = AuxiliaryField2D(𝒢)
vu = AuxiliaryField2D(𝒢)
vv = AuxiliaryField2D(𝒢)

# initialize conditions
@. u.ϕ = 1.0
@. v.ϕ = 0.0

# parameters
stoptime = 2.
ν  = 1.0e-1
c² = 1.0

# determine timestep
umax = maximum(abs.(u.ϕ))
vmax = maximum(abs.(v.ϕ))
cmax = maximum([umax,vmax])
Δx = minspacing2D(𝒢)
CFL = 0.25
dt  = CFL * minimum([Δx/cmax, Δx^2/ν])
println("Time step is $dt")

# solve equations
fields = (u, v)
auxil  = (uˣ, uʸ, vˣ, vʸ, uu, uv, vu, vv)
params = (𝒢, ν, c²)
rhs!   = solveChorinNS!
Nsteps = ceil(Int, stoptime / dt)
println("Number of steps is $Nsteps")

solutions = rk_solver!(rhs!, fields, params, dt, Nsteps; auxil = auxil)

gr()
theme(:default)
step = floor(Int, Nsteps / 50)
step = 1

fieldNames = ["u", "v"]

@animate for t in 1:step:Nsteps
    plots = []
    for (i, sol) in enumerate(solutions)
        ploti = surface(x[:],y[:],sol[t][:], title = fieldNames[i], camera = (0,90))
        push!(plots, ploti)
    end
    display(plot(plots...))
end
