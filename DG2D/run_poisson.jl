using BandedMatrices

include("mesh2D.jl")
include("dg_advection.jl")
include("triangles.jl")
include("../DG2D/dg_poisson.jl")


# simulation parameters and grid
n = 7
params = [1.0]

FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
mesh = garbage_triangle(n, filename)
field = dg_garbage_triangle(mesh)
#dg_poisson!(Δu, u, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)


#∇² = poisson_setup(Δu, u, ι, params, mesh, bc_u!, bc, bc_φ!, dbc)
