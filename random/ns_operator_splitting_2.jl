using Plots
using BenchmarkTools
include("../DG2D/dg_navier_stokes.jl")
include("../DG2D/mesh2D.jl")
include("../DG2D/utils2D.jl")
include("../random/navier_stokes_structs.jl")
include("../DG2D/dg_poisson.jl")
include("../DG2D/dg_helmholtz.jl")
include("../DG2D/triangles.jl")



# define polynomial order
n = 3

# load grids
# FileName = "pvortex4A01.neu"
FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)
ι = ns_fields(mesh)

# construct boundary data
Nv, VX, VY, K, EtoV, bctype, bc_name = meshreader_gambit_bc_2D(filename)
mapT, vmapT, bc_label = build_bc_maps(mesh, bctype, bc_name)

# set time and time step and viscocity
t = 0
const Δt = 1e-2
const ν  = 1e-2

# set up the Helmholtz and Poisson operators

# calculate boundary condition map
bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)
bc_p, dbc_p = calculate_pearson_bc_p(mesh) #homogenous for pressure

# set up operators for u and v
τ = compute_τ(mesh)
γ = 1 / (ν * Δt)
params_vel = [τ, γ]
Hᵘ, bᵘ = helmholtz_setup_bc(field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
Hᵛ, bᵛ = helmholtz_setup_bc(field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)
Hᵘ = (Hᵘ + Hᵘ')/2
Hᵛ = (Hᵛ + Hᵛ')/2
dropϵzeros!(Hᵘ)
dropϵzeros!(Hᵛ)
chol_Hᵘ = cholesky(-Hᵘ)
chol_Hᵛ = cholesky(-Hᵛ)

# set up operators for p
params = [τ]
# set up matrix and affine component
#Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_wierd, bc_∇!, dbc_wierd)
Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_p, bc_∇!, dbc_p)
#slight regularization
Δᵖ = (Δᵖ + Δᵖ' ) ./ 2
dropϵzeros!(Δᵖ)
chol_Δᵖ = cholesky(-Δᵖ)
