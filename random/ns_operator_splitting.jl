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
n = 10

# load grids
FileName = "pvortex4A01.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)
ι = ns_fields(mesh)

# construct boundary data
Nv, VX, VY, K, EtoV, bctype, bc_name = meshreader_gambit_bc_2D(filename)
mapT, vmapT, bc_label = build_bc_maps(mesh, bctype, bc_name)

# construct exact solution
t = 0
ν = 1e-2 # should declare as constant
u_exact = similar(mesh.x)
v_exact = similar(mesh.x)
p_exact = similar(mesh.x)
pearson_vortex!(u_exact, v_exact, p_exact, mesh, ν, t)

# construct boundary conditions for exact solution
ubc = zeros(size(mesh.nx))
vbc = zeros(size(mesh.nx))
pbc = zeros(size(mesh.nx))
undtbc = zeros(size(mesh.nx))



println("------------")


# check pointwise convergece of pressure equation
# use exact boundary conditions as the check
# will also check the divergence and everything else

# ∂ᵗ u = -u⋅∇u + ν Δu - ∇p
# ∇⋅u = 0
# from whence
# Δp = ∇ ⋅ ( -u⋅∇u ), we are checking this equation


# compute operators for helmholtz
t = 0
Δt = 0.1

u_exact = eval_grid(u_analytic, mesh, t)
v_exact = eval_grid(v_analytic, mesh, t)

# ∂ˣ
∂ˣu_exact = eval_grid(∂ˣu_analytic, mesh, t)
∂ˣv_exact = eval_grid(∂ˣv_analytic, mesh, t)
# ∂ʸ
∂ʸu_exact = eval_grid(∂ʸu_analytic, mesh, t)
∂ʸv_exact = eval_grid(∂ʸv_analytic, mesh, t)

dirichlet_u_bc = u_exact[mesh.vmapB];
neumann_u_bcx = ∂ˣu_exact[mesh.vmapB];
neumann_u_bcy = ∂ʸu_exact[mesh.vmapB];
neumann_v_bcx = ∂ˣv_exact[mesh.vmapB];
neumann_v_bcy = ∂ʸv_exact[mesh.vmapB];

bc_u = (mesh.vmapB, mesh.mapB, dirichlet_u_bc)
dbc_u = ([],[],0.0,0.0)
#dbc_u = (mesh.vmapB, mesh.mapB, neumann_u_bcx, neumann_u_bcy ) # for testing
dirichlet_v_bc = v_exact[mesh.vmapB];
bc_v = (mesh.vmapB, mesh.mapB, dirichlet_v_bc)
dbc_v = ([],[],0.0,0.0)
#dbc_v = (mesh.vmapB, mesh.mapB, neumann_v_bcx, neumann_v_bcy) #for testing

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

# preliminary orszag thing
u⁰ = eval_grid(u_analytic, mesh, 0)
v⁰ = eval_grid(v_analytic, mesh, 0)
@. ι.u.ϕ = u⁰
@. ι.v.ϕ = v⁰

# compute boundary conditions for u and v
t += Δt

u_exact = eval_grid(u_analytic, mesh, t)
v_exact = eval_grid(v_analytic, mesh, t)

dirichlet_u_bc = u_exact[mesh.vmapB];
bc_u = (mesh.vmapB, mesh.mapB, dirichlet_u_bc)
dbc_u = ([],[],0.0,0.0)
dirichlet_v_bc = v_exact[mesh.vmapB];
bc_v = (mesh.vmapB, mesh.mapB, dirichlet_v_bc)
dbc_v = ([],[],0.0,0.0)

# get affine part of operator
zero_value = 0.0 * u_exact

dg_helmholtz_bc!(bᵘ, zero_value, field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
dg_helmholtz_bc!(bᵛ, zero_value, field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)

# now set up pressure
# location of boundary grid points for dirichlet bc
bc = ([],[],0.0)
dbc = (mesh.vmapB, mesh.mapB, 0.0, 0.0)

bc = ([mesh.vmapB[1]], [mesh.mapB[1]], 0.0)
dbc = (mesh.vmapB[2:end], mesh.mapB[2:end], 0.0, 0.0)


# set up τ matrix
τ = compute_τ(mesh)
params = [τ]

# set up matrix and affine component
#Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_wierd, bc_∇!, dbc_wierd)
Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc, bc_∇!, dbc)
#slight regularization
Δᵖ = (Δᵖ + Δᵖ' ) ./ 2
dropϵzeros!(Δᵖ)
chol_Δᵖ = cholesky(-Δᵖ)
#Δᵖ  = lu(-Δᵖ)


# this is where the method actually starts
# prior to this is set up

# step 1, compute the advection step
# first u
sym_advec!(ι.u.φⁿ, u⁰, v⁰, u⁰, mesh)
@. ι.u.φⁿ *= -Δt
@. ι.u.φⁿ += u⁰
# then v
sym_advec!(ι.v.φⁿ, u⁰, v⁰, v⁰, mesh)
@. ι.v.φⁿ *= -Δt
@. ι.v.φⁿ += v⁰

# step two, project
#first compute boundary conditions for pressure
fu¹ = fv¹ = 0.0

∂pˣ, ∂pʸ = compute_pressure_terms(u⁰, v⁰, ν, fu¹, fv¹, t-Δt, mesh)

#@. ∂pˣ *= 0.0
#@. ∂pʸ *= 0.0
#=
bc = ([],[],0.0)
dbc = (mesh.vmapB, mesh.mapB, ∂pˣ[mesh.vmapB], ∂pʸ[mesh.vmapB])
=#
bc_p = ([mesh.vmapB[1]], [mesh.mapB[1]], 0.0)
dbc_p = (mesh.vmapB[2:end], mesh.mapB[2:end], ∂pˣ[mesh.vmapB[2:end]], ∂pʸ[mesh.vmapB[2:end]])
dg_poisson_bc!(bᵖ, zero_value, field, params_vel, mesh, bc!, bc_p, bc_∇!, dbc_p)

# take the divergence of the solution
rhsᵖ = similar(ι.p.ϕ)
∇⨀!(rhsᵖ , ι.u.φⁿ, ι.v.φⁿ, mesh)

# construct the right hand side for poissons equation
frhsᵖ = mesh.J .* (mesh.M * rhsᵖ) - bᵖ
@. frhsᵖ *= -1.0
# solve the linear system
p = reshape(chol_Δᵖ \ frhsᵖ[:], size(mesh.x));
# compute the gradient
∇!(ι.p.∂ˣ,ι.p.∂ʸ, p, mesh)
# project
ũ = ι.u.φⁿ - ι.p.∂ˣ
ṽ = ι.v.φⁿ - ι.p.∂ʸ

# now for the diffusive step

# get bc
u_exact = eval_grid(u_analytic, mesh, t)
v_exact = eval_grid(v_analytic, mesh, t)
dirichlet_u_bc = u_exact[mesh.vmapB];
bc_u = (mesh.vmapB, mesh.mapB, dirichlet_u_bc)
dbc_u = ([],[],0.0,0.0)
dirichlet_v_bc = v_exact[mesh.vmapB];
bc_v = (mesh.vmapB, mesh.mapB, dirichlet_v_bc)
dbc_v = ([],[],0.0,0.0)

# set up affine part
dg_helmholtz_bc!(bᵘ, zero_value, field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
dg_helmholtz_bc!(bᵛ, zero_value, field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)

#
rhsᵘ = -1 .* mesh.J .* (mesh.M * ũ ./ (ν*Δt) ) - bᵘ
rhsᵘ *= -1.0 #cholesky nonsense
# then v
rhsᵛ = -1 .* mesh.J .* (mesh.M * ṽ ./ (ν*Δt)) - bᵛ
rhsᵛ *= -1.0 #cholesky nonsense

# step one solve helmholtz equation for velocity field
u¹ = reshape(chol_Hᵘ \ rhsᵘ[:], size(mesh.x) )
v¹ = reshape(chol_Hᵛ \ rhsᵛ[:], size(mesh.x) )

# now set old values equal to new values
@. u⁰ = u¹
@. v⁰ = v¹

# check
u_exact = eval_grid(u_analytic, mesh, t)
v_exact = eval_grid(v_analytic, mesh, t)


println("before satisfying boundary conditions")
u_error = rel_error(u_exact, ũ)
v_error = rel_error(v_exact, ṽ)
println("The relative error is $(u_error)")
println("The relative error is $(v_error)")
println("with satisfying boudnary conditions")
u_error = rel_error(u_exact, u¹)
v_error = rel_error(v_exact, v¹)
println("The relative error is $(u_error)")
println("The relative error is $(v_error)")

println("with bc and 1 norm")
u_error = rel_1_error(u_exact, u¹)
v_error = rel_1_error(v_exact, v¹)
println("The relative error is $(u_error)")
println("The relative error is $(v_error)")

println("relative error in boundary conditions")
println("before")
println(rel_error(u_exact[mesh.vmapB], ũ[mesh.vmapB]))
println("after")
println(rel_error(u_exact[mesh.vmapB], u¹[mesh.vmapB]))
∇⨀!(rhsᵖ , u¹, v¹, mesh)
println("The maximum incompressibility is now $(maximum(abs.(rhsᵖ)))")

∇⨀!(rhsᵖ , ũ, ṽ, mesh)
println("The maximum incompressibility before was $(maximum(abs.(rhsᵖ)))")


rightwall = findall(mesh.x[:] .≈ maximum(mesh.x))
leftwall = findall(mesh.x[:] .≈ minimum(mesh.x))

topwall = findall(mesh.y[:] .≈ maximum(mesh.y))
bottomwall = findall(mesh.y[:] .≈ minimum(mesh.y))

println("------------")

#=
gr()
camera_top = 90 #this is a very hacky way to get a 2D contour plot
camera_side = 0
x = mesh.x;
y = mesh.y;
u = copy(u_exact .- u¹);
v = copy(v_exact .- v¹);
p1 = surface(x[:],y[:],u[:], camera = (camera_side,camera_top))
p2 = surface(x[:],y[:],v[:], camera = (camera_side,camera_top))
plot(p1,p2)
=#
