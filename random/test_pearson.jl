# To make easy comparisons to the Matlab code

include("../DG2D/dg_navier_stokes.jl")
include("../DG2D/dg_poisson.jl")
include("../DG2D/triangles.jl")
include("../DG2D/mesh2D.jl")

# define polynomial order
n = 4

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

# will always need this
inflow_index = findall(bc_label .== "In")
mapI = mapT[inflow_index][1]
vmapI = vmapT[inflow_index][1]
outflow_index = findall(bc_label .== "Out")
mapO = mapT[outflow_index][1]
vmapO = vmapT[outflow_index][1]

# convenience variables
xO = mesh.x[vmapO];
yO = mesh.y[vmapO];
nxO = mesh.nx[mapO];
nyO = mesh.ny[mapO];
xI = mesh.x[vmapI];
yI = mesh.y[vmapI];
nxI = mesh.nx[mapI];
nyI = mesh.ny[mapI];

# dirichlet boundary conditions on the inflow
@. ubc[mapI] = u_exact[vmapI];
@. vbc[mapI] = v_exact[vmapI];
@. pbc[mapI] = p_exact[vmapI];
@. undtbc[mapI] = (-nxI * sin(2*pi*yI)+ nyI * sin(2*pi*xI) ) .* exp(-ν*4*π^2*t);

# dirichlet boundary conditions for the pressure at the outflow
@. pbc[mapO] = p_exact[vmapO];

# neuman boundary conditions for the
@. ubc[mapO] = nyO *( ( 2*π ) * (-cos(2*π*yO) * exp(-ν*4*π^2*t) ) );
@. vbc[mapO] = nxO *( ( 2*π ) * ( cos(2*π*xO) * exp(-ν*4*π^2*t) ) );



#dirichlet
function bc!(ϕ, mesh, bc)
    @. ϕ.fⁿ[bc[2]] = ϕ.u[bc[1]]  - bc[3]
    return nothing
end
#neumann
function bc_∇!(ϕ, mesh, bc)
    @. ϕ.fˣ[bc[2]] = ϕ.φˣ[bc[1]] - bc[3]
    @. ϕ.fʸ[bc[2]] = ϕ.φʸ[bc[1]] - bc[4]
    return nothing
end

#=
# set up pressure laplacian
#boundary condition indices
bc = ([vmapI vmapO][:], [mapI mapO][:], pbc)
# location of boundary grid points for neumann bc
dbc = ([], [], pbc)
Δᵖ, bᵖ = poisson_setup_bc(ι.p, params, mesh, bc_p!, bc, bc_∇p!, dbc)
Δᵖ = cholesky(Δᵖ)

# set up u-velocity laplacian
#boundary condition indices
bc = (vmapI, mapI, ubc)
# location of boundary grid points for neumann bc
dbc = (vmapO, mapO, ubc)
#for future updates one just needs to evaluate poisson_bc or helmholtz_bc at zero
Δᵘ, bᵘ = helhomtz_setup_bc(ι.u, params, mesh, bc_u!, bc, bc_∇u!, dbc)
Δᵘ = cholesky(Δᵘ)
# set up v-velocity laplacian
#boundary condition indices
bc = (vmapI, mapI, vbc)
# location of boundary grid points for neumann bc
dbc = (vmapO, mapO, ubc)
Δᵛ, bᵛ = helhomtz_setup_bc(ι.v, params, mesh, bc_v!, bc, bc_∇v!, dbc)
Δᵛ = cholesky(Δᵛ)
=#




# functions
u_analytic(x,y,t) = -sin(2 * π * y ) * exp( - ν * 4 * π^2 * t);
v_analytic(x,y,t) =  sin(2 * π * x ) * exp( - ν * 4 * π^2 * t);
p_analytic(x,y,t) = -cos(2 * π * x ) * cos(2 *π * y) * exp( - ν * 8 *π^2 * t);

#∂ˣ
∂ˣu_analytic(x,y,t) = 0.0;
∂ˣv_analytic(x,y,t) =  2 * π * cos(2 *π * x ) * exp( - ν * 4 * pi^2 * t);
∂ˣp_analytic(x,y,t) = 2 * π * sin(2 *π * x ) * cos(2 *π * y) * exp( - ν * 8 * π^2 * t);

#∂ʸ
∂ʸu_analytic(x,y,t) = - 2 * π * cos(2 *π * y ) * exp( - ν * 4 * pi^2 * t);
∂ʸv_analytic(x,y,t) =  0.0;
∂ʸp_analytic(x,y,t) = 2 * π * cos(2 *π * x ) * sin(2 *π * y) * exp( - ν * 8 * π^2 * t);

u∇ux_analytic(x,y,t) = u_analytic(x,y,t) * ∂ˣu_analytic(x,y,t) + v_analytic(x,y,t) * ∂ʸu_analytic(x,y,t)
u∇uy_analytic(x,y,t) = u_analytic(x,y,t) * ∂ˣv_analytic(x,y,t) + v_analytic(x,y,t) * ∂ʸv_analytic(x,y,t)

function eval_grid(phield, mesh, t)
    tmp = [phield(mesh.x[i],mesh.y[i], t) for i in 1:length(mesh.x) ]
    return reshape(tmp, size(mesh.x))
end

eval_grid(u_analytic, mesh, 0);


# potential struct for navier_stokes

struct dg_field{T}
    ϕ::T
    ϕ̇::T
    ∂ˣ::T
    ∂ʸ::T
    ∂ⁿ::T
    φˣ::T
    φʸ::T
    φⁿ::T
    fˣ::T
    fʸ::T
    fⁿ::T
    """
    dg_field(mesh)

    # Description

        initialize dg struct

    # Arguments

    -   `mesh`: a mesh to compute on

    # Return Values:

    -   `ϕ` : the field to be computed,
    -   `ϕ̇`: numerical solutions for the field
    -   `∂ˣ`: x-component of derivative
    -   `∂ʸ`: y-component of derivative
    -   `∂ⁿ`: normal component of derivative
    -   `φˣ`: x-component of flux
    -   `φʸ`: y-component of flux
    -   `φⁿ`: normal component of flux
    -   `fˣ`: the numerical flux on face in the x-direction for the computation
    -   `fʸ`: the numerical flux on face in the y-direction for the computation
    -   `fⁿ`: the numerical flux on face in the normal direction for the computation

    """
    function dg_field(mesh)
        # set up the solution
        ϕ   = similar(mesh.x)
        ϕ̇   = similar(mesh.x)
        ∂ˣ  = similar(mesh.x)
        ∂ʸ  = similar(mesh.x)
        ∂ⁿ  = similar(mesh.x)
        φˣ  = similar(mesh.x)
        φʸ  = similar(mesh.x)
        φⁿ  = similar(mesh.x)
        fˣ  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        fʸ  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        fⁿ  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        return new{typeof(ϕ)}(ϕ, ϕ̇, ∂ˣ, ∂ʸ, ∂ⁿ, φˣ, φʸ, φⁿ, fˣ, fʸ, fⁿ)
    end
end

struct ns_fields{T}
    u::T
    v::T
    p::T
    """
    ns_field(mesh)

    # Description

        initialize dg struct

    # Arguments

    -   `mesh`: a mesh to compute on

    # Return Values:

    -   `u` : the u-velocity component struct
    -   `v` : the v-velocity component struct
    -   `p` : the pressure struct

    """
    function ns_fields(mesh)
        # set up the solution
        u = dg_field(mesh)
        v = dg_field(mesh)
        p = dg_field(mesh)
        return new{typeof(u)}(u, v, p)
    end
end
#=

# set up functions to evaluate boundary conditions
#dirichlet
function bc_p!(ι, mesh, bc)
    @. ι.p.fⁿ[bc[2]] = ι.p.ϕ[bc[1]]  - bc[3]
    return nothing
end
#neumann
function bc_∇p!(ι, mesh, bc)
    @. ι.p.fˣ[bc[2]] = ι.p.φˣ[bc[1]] - bc[3]
    @. ι.p.fʸ[bc[2]] = ι.p.φʸ[bc[1]] - bc[4]
    return nothing
end

#dirichlet
function bc_u!(ι, mesh, bc)
    @. ι.u.fⁿ[bc[2]] = ι.u.ϕ[bc[1]] - bc[3]
    return nothing
end
#neumann

function bc_∇u!(ι, mesh, bc)
    @. ι.u.fˣ[bc[2]] = ι.u.φˣ[bc[1]] - bc[3]
    @. ι.u.fʸ[bc[2]] = ι.u.φʸ[bc[1]] - bc[4]
    return nothing
end

#dirichlet
function bc_v!(ι, mesh, bc)
    @. ι.v.fⁿ[bc[2]] = ι.v.ϕ[bc[1]] - bc[3]
    return nothing
end
#neumann
function bc_∇v!(ι, mesh, bc)
    @. ι.v.fˣ[bc[2]] = ι.v.φˣ[bc[1]] - bc[3]
    @. ι.v.fʸ[bc[2]] = ι.v.φʸ[bc[1]] - bc[4]
    return nothing
end
=#

# check pointwise convergece of pressure equation
# use exact boundary conditions as the check
# will also check the divergence and everything else

# ∂ᵗ u = -u⋅∇u + ν Δu - ∇p
# ∇⋅u = 0
# from whence
# Δp = ∇ ⋅ ( -u⋅∇u ), we are checking this equation

println("the size of the solution is $(length(mesh.x))")
println("------------------")
# first compute the advective term
t = 0
# u component set
tmp = eval_grid(u_analytic, mesh, t)
@. ι.u.ϕ = tmp
# v component set
tmp = eval_grid(v_analytic, mesh, t)
@. ι.v.ϕ = tmp
# p component set
tmp = eval_grid(p_analytic, mesh, t)
@. ι.p.ϕ = tmp

# compute advection
sym_advec!(ι.u.φⁿ, ι.u.ϕ, ι.v.ϕ, ι.u.ϕ, mesh)
sym_advec!(ι.v.φⁿ, ι.u.ϕ, ι.v.ϕ, ι.v.ϕ, mesh)

# compute advection analytically
advecu = eval_grid(u∇ux_analytic, mesh, t)
advecv = eval_grid(u∇uy_analytic, mesh, t)

# state
relu = rel_error(advecu, ι.u.φⁿ)
relv = rel_error(advecv, ι.v.φⁿ)
println("The error in computing the advection for u is $(relu)")
println("The error in computing the advection for v is $(relv)")

# compute divergence of advection
rhs = similar(ι.p.ϕ)
∇⨀!(rhs , ι.u.φⁿ, ι.v.φⁿ, mesh)
@. rhs *= -1.0 # since its the negative divergence that shows up

# set up boundary conditions for pressure
# location of boundary grid points for dirichlet bc
dirichlet_pressure_bc = ι.p.ϕ[mesh.vmapB];
bc = (mesh.vmapB, mesh.mapB, dirichlet_pressure_bc)
dbc = ([],[],0.0,0.0)

# set up τ matrix
τ = compute_τ(mesh)
params = [τ]

# set up matrix and affine component
Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc, bc_∇!, dbc)

# set up appropriate rhs
frhsᵖ = mesh.J .* (mesh.M * rhs) - bᵖ
@. frhsᵖ *= -1.0
# cholesky decomposition
Δᵖ = -(Δᵖ + Δᵖ')/2
Δᵖ = cholesky(Δᵖ)

# compute answer
num_solᵖ = Δᵖ \ frhsᵖ[:];

# compute analytic answer
# p component set
tmp = eval_grid(p_analytic, mesh, t)
@. ι.p.ϕ = tmp

# check answer
w2inf = maximum(abs.(ι.p.ϕ[:] .- num_solᵖ)) / maximum(abs.(ι.p.ϕ))
println("The relative error in computing the solution is $(w2inf)")
println("----------------")


Δt = 0.1
# compute operators for helmholtz
t = 0

u_exact = eval_grid(u_analytic, mesh, t)
v_exact = eval_grid(v_analytic, mesh, t)

dirichlet_u_bc = u_exact[mesh.vmapB];
bc_u = (mesh.vmapB, mesh.mapB, dirichlet_u_bc)
dbc_u = ([],[],0.0,0.0)
dirichlet_v_bc = v_exact[mesh.vmapB];
bc_v = (mesh.vmapB, mesh.mapB, dirichlet_v_bc)
dbc_v = ([],[],0.0,0.0)

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

bc_wierd = ([mesh.vmapB[1]], [mesh.mapB[1]], 0.0)
dbc_wierd = (mesh.vmapB[2:end], mesh.mapB[2:end], 0.0, 0.0)


# set up τ matrix
τ = compute_τ(mesh)
params = [τ]

# set up matrix and affine component
Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_wierd, bc_∇!, dbc_wierd)
sΔᵖ, sbᵖ = poisson_setup_bc(field, params, mesh, bc!, bc, bc_∇!, dbc)
#slight regularization
Δᵖ = (Δᵖ + Δᵖ' ) ./ 2
dropϵzeros!(Δᵖ)
chol_Δᵖ = cholesky(-Δᵖ)
#Δᵖ  = lu(-Δᵖ)


# compute forcing term for helmholtz
# first u
sym_advec!(ι.u.φⁿ, ι.u.ϕ, ι.v.ϕ, ι.u.ϕ, mesh)
@. ι.u.φⁿ *= -1 / ν
@. ι.u.φⁿ -= u⁰ / ( ν * Δt )
rhsᵘ = mesh.J .* (mesh.M * ι.u.φⁿ) - bᵘ
rhsᵘ *= -1.0 #cholesky nonsense
# then v
sym_advec!(ι.v.φⁿ, ι.u.ϕ, ι.v.ϕ, ι.v.ϕ, mesh)
@. ι.v.φⁿ *= -1 / ν
@. ι.v.φⁿ -= v⁰ / ( ν * Δt )
rhsᵛ = mesh.J .* (mesh.M * ι.v.φⁿ) - bᵛ
rhsᵛ *= -1.0 #cholesky nonsense

# step one solve helmholtz equation for velocity field
ũ = reshape(chol_Hᵘ \ rhsᵘ[:], size(mesh.x) )
ṽ = reshape(chol_Hᵛ \ rhsᵛ[:], size(mesh.x) )

# step two, project
rhsᵖ = similar(ι.p.ϕ)
∇⨀!(rhsᵖ , ũ, ṽ, mesh)

frhsᵖ = mesh.J .* (mesh.M * rhsᵖ) - bᵖ
@. frhsᵖ *= -1.0
p = reshape(chol_Δᵖ \ frhsᵖ[:], size(mesh.x));
∇!(ι.p.∂ˣ,ι.p.∂ʸ, p, mesh)
u¹ = ũ - ι.p.∂ˣ
v¹ = ṽ - ι.p.∂ʸ

# now set old values equal to new values
@. u⁰ = u¹
@. v⁰ = v¹

# check
u_exact = eval_grid(u_analytic, mesh, t)
v_exact = eval_grid(v_analytic, mesh, t)


println("without incompressibility")
u_error = rel_error(u_exact, ũ)
v_error = rel_error(v_exact, ṽ)
println("The relative error is $(u_error)")
println("The relative error is $(v_error)")
println("with incompressibility")
u_error = rel_error(u_exact, u¹)
v_error = rel_error(v_exact, v¹)
println("The relative error is $(u_error)")
println("The relative error is $(v_error)")

println("with incompressibility and 1 norm")
u_error = rel_1_error(u_exact, u¹)
v_error = rel_1_error(v_exact, v¹)
println("The relative error is $(u_error)")
println("The relative error is $(v_error)")

println("relative error in boundary conditions")
println(rel_error(u_exact[mesh.vmapB], ũ[mesh.vmapB]))
println(rel_error(u_exact[mesh.vmapB], u¹[mesh.vmapB]))
∇⨀!(rhsᵖ , u¹, v¹, mesh)
println("The maximum incompressibility is $(maximum(abs.(rhsᵖ)))")

∇⨀!(rhsᵖ , ũ, ṽ, mesh)
println("The maximum incompressibility before was $(maximum(abs.(rhsᵖ)))")


rightwall = findall(mesh.x[:] .≈ maximum(mesh.x))
leftwall = findall(mesh.x[:] .≈ minimum(mesh.x))

topwall = findall(mesh.y[:] .≈ maximum(mesh.y))
bottomwall = findall(mesh.y[:] .≈ minimum(mesh.y))


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
