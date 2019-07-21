# To make easy comparisons to the Matlab code
include("../DG2D/dg_navier_stokes.jl")
include("../DG2D/dg_poisson.jl")
include("../DG2D/triangles.jl")
include("../DG2D/mesh2D.jl")

# define polynomial order
n = 1

# load grids
FileName = "pvortex4A01.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)

# construct boundary data
Nv, VX, VY, K, EtoV, bctype, bc_name = meshreader_gambit_bc_2D(filename)
mapT, vmapT, bc_label = build_bc_maps(mesh, bctype, bc_name)

# construct exact solution
t = 0
ν = 1e-2
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
@. ubc[mapO] = nyO *( ( 2*π) * (-cos(2*π*yO) * exp(-ν*4*π^2*t)));
@. vbc[mapO] = nxO *(( 2*π)*( cos(2*π*xO)*exp(-ν*4*π^2*t)));


# set up functions to evaluate boundary conditions
#homogenous dirichlet
function bc_p!(du, u, bc)
    @. du[bc[2]] = u[bc[1]]  - bc[3]
    return nothing
end
#homogenous neumann
function bc_∇p!(fˣ, fʸ, φˣ, φʸ, bc)
    @. fˣ[bc[2]] = φˣ[bc[1]] - bc[3]
    @. fʸ[bc[2]] = φʸ[bc[1]] - bc[4]
    return nothing
end

#homogenous dirichlet
function bc_u!(du, u, bc)
    @. du[bc[2]] = 2 * u[bc[1]] - bc[3]
    return nothing
end
#homogenous neumann
function bc_∇u!(fˣ, fʸ, φˣ, φʸ, bc)
    @. fˣ[bc[2]] = 2 * φˣ[bc[1]]
    @. fʸ[bc[2]] = 2 * φʸ[bc[1]]
    return nothing
end

#homogenous dirichlet
function bc_v!(du, u, bc)
    @. du[bc[2]] = 2 * u[bc[1]] - bc[3]
    return nothing
end
#homogenous neumann
function bc_∇v!(fˣ, fʸ, φˣ, φʸ, bc)
    @. fˣ[bc[2]] = 2 * φˣ[bc[1]]
    @. fʸ[bc[2]] = 2 * φʸ[bc[1]]
    return nothing
end

# set up pressure laplacian
#boundary condition indices
bc = ([vmapI vmapO][:], [mapI mapO][:], pbc)
# location of boundary grid points for neumann bc
dbc = ([], [], pbc)
Δᵖ = poisson_setup(ι, params, mesh, bc_p!, bc, bc_∇p!, dbc)
Δᵖ = cholesky(Δᵖ)

# set up u-velocity laplacian
#boundary condition indices
bc = (vmapI, mapI, ubc)
# location of boundary grid points for neumann bc
dbc = (vmapO, mapO, ubc)
Δᵘ = poisson_setup(ι, params, mesh, bc_u!, bc, bc_∇u!, dbc)
Δᵘ = cholesky(Δᵘ)
# set up v-velocity laplacian
#boundary condition indices
bc = (vmapI, mapI, vbc)
# location of boundary grid points for neumann bc
dbc = (vmapO, mapO, ubc)
Δᵛ = poisson_setup(ι, params, mesh, bc_v!, bc, bc_∇v!, dbc)
Δᵛ = cholesky(Δᵛ)
