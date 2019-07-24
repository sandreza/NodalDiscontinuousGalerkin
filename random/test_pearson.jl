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
ι = ns_fields(mesh)

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



###

# functions
u_analytic(x,y,t) = -sin(2 * pi * y ) * exp( - ν * 4 * pi^2 * t);
v_analytic(x,y,t) =  sin(2 * pi * x ) * exp( - ν * 4 * pi^2 * t);
p_analytic(x,y,t) = -cos(2 * pi * x ) * cos(2 *π * y) * exp( - ν * 8 *π^2 * t);

#∂ˣ
∂ˣu_analytic(x,y,t) = 0.0;
∂ˣv_analytic(x,y,t) =  2 * π * cos(2 *π * x ) * exp( - ν * 4 * pi^2 * t);
∂ˣp_analytic(x,y,t) = 2 * π * sin(2 *π * x ) * cos(2 *π * y) * exp( - ν * 8 * π^2 * t);

#∂ʸ
∂ʸu_analytic(x,y,t) = - 2 * π * cos(2 *π * y ) * exp( - ν * 4 * pi^2 * t);
∂ʸv_analytic(x,y,t) =  0.0;
∂ʸp_analytic(x,y,t) = 2 * π * cos(2 *π * x ) * sin(2 *π * y) * exp( - ν * 8 * π^2 * t);

function eval_grid(phield, mesh, t)
    tmp = [phield(mesh.x[i],mesh.y[i], t) for i in 1:length(mesh.x) ]
    return reshape(tmp, size(mesh.x))
end

eval_grid(u_analytic, mesh, 0);
###

###
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
###

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
end \phi
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
###
