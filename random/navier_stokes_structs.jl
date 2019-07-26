
#=
include("../DG2D/dg_navier_stokes.jl")
include("../random/navier_stokes_structs.jl")
include("../DG2D/dg_poisson.jl")
include("../DG2D/dg_helmholtz.jl")
include("../DG2D/triangles.jl")
include("../DG2D/mesh2D.jl")
include("../DG2D/utils2D.jl")
=#

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




# exact answer pearson_vortex

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

#∂ᵗ
∂ᵗu_analytic(x,y,t) = -sin(2 * π * y ) * exp( - ν * 4 * π^2 * t) * (- ν * 4 * π^2);
∂ᵗv_analytic(x,y,t) =  sin(2 * π * x ) * exp( - ν * 4 * π^2 * t) * (- ν * 4 * π^2);
∂ᵗp_analytic(x,y,t) = -cos(2 * π * x ) * cos(2 *π * y) * exp( - ν * 8 *π^2 * t) * ( - ν * 8 *π^2 );

u∇ux_analytic(x,y,t) = u_analytic(x,y,t) * ∂ˣu_analytic(x,y,t) + v_analytic(x,y,t) * ∂ʸu_analytic(x,y,t)
u∇uy_analytic(x,y,t) = u_analytic(x,y,t) * ∂ˣv_analytic(x,y,t) + v_analytic(x,y,t) * ∂ʸv_analytic(x,y,t)

function eval_grid(phield, mesh, t)
    tmp = [phield(mesh.x[i],mesh.y[i], t) for i in 1:length(mesh.x) ]
    return reshape(tmp, size(mesh.x))
end



# super inefficient, only need points on boundary yet things are evaluated everywhere
function compute_pressure_terms(u⁰, v⁰, ν, fu¹, fv¹, t¹, mesh)
    ∂ᵗu¹ = eval_grid(∂ᵗu_analytic, mesh, t¹)
    ∂ᵗv¹ = eval_grid(∂ᵗv_analytic, mesh, t¹)
    𝒩u = similar(u⁰)
    sym_advec!(𝒩u , u⁰, v⁰, u⁰, mesh)
    𝒩v = similar(v⁰)
    sym_advec!(𝒩v , u⁰, v⁰, v⁰, mesh)
    tmpu, tmpv = ∇⨂∇⨂(u⁰, v⁰, mesh)
    tmpu *= ν
    tmpv *= ν
    px = @. ∂ᵗu¹ + 𝒩u + tmpu - fu¹
    py = @. ∂ᵗv¹ + 𝒩v + tmpv - fv¹
    return -px, -py
end

#stuff I probably won't need
#=
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



=#


# potential struct for navier_stokes


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

# for checking correctness of operators
#=

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



=#


#=
inflow_index = findall(bc_label .== "In")
mapI = mapT[inflow_index][1]
vmapI = vmapT[inflow_index][1]
outflow_index = findall(bc_label .== "Out")
mapO = mapT[outflow_index][1]
vmapO = vmapT[outflow_index][1]
=#
