
#=
include("../DG2D/dg_navier_stokes.jl")
include("../random/navier_stokes_structs.jl")
include("../DG2D/dg_poisson.jl")
include("../DG2D/dg_helmholtz.jl")
include("../DG2D/triangles.jl")
include("../DG2D/mesh2D.jl")
include("../DG2D/utils2D.jl")
=#

function eval_grid(phield, mesh, t)
    tmp = [phield(mesh.x[i],mesh.y[i], t) for i in 1:length(mesh.x) ]
    return reshape(tmp, size(mesh.x))
end

struct dg_field{T}
    ϕ::T
    ϕ⁺::T
    ϕ⁻::T
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
    fx⁺::T
    fx⁻::T
    fy⁺::T
    fy⁻::T
    u::T
    u̇::T
    """
    dg_field(mesh)

    # Description

        initialize dg struct

    # Arguments

    -   `mesh`: a mesh to compute on

    # Return Values:

    -   `ϕ` : the field to be computed,
    -   `ϕ⁺` : the field to be computed, exterior nodes
    -   `ϕ⁻` : the field to be computed, interior nodes
    -   `ϕ̇`: numerical solutions for the field
    -   `∂ˣ`: x-component of derivative
    -   `∂ʸ`: y-component of derivative
    -   `∂ⁿ`: normal component of derivative
    -   `φˣ`: x-component of flux
    -   `φʸ`: y-component of flux
    -   `φⁿ`: normal component of flux
    -   `fˣ`: the numerical jump in flux on face in the x-direction for the computation
    -   `fʸ`: the numerical jump in flux on face in the y-direction for the computation
    -   `fx⁺`: the numerical flux on interior face in the x-direction for the computation
    -   `fy⁺`: the numerical flux on interior face in the y-direction for the computation
    -   `fx⁻`: the numerical flux on interior face in the x-direction for the computation
    -   `fy⁻`: the numerical flux on exterior face in the y-direction for the computation
    -   `fⁿ`: the numerical jump in flux on face in the normal direction for the computation
    -   `u`: for interaction with old structs
    -   `u̇`: for interaction with old structs

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
        ϕ⁺  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        ϕ⁻  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        fx⁺  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        fx⁻  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        fy⁺  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        fy⁻  = zeros(mesh.nfp * mesh.nFaces, mesh.K)
        u   = similar(mesh.x)
        u̇   = similar(mesh.x)
        return new{typeof(ϕ)}(ϕ, ϕ⁺, ϕ⁻, ϕ̇, ∂ˣ, ∂ʸ, ∂ⁿ, φˣ, φʸ, φⁿ, fˣ, fʸ, fⁿ, fx⁺, fx⁻, fy⁺, fy⁻, u, u̇)
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
# might need to change to ϕ.ϕ

function bc!(ϕ, mesh, bc)
    @. ϕ.fⁿ[bc[2]] = ϕ.u[bc[1]]  - bc[3]
    return nothing
end

function bc2!(ϕ, mesh, bc)
    @. ϕ.fⁿ[bc[2]] = ϕ.ϕ[bc[1]]  - bc[3]
    return nothing
end
#neumann
function bc_∇!(ϕ, mesh, bc)
    @. ϕ.fˣ[bc[2]] = ϕ.φˣ[bc[1]] - bc[3]
    @. ϕ.fʸ[bc[2]] = ϕ.φʸ[bc[1]] - bc[4]
    return nothing
end


# initial condition for stommel gyr
Ψ_stommel(x,y,t) = sin(π * x)^2 * sin(π * y )^2  ;
u_stommel(x,y,t) =  sin(π * x)^2 * sin(π * y ) * cos(π * y) ;
v_stommel(x,y,t) = - sin(π * x) * cos(π * x) * sin(π * y )^2;
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

#Δ
Δu_analytic(x,y,t) = (2 * π )^2 * sin(2 * π * y ) * exp( - ν * 4 * π^2 * t);
Δv_analytic(x,y,t) =  - (2 * π )^2 * sin(2 * π * x ) * exp( - ν * 4 * π^2 * t);
Δp_analytic(x,y,t) = ((2 * π )^2 + (2 * π )^2 ) * cos(2 * π * x ) * cos(2 *π * y) * exp( - ν * 8 *π^2 * t);

#∂ᵗ
∂ᵗu_analytic(x,y,t) = -sin(2 * π * y ) * exp( - ν * 4 * π^2 * t) * (- ν * 4 * π^2);
∂ᵗv_analytic(x,y,t) =  sin(2 * π * x ) * exp( - ν * 4 * π^2 * t) * (- ν * 4 * π^2);
∂ᵗp_analytic(x,y,t) = -cos(2 * π * x ) * cos(2 *π * y) * exp( - ν * 8 *π^2 * t) * ( - ν * 8 * π^2 );

u∇ux_analytic(x,y,t) = u_analytic(x,y,t) * ∂ˣu_analytic(x,y,t) + v_analytic(x,y,t) * ∂ʸu_analytic(x,y,t)
u∇uy_analytic(x,y,t) = u_analytic(x,y,t) * ∂ˣv_analytic(x,y,t) + v_analytic(x,y,t) * ∂ʸv_analytic(x,y,t)


#=
"""
∇⨂∇⨂(ns, mesh)

# Description

- compute curl curl of velocity field and include the lift terms


"""
function ∇⨂∇⨂(ns, ω, mesh)
    # compute ∇


    return tmpu, tmpv
end
=#

# super inefficient, only need points on boundary yet things are evaluated everywhere
function compute_pressure_terms(u⁰, v⁰, ν, fu¹, fv¹, t⁰, mesh)
    ∂ᵗu¹ = eval_grid(∂ᵗu_analytic, mesh, t⁰)
    ∂ᵗv¹ = eval_grid(∂ᵗv_analytic, mesh, t⁰)
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

function compute_ghost_points!(ns, bc_u, bc_v, mesh)
    # compute interior and exterior points for u
    @. ns.u.ϕ⁺[:] = ns.u.ϕ[mesh.vmapP]
    @. ns.u.ϕ⁻[:] = ns.u.ϕ[mesh.vmapM]
    # set the external flux equal to the boundary condition flux
    # this is because we are using a rusonov flux
    @. ns.u.ϕ⁺[mesh.mapB] = bc_u[3]
    # compute interior and exterior points for v
    @. ns.v.ϕ⁺[:] = ns.v.ϕ[mesh.vmapP]
    @. ns.v.ϕ⁻[:] = ns.v.ϕ[mesh.vmapM]
    # set the external flux equal to the boundary condition flux
    # this is because we are using a rusonov flux
    @. ns.v.ϕ⁺[mesh.mapB] = bc_v[3]
    return nothing
end

function compute_surface_fluxes!(ns, mesh)
    # exterior fluxes for u
    @. ns.u.fx⁺ = ns.u.ϕ⁺ * ns.u.ϕ⁺
    @. ns.u.fy⁺ = ns.v.ϕ⁺ * ns.u.ϕ⁺
    # interior fluxes for u
    @. ns.u.fx⁻ = ns.u.ϕ⁻ * ns.u.ϕ⁻
    @. ns.u.fy⁻ = ns.v.ϕ⁻ * ns.u.ϕ⁻
    # exterior fluxes for v
    @. ns.v.fx⁺ = ns.u.ϕ⁺ * ns.v.ϕ⁺
    @. ns.v.fy⁺ = ns.v.ϕ⁺ * ns.v.ϕ⁺
    # interior fluxes for v
    @. ns.v.fx⁻ = ns.u.ϕ⁻ * ns.v.ϕ⁻
    @. ns.v.fy⁻ = ns.v.ϕ⁻ * ns.v.ϕ⁻

    return nothing
end

function compute_maximum_face_velocity(ns, mesh)
    # compute normal velocities
    tmp⁺ = abs.( mesh.nx .* ns.u.ϕ⁺ + mesh.ny .* ns.v.ϕ⁺ )
    tmp⁻ = abs.( mesh.nx .* ns.u.ϕ⁻ + mesh.ny .* ns.v.ϕ⁻ )
    maxtmp = [ maximum([tmp⁻[i] tmp⁺[i]]) for i in 1:length(tmp⁺) ]
    maxface = maximum(reshape(maxtmp,mesh.nfp, mesh.nFaces *  mesh.K), dims = 1);
    maxtmp = reshape(maxtmp, mesh.nfp, mesh.nFaces * mesh.K)
    for j in 1:(mesh.nFaces * mesh.K)
            @. maxtmp[:, j] = maxface[j]
    end
    return reshape(maxtmp, size(ns.u.ϕ⁺))
end

function compute_lift_terms(ns, mesh, maxvel)
    # compute surface flux for u
    @. ns.u.fⁿ = mesh.nx * ( ns.u.fx⁺ - ns.u.fx⁻) + mesh.ny * ( ns.u.fy⁺ - ns.u.fy⁻) + maxvel * (ns.u.ϕ⁻ - ns.u.ϕ⁺)
    # compute lift term for u
    tmpu = mesh.lift * ( mesh.fscale .* ns.u.fⁿ) * 0.5
    # compute surface flux for v
    @. ns.v.fⁿ = mesh.nx * ( ns.v.fx⁺ - ns.v.fx⁻) + mesh.ny * ( ns.v.fy⁺ - ns.v.fy⁻) + maxvel * (ns.v.ϕ⁻ - ns.v.ϕ⁺)
    tmpv = mesh.lift * ( mesh.fscale .* ns.v.fⁿ) * 0.5
    return tmpu, tmpv
end

function compute_div_lift_terms(ι, mesh)
    # compute surface flux for u
    diffs = @. mesh.nx[:] * (ι.u.φⁿ[mesh.vmapP]-ι.u.φⁿ[mesh.vmapM]) + mesh.ny[:] * (ι.v.φⁿ[mesh.vmapP]-ι.v.φⁿ[mesh.vmapM])
    diffs = reshape(diffs, mesh.nFaces *mesh.nfp, mesh.K)
    # compute lift term
    div_lift = mesh.lift * ( mesh.fscale .* diffs) * 0.5

    return div_lift
end

function compute_pressure_lift_terms(ι, mesh)
    # compute surface flux for u
    diffsx = @. mesh.nx[:] * (ι.p.ϕ[mesh.vmapP]-ι.p.ϕ[mesh.vmapM])
    diffsx = reshape(diffsx, mesh.nFaces *mesh.nfp, mesh.K)

    diffsy = @. mesh.ny[:] * (ι.p.ϕ[mesh.vmapP]-ι.p.ϕ[mesh.vmapM])
    diffsy = reshape(diffsy, mesh.nFaces *mesh.nfp, mesh.K)
    # compute lift terms y
    px_lift = mesh.lift * ( mesh.fscale .* diffsx) * 0.5
    py_lift = mesh.lift * ( mesh.fscale .* diffsy) * 0.5

    return px_lift, py_lift
end


#these enter in as a right hand side to the appropriate equations
function calculate_pearson_bc_vel(mesh, t)
    # it is assumed that t refers to time t¹

    # compute u and v boundary conditions (since it is time dependent)
    u_exact = eval_grid(u_analytic, mesh, t)
    v_exact = eval_grid(v_analytic, mesh, t)
    dirichlet_u_bc = u_exact[mesh.vmapB];
    bc_u = (mesh.vmapB, mesh.mapB, dirichlet_u_bc)
    dbc_u = ([],[],0.0,0.0)
    dirichlet_v_bc = v_exact[mesh.vmapB];
    bc_v = (mesh.vmapB, mesh.mapB, dirichlet_v_bc)
    dbc_v = ([],[],0.0,0.0)

    return bc_u, dbc_u, bc_v, dbc_v
end

function calculate_stommel_bc_vel(mesh, t)
    # it is assumed that t refers to time t¹
    # compute u and v boundary conditions (since it is time dependent)
    bc_u = (mesh.vmapB, mesh.mapB, 0.0)
    dbc_u = ([],[],0.0,0.0)
    bc_v = (mesh.vmapB, mesh.mapB, 0.0)
    dbc_v = ([],[],0.0,0.0)
    return bc_u, dbc_u, bc_v, dbc_v
end

function calculate_pearson_bc_p(mesh, t, Δt, ν, u⁰, v⁰)
    # it is assumed that t refers to time t¹

    # compute pressure boundary conditions
    # note that this is a computation over the entire domain
    # we can use this to form the residual to see how well we are satisfying the PDE

    ∂pˣ, ∂pʸ = compute_pressure_terms(u⁰, v⁰, ν, 0.0, 0.0, t, mesh)
    @. ∂pˣ *= Δt
    @. ∂pʸ *= Δt
    # just to make invertible
    bc_p = ([], [], 0.0)
    dbc_p = (mesh.vmapB[1:end], mesh.mapB[1:end], ∂pˣ[mesh.vmapB[1:end]], ∂pʸ[mesh.vmapB[1:end]])
    return bc_p, dbc_p
end

function calculate_pearson_bc_p(mesh)
    # it is assumed that t refers to time t¹

    # compute pressure boundary conditions
    # note that this is a computation over the entire domain
    # we can use this to form the residual to see how well we are satisfying the PDE

    #∂pˣ, ∂pʸ = compute_pressure_terms(u⁰, v⁰, ν, 0.0, 0.0, t-Δt, mesh)
    ∂pˣ = zeros(size(mesh.x))
    ∂pʸ = zeros(size(mesh.x))
    # nuemann boundary conditions for pressure
    bc_p = ([], [], 0.0)
    dbc_p = (mesh.vmapB[1:end], mesh.mapB[1:end], ∂pˣ[mesh.vmapB[1:end]], ∂pʸ[mesh.vmapB[1:end]])
    return bc_p, dbc_p
end


function ns_advection!(ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)
    # compute u,v surface contributions
    compute_ghost_points!(ι, bc_u, bc_v, mesh)
    compute_surface_fluxes!(ι, mesh)
    maxvel = compute_maximum_face_velocity(ι, mesh)
    ∮u, ∮v = compute_lift_terms(ι, mesh, maxvel)
    # println("the jump in flux is $(maximum(abs.(∮u)))")
    # now compute contributions fo each field
    # first u
    sym_advec!(ι.u.φⁿ, u⁰, v⁰, u⁰, mesh)
    @. ι.u.φⁿ += ∮u
    @. ι.u.φⁿ *= -Δt
    @. ι.u.φⁿ += u⁰
    # then v
    sym_advec!(ι.v.φⁿ, u⁰, v⁰, v⁰, mesh)
    @. ι.v.φⁿ += ∮v
    @. ι.v.φⁿ *= -Δt
    @. ι.v.φⁿ += v⁰
    return nothing
end

function ns_stommel!(f, ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)
    @. ι.u.φⁿ += Δt * ( -f * v⁰ + sin(π * mesh.y) )
    @. ι.v.φⁿ += Δt * f * u⁰
end

function ns_stommel_β!(f, β, ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)
    @. ι.u.φⁿ += Δt *  ( (f + β * (mesh.y + 1.0) ) * v⁰ + sin(π * mesh.y ./ 2.0 ) )
    @. ι.v.φⁿ += -Δt * (  f + β * (mesh.y + 1.0) ) * u⁰
end

function ns_projection!(ι, bc_p, dbc_p, chol_Δᵖ, ũ, ṽ, bᵖ, params_vel)
    zero_value = zeros(size(mesh.x))
    dg_poisson_bc!(bᵖ, zero_value, field, params_vel, mesh, bc!, bc_p, bc_∇!, dbc_p)

    # take the divergence of the solution
    rhsᵖ = similar(ι.p.ϕ)
    ∇⨀!(rhsᵖ, ι.u.φⁿ, ι.v.φⁿ, mesh)
    println("The maximum incompressibility of the nonlinear part is")
    println(maximum(abs.(rhsᵖ)))
    #construct appropriate lift!
    ∮∇⨀u = compute_div_lift_terms(ι, mesh)
    @. rhsᵖ += ∮∇⨀u

    # construct the right hand side for poissons equation
    frhsᵖ = mesh.J .* (mesh.M * rhsᵖ) - bᵖ
    @. frhsᵖ *= -1.0
    # since we are imposing average of p is zero
    rhs_p = zeros(length(frhsᵖ)+1)
    tmp = length(frhsᵖ)
    @. rhs_p[1:tmp] = frhsᵖ[:]
    # solve the linear system
    p = reshape( (chol_Δᵖ \ rhs_p[:])[1:length(mesh.x)], size(mesh.x));
    @. ι.p.ϕ = p
    # compute the gradient
    ∇!(ι.p.∂ˣ,ι.p.∂ʸ, p, mesh)

    # compute pressure lift terms
    px_lift, py_lift = compute_pressure_lift_terms(ι, mesh)
    # project
    @. ũ = ι.u.φⁿ - ι.p.∂ˣ - px_lift
    @. ṽ = ι.v.φⁿ - ι.p.∂ʸ - py_lift

    ∇⨀!(rhsᵖ, ũ, ṽ, mesh)
    println("The maximum incompressibility of the nonlinear part is now")
    println(maximum(abs.(rhsᵖ)))

    if second_order
        tmpˣ, tmpʸ = ∇⨂∇⨂(ι.u.φⁿ, ι.v.φⁿ, mesh)
        @. ι.u.u̇ = tmpˣ
        @. ι.v.u̇ = tmpʸ
    end

    return nothing
end

function ns_curl_curl!(ι, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹, ũ, ṽ, chol_Δᵘ, chol_Δᵛ, mesh)
    zero_value = zeros(size(mesh.x))
    tmpˣ, tmpʸ = ∇⨂∇⨂(ι.u.φⁿ, ι.v.φⁿ, mesh)
    @. ι.u.u̇ = tmpˣ
    @. ι.v.u̇ = tmpʸ
    rhsᵖ = similar(ι.p.ϕ)
    ∇⨀!(rhsᵖ, tmpˣ, tmpʸ, mesh)
    #println("The maximum incompressibility of the nonlinear part is")
    #println(maximum(abs.(rhsᵖ)))

    zero_value = zeros(size(mesh.x))
    # set up affine part
    dg_poisson_bc!(bᵘ, zero_value, field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
    dg_poisson_bc!(bᵛ, zero_value, field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)

    #
    rhsᵘ = 1 .* mesh.J .* (mesh.M * tmpˣ) + bᵘ

    # then v
    rhsᵛ = 1 .* mesh.J .* (mesh.M * tmpʸ) + bᵛ


    # step one solve helmholtz equation for velocity field
    tmpu¹ = reshape(chol_Δᵘ \ rhsᵘ[:], size(mesh.x) )
    tmpv¹ = reshape(chol_Δᵛ \ rhsᵛ[:], size(mesh.x) )
    @. ũ = tmpu¹
    @. ṽ = tmpv¹

    ∇⨀!(rhsᵖ, ũ, ṽ, mesh)
    #println("The maximum incompressibility of the nonlinear part is now")
    #println(maximum(abs.(rhsᵖ)))

end

function ns_diffuse!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹, ũ, ṽ, params_vel)
    zero_value = zeros(size(mesh.x))
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
    tmpu¹ = reshape(chol_Hᵘ \ rhsᵘ[:], size(mesh.x) )
    tmpv¹ = reshape(chol_Hᵛ \ rhsᵛ[:], size(mesh.x) )
    @. u¹ = tmpu¹
    @. v¹ = tmpv¹
    return nothing
end


function ns_diffuse_2!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹, ũ, ṽ, params_vel)
    zero_value = zeros(size(mesh.x))
    # set up affine part
    dg_helmholtz_bc!(bᵘ, zero_value, field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
    dg_helmholtz_bc!(bᵛ, zero_value, field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)

    #
    rhsᵘ = -1 .* mesh.J .* (mesh.M *  ( ũ ./ (ν*Δt/2) .- ι.u.u̇ ) ) - bᵘ
    rhsᵘ *= -1.0 #cholesky nonsense
    # then v
    rhsᵛ = -1 .* mesh.J .* (mesh.M * (ṽ ./ (ν*Δt/2) .- ι.v.u̇) ) - bᵛ
    rhsᵛ *= -1.0 #cholesky nonsense

    # step one solve helmholtz equation for velocity field
    tmpu¹ = reshape(chol_Hᵘ \ rhsᵘ[:], size(mesh.x) )
    tmpv¹ = reshape(chol_Hᵛ \ rhsᵛ[:], size(mesh.x) )
    @. u¹ = tmpu¹
    @. v¹ = tmpv¹
    return nothing
end

function ns_pearson_check(ι, mesh, t, u¹, v¹, ũ, ṽ)
    println("-------------------------")
    u_exact = eval_grid(u_analytic, mesh, t)
    v_exact = eval_grid(v_analytic, mesh, t)
    px_exact = eval_grid(∂ˣp_analytic, mesh, t)
    py_exact = eval_grid(∂ʸp_analytic, mesh, t)

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
    tmp = similar(u¹)
    ∇⨀!(tmp , u¹, v¹, mesh)
    println("The maximum incompressibility is now $(maximum(abs.(tmp)))")
    ∇⨀!(tmp , ũ, ṽ, mesh)
    println("The maximum incompressibility before was $(maximum(abs.(tmp)))")

    println("the relative error in computing the pressure gradient is ")
    px_error = rel_error(px_exact, ι.p.∂ˣ ./ Δt)
    py_error = rel_error(py_exact, ι.p.∂ʸ ./ Δt)
    println("The px relative error is $(px_error)")
    println("The py relative error is $(py_error)")
    println(" ")

    println("the maximum discontinuity across gridpoints for u is ")
    jump_max = maximum(abs.(ι.u.ϕ[mesh.vmapP] .- ι.u.ϕ[mesh.vmapM]))
    println(jump_max)
    println("the maximum discontinuity across gridpoints for p is ")
    jump_max = maximum(abs.(ι.p.ϕ[mesh.vmapP] .- ι.p.ϕ[mesh.vmapM]))
    println(jump_max)
    println("the maximum discontinuity across gridpoints for v is ")
    jump_max = maximum(abs.(ι.v.ϕ[mesh.vmapP] .- ι.v.ϕ[mesh.vmapM]))
    println(jump_max)



    xlift = @. mesh.nx[:] * (ι.p.∂ˣ[mesh.vmapP]-ι.p.∂ˣ[mesh.vmapM])
    ylift = @. mesh.ny[:] * (ι.p.∂ʸ[mesh.vmapP]-ι.p.∂ʸ[mesh.vmapM])
    xlift = reshape(xlift, mesh.nFaces * mesh.nfp, mesh.K)
    ylift = reshape(ylift, mesh.nFaces * mesh.nfp, mesh.K)

    ∂ˣq = ι.p.∂ˣ + mesh.lift * ( mesh.fscale .* xlift) * 0.5
    ∂ʸq = ι.p.∂ʸ + mesh.lift * ( mesh.fscale .* ylift) * 0.5

    discontinuity_error = maximum(abs.(∂ˣq[mesh.vmapP] - ∂ˣq[mesh.vmapM]))
    println("The maximum discontinuity in qx is $(discontinuity_error)")

    discontinuity_error = maximum(abs.(∂ʸq[mesh.vmapP] - ∂ʸq[mesh.vmapM]))
    println("The maximum discontinuity in qy is $(discontinuity_error)")
    println("-------------------------")
    return nothing
end

function ns_timestep!(u⁰, v⁰, u¹, v¹, ũ, ṽ, ν, Δt, ι, mesh, bᵘ, bᵛ, bᵖ, t_list)
    t = t_list[1]
    # step 1: Advection
    @. ι.u.ϕ = u⁰
    @. ι.v.ϕ = v⁰
    bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)
    ns_advection!(ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)
    # if you mess up the boundary values you get errors

    # step 2: Pressure projection
    bc_p, dbc_p = calculate_pearson_bc_p(mesh)
    ns_projection!(ι, bc_p, dbc_p, chol_Δᵖ, ũ, ṽ, bᵖ, params_vel)
    # now consider next time-step
    @. t_list += Δt
    t = t_list[1]

    # step 3: Diffuse
    bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)
    ns_diffuse!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹,  ũ, ṽ, params_vel)

    # step 4: set new value of velocity
    @. u⁰ = u¹
    @. v⁰ = v¹
    return nothing
end


function ns_timestep_other!(u⁰, v⁰, u¹, v¹, ũ, ṽ, ν, Δt, ι, mesh, bᵘ, bᵛ, bᵖ, t_list)
    t = t_list[1]
    # step 1: Advection
    @. ι.u.ϕ = u⁰
    @. ι.v.ϕ = v⁰
    bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)
    ns_advection!(ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)
    # if you mess up the boundary values you get errors

    # step 2: Pressure projection
    bc_p, dbc_p = calculate_pearson_bc_p(mesh, t, Δt, ν, u⁰, v⁰)
    fu¹ = 0.0
    fv¹ = 0.0
    compute_pressure_terms(u⁰, v⁰, ν, fu¹, fv¹, t, mesh)
    #ns_projection!(ι, bc_p, dbc_p, chol_Δᵖ, ũ, ṽ, bᵖ, params_vel)
    ns_curl_curl!(ι, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹, ũ, ṽ, chol_Δᵘ, chol_Δᵛ, mesh)
    # now consider next time-step
    @. t_list += Δt
    t = t_list[1]

    # step 3: Diffuse
    bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)
    if second_order
        ns_diffuse_2!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹,  ũ, ṽ, params_vel)
    else
        ns_diffuse!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹,  ũ, ṽ, params_vel)
    end

    # step 4: set new value of velocity
    @. u⁰ = u¹
    @. v⁰ = v¹
    return nothing
end

function modify_pressure_Δ(Δᵖ)
    m,n = size(Δᵖ)
    nΔᵖ = spzeros(m+1,n+1)
    @. nΔᵖ[1:n, 1:m] = Δᵖ
    # pad with enforcement of Lagrange multipliers
    @. nΔᵖ[1:n,m+1] = 1.0
    @. nΔᵖ[n+1,1:m] = 1.0
    dropϵzeros!(nΔᵖ)
    maximum(abs.((nΔᵖ - nΔᵖ' ) ./ 2))
    nΔᵖ = (nΔᵖ + nΔᵖ' ) ./ 2
    dropϵzeros!(nΔᵖ)
    lu_Δᵖ = lu(-nΔᵖ)
    return lu_Δᵖ
end



function ns_timestep_stommel!(f, u⁰, v⁰, u¹, v¹, ũ, ṽ, ν, Δt, ι, mesh, bᵘ, bᵛ, bᵖ, t_list, bc_u, bc_v, dbc_u, dbc_v)
    t = t_list[1]
    # step 1: Advection
    #@. ι.u.ϕ = u⁰
    #@. ι.v.ϕ = v⁰

    #bc_u, dbc_u, bc_v, dbc_v = calculate_stommel_bc_vel(mesh, t)
    #ns_advection!(ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)
    @. ι.u.φⁿ = u⁰
    @. ι.v.φⁿ = v⁰
    # ns_stommel!(f, ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)
    ns_stommel_β!(f, f/3.0, ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)
    # if you mess up the boundary values you get errors

    # step 2: Pressure projection
    #bc_p, dbc_p = calculate_pearson_bc_p(mesh)
    #ns_projection!(ι, bc_p, dbc_p, chol_Δᵖ, ũ, ṽ, bᵖ, params_vel)
    ns_curl_curl!(ι, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹, ũ, ṽ, chol_Δᵘ, chol_Δᵛ, mesh)
    #@. ũ = ι.u.ϕ
    #@. ṽ = ι.v.ϕ

    # now consider next time-step
    @. t_list += Δt
    t = t_list[1]

    # step 3: Diffuse
    bc_u, dbc_u, bc_v, dbc_v = calculate_stommel_bc_vel(mesh, t)
    if second_order
        ns_diffuse_2!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹,  ũ, ṽ, params_vel)
    else
        ns_diffuse!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹,  ũ, ṽ, params_vel)
    end

    # step 4: set new value of velocity
    @. u⁰ = u¹
    @. v⁰ = v¹
    return nothing
end




function ns_curl_curl2!(ι, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹, ũ, ṽ, chol_Δᵘ, chol_Δᵛ, mesh)
    # first compute ω
    #@. ι.v.ϕ = ι.v.φⁿ
    #@. ι.u.ϕ = ι.u.φⁿ
    #tmp1 = ∂ˣ_∮(ι.v, mesh, bc2!, bc_v)
    #tmp2 = ∂ʸ_∮(ι.u, mesh, bc2!, bc_u)
    # could try to include b.c here
    tmp1 = ∂ˣ_∮(ι.v.φⁿ, mesh)
    tmp2 = ∂ʸ_∮(ι.u.φⁿ, mesh)
    ω = tmp1 - tmp2
    # check incompressibility

    tmp1 = ∂ˣ_∮(ω, mesh)
    tmp2 = ∂ʸ_∮(ω, mesh)

    tmpˣ =   tmp2
    tmpʸ = - tmp1

    # save the "laplacian of u and v"
    @. ι.u.u̇ = tmpˣ
    @. ι.v.u̇ = tmpʸ

    rhsᵖ = similar(ι.p.ϕ)
    ∇⨀!(rhsᵖ, tmpˣ, tmpʸ, mesh)


    println("The maximum incompressibility of the nonlinear part is")
    println(maximum(abs.(rhsᵖ)))

    zero_value = zeros(size(mesh.x))
    # set up affine part
    dg_poisson_bc!(bᵘ, zero_value, field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
    dg_poisson_bc!(bᵛ, zero_value, field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)

    #
    rhsᵘ = 1 .* mesh.J .* (mesh.M * tmpˣ) + bᵘ

    # then v
    rhsᵛ = 1 .* mesh.J .* (mesh.M * tmpʸ) + bᵛ


    # step one solve helmholtz equation for velocity field
    tmpu¹ = reshape(chol_Δᵘ \ rhsᵘ[:], size(mesh.x) )
    tmpv¹ = reshape(chol_Δᵛ \ rhsᵛ[:], size(mesh.x) )
    @. ũ = tmpu¹
    @. ṽ = tmpv¹

    ∇⨀!(rhsᵖ, ũ, ṽ, mesh)
    println("The maximum incompressibility of the nonlinear part is now")
    println(maximum(abs.(rhsᵖ)))

end


function ∂ˣ_∮(ι, mesh, bc_ϕ!, bc)
    # form field differnces at faces
    @. ι.fⁿ[:] =  (ι.ϕ[mesh.vmapM] - ι.ϕ[mesh.vmapP]) / 2 #central flux
    # enforce bc
    bc_ϕ!(ι, mesh, bc)
    # compute normal component in the x-direction
    @. ι.fˣ = mesh.nx * ι.fⁿ
    # compute lift term
    liftx = mesh.lift * (mesh.fscale .* ι.fˣ )
    # compute partial with respect to x
    ∇!(ι.∂ˣ, ι.∂ʸ, ι.ϕ, mesh)
    return ι.∂ˣ + liftx
end#

function ∂ʸ_∮(ι, mesh, bc_ϕ!, bc)
    # form field differnces at faces
    @. ι.fⁿ[:] =  (ι.ϕ[mesh.vmapM] - ι.ϕ[mesh.vmapP]) / 2 #central flux
    # enforce bc
    bc_ϕ!(ι, mesh, bc)
    # compute normal component in the x-direction
    @. ι.fʸ = mesh.ny * ι.fⁿ
    # compute lift term
    lifty = mesh.lift * (mesh.fscale .* ι.fʸ )
    # compute partial with respect to x
    ∇!(ι.∂ˣ, ι.∂ʸ, ι.ϕ, mesh)
    return ι.∂ʸ + lifty
end

# no boundary conditions
function ∂ˣ_∮(ϕ, mesh)
    # form field differnces at faces
    fⁿ =  (ϕ[mesh.vmapM] - ϕ[mesh.vmapP]) ./ 2 #central flux
    fˣ = mesh.nx .* reshape(fⁿ, size(mesh.nx))
    # compute lift term
    liftx = mesh.lift * (mesh.fscale .* fˣ )
    # compute partial with respect to x
    ∂ˣ = similar(ϕ)
    ∂ʸ = similar(ϕ)
    ∇!(∂ˣ, ∂ʸ, ϕ, mesh)
    return ∂ˣ + liftx
end#

#no boundary conditions
function ∂ʸ_∮(ϕ, mesh)
    # form field differnces at faces
    fⁿ=  (ϕ[mesh.vmapM] - ϕ[mesh.vmapP]) ./ 2 #central flux
    fʸ = mesh.ny .* reshape(fⁿ, size(mesh.ny))
    # compute lift term
    lifty = mesh.lift * (mesh.fscale .* fʸ )
    # compute partial with respect to x
    ∂ˣ = similar(ϕ)
    ∂ʸ = similar(ϕ)
    ∇!(∂ˣ, ∂ʸ, ϕ, mesh)
    return ∂ʸ + lifty
end

function solve_Ψ(u, v, mesh, Δ)
    ω = ∇⨂(u,v,mesh)
    rhs = mesh.J .* (mesh.M * ω)
    Ψ = Δ \ rhs[:]
    return Ψ

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
