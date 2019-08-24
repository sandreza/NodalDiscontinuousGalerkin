#using Plots
using BenchmarkTools
include("../DG2D/dg_navier_stokes.jl")
include("../DG2D/mesh2D.jl")
include("../DG2D/utils2D.jl")
include("../random/navier_stokes_structs.jl")
include("../DG2D/dg_poisson.jl")
include("../DG2D/dg_helmholtz.jl")
include("../DG2D/triangles.jl")

# taking away the lift term in the pressure solve seems to increase accuracy
# only for large elements actually, it seems like the lift operators are
# super important when one has a large number of elements, less so with less elements

# define polynomial order, n=11 is about the right size
n = 10 #even good odd okay
plotting = false
timings = false
stommel = true
const second_order = true
const debug = false
# load grids
#FileName = "pvortex4A01.neu"
FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName

# set up structs
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)

ω = dg_field(mesh)
ωˣ = dg_field(mesh)
ωʸ = dg_field(mesh)

ι = ns_fields(mesh)

# construct boundary data
Nv, VX, VY, K, EtoV, bctype, bc_name = meshreader_gambit_bc_2D(filename)
mapT, vmapT, bc_label = build_bc_maps(mesh, bctype, bc_name)

# set time and time step and viscocity
t = 0.0
t_list = [t]
const Δt = 1e-3 # 1e-3 seems good for convergence tests
const ν  = 1e-2

# set up the Helmholtz and Poisson operators

# calculate boundary condition map
bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)
bc_p, dbc_p = calculate_pearson_bc_p(mesh) #homogenous for pressure

if stommel
    bc_u, dbc_u, bc_v, dbc_v = calculate_stommel_bc_vel(mesh, t)
end

if second_order
    τ = compute_τ(mesh)
    γ = 1 / (ν * Δt / 2) #crank nicholson
    const params_vel = [τ, γ]
    Hᵘ, bᵘ = helmholtz_setup_bc(field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
    Hᵛ, bᵛ = helmholtz_setup_bc(field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)
    Hᵘ = (Hᵘ + Hᵘ') / 2
    Hᵛ = (Hᵛ + Hᵛ') / 2
    dropϵzeros!(Hᵘ)
    dropϵzeros!(Hᵛ)
    chol_Hᵘ = cholesky(-Hᵘ)
    chol_Hᵛ = cholesky(-Hᵛ)

    # set up operators for p
    params = [τ]
    # set up matrix and affine component
    #Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_wierd, bc_∇!, dbc_wierd)
    Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_p, bc_∇!, dbc_p)
    # augments for neumann boundary conditions to make solution unique
    chol_Δᵖ = modify_pressure_Δ(Δᵖ)
    #this is actually the LU factorization, misnamed for sure...
    #lu makes things about a factor of 60 slower ...

    Δᵘ, bᵘ = poisson_setup_bc(field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
    Δᵛ, bᵛ = poisson_setup_bc(field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)
    Δᵘ = (Δᵘ + Δᵘ') / 2
    Δᵛ = (Δᵛ + Δᵛ') / 2
    dropϵzeros!(Δᵘ)
    dropϵzeros!(Δᵛ)
    chol_Δᵘ = cholesky(-Δᵘ)
    chol_Δᵛ = cholesky(-Δᵛ)
else
    # set up operators for u and v
    τ = compute_τ(mesh)
    γ = 1 / (ν * Δt)
    const params_vel = [τ, γ]
    Hᵘ, bᵘ = helmholtz_setup_bc(field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
    Hᵛ, bᵛ = helmholtz_setup_bc(field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)
    Hᵘ = (Hᵘ + Hᵘ') / 2
    Hᵛ = (Hᵛ + Hᵛ') / 2
    dropϵzeros!(Hᵘ)
    dropϵzeros!(Hᵛ)
    chol_Hᵘ = cholesky(-Hᵘ)
    chol_Hᵛ = cholesky(-Hᵛ)

    # set up operators for p
    params = [τ]
    # set up matrix and affine component
    #Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_wierd, bc_∇!, dbc_wierd)
    Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_p, bc_∇!, dbc_p)
    # augments for neumann boundary conditions to make solution unique
    chol_Δᵖ = modify_pressure_Δ(Δᵖ)
    #this is actually the LU factorization, misnamed for sure...
    #lu makes things about a factor of 60 slower ...

    Δᵘ, bᵘ = poisson_setup_bc(field, params_vel, mesh, bc!, bc_u, bc_∇!, dbc_u)
    Δᵛ, bᵛ = poisson_setup_bc(field, params_vel, mesh, bc!, bc_v, bc_∇!, dbc_v)
    Δᵘ = (Δᵘ + Δᵘ') / 2
    Δᵛ = (Δᵛ + Δᵛ') / 2
    dropϵzeros!(Δᵘ)
    dropϵzeros!(Δᵛ)
    chol_Δᵘ = cholesky(-Δᵘ)
    chol_Δᵛ = cholesky(-Δᵛ)
end



###
# now that the problem has been set up
# we set up the initial condition
u⁰ = eval_grid(u_analytic, mesh, t)
v⁰ = eval_grid(v_analytic, mesh, t)
# projection velocities
ũ = similar(u⁰)
ṽ = similar(v⁰)
# velocity at the next time-step
u¹ = similar(u⁰)
v¹ = similar(v⁰)

# storing right hand sides
bᵘ = similar(u⁰)
bᵛ = similar(u⁰)
bᵖ = similar(u⁰)

∇⨀!(bᵖ, u⁰, v⁰, mesh)
println("The analytically computed divergence is ")
println(maximum(abs.(bᵖ)))

####
# main loop
# step 1: Advection
@. ι.u.ϕ = u⁰
@. ι.v.ϕ = v⁰
bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)
ns_advection!(ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt)

# step 2: Pressure projection
bc_p, dbc_p = calculate_pearson_bc_p(mesh)
ns_projection!(ι, bc_p, dbc_p, chol_Δᵖ, ũ, ṽ, bᵖ, params_vel)
tmp = copy(ṽ)
# now consider next time-step
t += Δt
# step 3: Diffuse
bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)
ns_diffuse!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹,  ũ, ṽ, params_vel)

# step 4: set new value of velocity
@. u⁰ = u¹
@. v⁰ = v¹

# check correctness
ns_pearson_check(ι, mesh, t, u¹, v¹, ũ, ṽ)

####

# WARNING lift may have wrong sign in the velocity equation
t = 0.0
t_list = [t]
# we set up the initial condition, pearson
u⁰ = eval_grid(u_analytic, mesh, t)
v⁰ = eval_grid(v_analytic, mesh, t)
if stommel
    u⁰ = eval_grid(u_stommel, mesh, t)
    v⁰ = eval_grid(v_stommel, mesh, t)
end
# projection velocities
ũ = similar(u⁰)
ṽ = similar(v⁰)
# velocity at the next time-step
u¹ = similar(u⁰)
v¹ = similar(v⁰)

# storing right hand sides
bᵘ = similar(u⁰)
bᵛ = similar(u⁰)
bᵖ = similar(u⁰)
f = 10^(-1) * 1
println("maxval + incomp")
println(maximum(abs.(u⁰)))
println(maximum(abs.(v⁰)))
divu  = similar(u⁰)
∇⨀!(divu, u⁰, v⁰, mesh)
thing = abs.(divu[:])
println(maximum(thing))
# check the timestep
times = 1:2000
for i in times
    # pressure on lin 293 and 294 is multiplied by zero for bc
    #ns_timestep!(u⁰, v⁰, u¹, v¹, ũ, ṽ, ν, Δt, ι, mesh, bᵘ, bᵛ, bᵖ, t_list)
    #ns_timestep_other!(u⁰, v⁰, u¹, v¹, ũ, ṽ, ν, Δt, ι, mesh, bᵘ, bᵛ, bᵖ, t_list)
    ns_timestep_stommel!(f, u⁰, v⁰, u¹, v¹, ũ, ṽ, ν, Δt, ι, mesh, bᵘ, bᵛ, bᵖ, t_list, bc_u, bc_v, dbc_u, dbc_v)
    #println("time is $(t_list[1])")
    #ns_pearson_check(ι, mesh, t_list[1], u¹, v¹, ũ, ṽ)
    if plotting
        divu = similar(mesh.x)
        ∇⨀!(divu, u¹, v¹, mesh)
        thing = log.( abs.(divu[:]))
        u_exact = eval_grid(u_analytic, mesh, t_list[1])
        thing = log.(abs.(u¹[:] - u_exact[:]))
        p3 = surface(mesh.x[:],mesh.y[:], thing , camera = (0,90))
        #surface(mesh.x[:],mesh.y[:], px_exact - ι.p.∂ˣ ./ Δt  , camera = (0,90))
        display(p3)
        sleep(0.1)

    end
    #=
    println("The maximum value of the solution is")
    println(maximum(abs.(u⁰)))
    println(maximum(abs.(v⁰)))
    println("---------------")
    divu  = similar(u¹)
    ∇⨀!(divu, ũ, ṽ, mesh)
    ∇⨀!(divu, u¹, v¹, mesh)
    thing = abs.(divu[:])
    println("The incompressibility is")
    println(maximum(thing))
    println("The average incompressibility is ")
    println(mean(thing))
    println("---------------")
    =#
end
println("---------------")
∇⨀!(divu, u¹, v¹, mesh)
thing = abs.(divu[:])
println("The incompressibility is")
println(maximum(thing))
println("The average incompressibility is ")
println(mean(thing))
println("---------------")

println(maximum(abs.(u⁰)))
println(maximum(abs.(v⁰)))
#ns_pearson_check(ι, mesh, t_list[1], u¹, v¹, ũ, ṽ)
#tmpfx = copy(ι.u.φⁿ)
#tmpfy = copy(ι.v.φⁿ)

if timings
    println("computing one time-step takes (fully computing pressure bc)")
    @btime ns_timestep_other!(u⁰, v⁰, u¹, v¹, ũ, ṽ, ν, Δt, ι, mesh, bᵘ, bᵛ, bᵖ, t_list);
    println("computing one time-step takes")
    @btime ns_timestep!(u⁰, v⁰, u¹, v¹, ũ, ṽ, ν, Δt, ι, mesh, bᵘ, bᵛ, bᵖ, t_list);
    println("the pressure solve takes")
    tmp = randn(length(mesh.x)+1);
    @btime chol_Δᵖ \ tmp;
    println("The diffusion solve takes")
    tmp = randn(length(mesh.x));
    @btime chol_Hᵘ \ tmp;
    println("The advection step takes")
    @btime ns_advection!(ι, bc_u, bc_v, mesh, u⁰, v⁰, Δt);
    println("The pressure step takes")
    @btime ns_projection!(ι, bc_p, dbc_p, chol_Δᵖ, ũ, ṽ, bᵖ, params_vel);
    println("The diffusion step takes")
    @btime ns_diffuse!(ι, mesh, bc_u, bc_v, dbc_u, dbc_v, ν, Δt, bᵘ, bᵛ, u¹, v¹,  ũ, ṽ, params_vel)
end


if plotting
    divu = similar(mesh.x)
    ∇⨀!(divu, u¹, v¹, mesh)
    thing = log.( abs.(divu[:])) / log(10)
    thing = v¹[:]
    thing = ũ
    thing = v⁰[:]
    thing = u⁰[:]
    thing = divu
    p3 = surface(mesh.x[:],mesh.y[:], thing , camera = (0,90))
    p3 = surface(mesh.x[:],mesh.y[:], thing , camera = (0,90), clims = (-10,-5))

    thing1 = v⁰[:]
    thing2 = u⁰[:]
    p1 = surface(mesh.x[:],mesh.y[:], thing1 , camera = (0,90))
    p2 = surface(mesh.x[:],mesh.y[:], thing2 , camera = (0,90))
    display(plot(p1,p2))

end
