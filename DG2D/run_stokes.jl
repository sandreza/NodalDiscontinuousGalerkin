
using BenchmarkTools
using LinearAlgebra
using Plots

const ν = 1e-2
#potentially need to debug reference pressure

include("mesh2D.jl")
include("dg_advection.jl")
include("../DG2D/triangles.jl")
include("../DG2D/dg_helmholtz.jl")
include("../random/navier_stokes_structs.jl")
include("../DG2D/dg_stokes.jl")

timings = false
plotting_matrix = false
eigenvalues = false
check_correctness = true
plotting_solution = true
# simulation parameters and grid
n = 1
FileName = "Maxwell0125.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)
ι = field
ns = ns_fields(mesh)

# location of boundary grid points for dirichlet bc
bc = (mesh.vmapB, mesh.mapB)
dbc = ([],[])


#compute tau and define γ
γ = 000.0
τ = compute_τ(mesh)
params = [τ, γ]

#
t = 0.0
# may need to augment with derivative boundary conditions
#bc_u, dbc_u, bc_v, dbc_v = calculate_pearson_bc_vel(mesh, t)

# boundary conditions for velocity, for whatever reason plugging in zero yields faster results. Could be the case that one sets zero to create the operator, then creates new boundary conditions later and evaluates the operator at that specific point
bc_u = (mesh.vmapB,mesh.mapB, 0.0)
#dbc_u = (mesh.vmapB,mesh.mapB, 0.0, 0.0)
dbc_u = ([],[],0.0,0.0)


bc_v = (mesh.vmapB,mesh.mapB, 0.0)
#dbc_v = (mesh.vmapB,mesh.mapB, 0.0,0.0)
dbc_v = ([],[],0.0,0.0)

# no boundary conditions for pressure
bc_p = ([],[],0.0)
dbc_p = ([],[],0.0,0.0)

bc_ϕ = (bc_u, bc_v, bc_p)
dbc_ϕ = (dbc_u, dbc_v, dbc_p)

m,n = size(mesh.x)
ϕ = zeros((m,n,3))
𝒮ϕ = zeros((m,n,3))

dg_stokes_bc!(𝒮ϕ, ϕ, ns, params, mesh, bc!, bc_ϕ, bc_∇!, dbc_ϕ)


𝒮, b = stokes_setup_bc(ϕ, ns, params, mesh, bc!, bc_ϕ, bc_∇!, dbc_ϕ)


# now modify the operator
n𝒮, nb = modify_stokes_operator(𝒮, b);
display(spy(n𝒮))
S = (𝒮 + 𝒮')/2
A = (𝒮 - 𝒮')/2

dropϵzeros!(S, 1e-14)
dropϵzeros!(A, 1e-14)

if eigenvalues
    λ = eigvals(Array(𝒮));

    # check to see that its not invertible
    min_λ = minimum(abs.(λ))
    println("the smallest eigenvalue (in absolute value) is $(min_λ) for the non-invertible operator")


    # check to see that its invertible
    nλ = eigvals(Array(n𝒮));
    min_nλ = minimum(abs.(nλ))
    println("the smallest eigenvalue (in absolute value) is $(min_nλ) for the invertible operator")
end

println("computing the LU decompositions")
lu_𝒮 = lu(𝒮)
lu_n𝒮 = lu(n𝒮)
println("done computing the LU decompositions")

m, n = size(𝒮)

if timings
    println("The amount of time to invert the lu decomposition for an $(m)x$(n) matrix is")
    @btime lu_𝒮 \ b[:];

    println("The amount of time to invert the lu decomposition for an $(m+1)x$(n+1) matrix is")
    @btime lu_n𝒮 \ nb[:];
end
# create a forcing function and check to see that it passes some sanity checks

# if the forcing is a gradient the the velocity fields should be zero
fx(x,y,t) = cos(2π * x) * sin(2π * y)
fy(x,y,t) = sin(2π * x) * cos(2π * y)
p_exact(x,y,t) = -sin(2π * x) * sin(2π * y) / (2π)

c1 = eval_grid(fx, mesh, 0)
c2 = eval_grid(fy, mesh, 0)
c3 = eval_grid(p_exact, mesh, 0)

c1 = mesh.J .* (mesh.M * c1)
c2 = mesh.J .* (mesh.M * c2)

m,n = size(mesh.x)
rhs = zeros(m,n,3)
@. rhs[:,:,1] = c1
@. rhs[:,:,2] = c2
@. rhs[:,:,3] = 0.0

# set rhs
@. nb[1:end-1] = rhs[:]

sol = lu_n𝒮 \ nb[:];

sol = reshape(sol[1:end-1], m, n, 3)

u_computed = copy(sol[:,:,1])
v_computed = copy(sol[:,:,2])
p_computed = copy(sol[:,:,3])

u_error = maximum(abs.(u_computed))
v_error = maximum(abs.(v_computed))
gauge = sum(c3) / length(c3) #
@. p_computed += gauge
p_error = rel_error(c3, p_computed)

println("the error for the u-velocity is $(u_error)")
println("the error for the v-velocity is $(v_error)")
println("the error for the pressure is $(p_error)")

println("The largest jump in pressure is ")
println(maximum(abs.(p_computed[mesh.vmapM]-p_computed[mesh.vmapP])))

tmp = similar(u_computed)
tmpx = similar(u_computed)
tmpy = similar(u_computed)
∇⨀!(tmp, u_computed, v_computed, mesh)
println("The amount of incompressibility is ")
println(maximum(abs.(tmp)))

∇!(tmpx, tmpy, u_computed, mesh)
println("The error in computing the x-derivative u-velocity is ")
println(maximum(abs.(tmpx)))
println("The error in computing the y-derivative u-velocity is ")
println(maximum(abs.(tmpx)))


∇!(tmpx, tmpy, v_computed, mesh)
println("The error in computing the x-derivative v-velocity is ")
println(maximum(abs.(tmpx)))
println("The error in computing the y-derivative v-velocity is ")
println(maximum(abs.(tmpx)))

#=


thing = log.(abs.(p_computed[:] -c3[:]))
thing = p_computed[:] -c3[:]
thing = u_computed
p3 = surface(mesh.x[:],mesh.y[:], thing , camera = (0,90))

∂ˣ_∮(ns.p, mesh, bc!, bc_p)


bc!(ns.p.fⁿ, ns.p.ϕ, bc)


A = 𝒮- 𝒮'
dropϵzeros!(A, 1e-14)
=#
