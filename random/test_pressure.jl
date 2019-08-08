# tests the pressure solve specifically
#using Plots
using BenchmarkTools

include("../DG2D/dg_navier_stokes.jl")
include("../DG2D/mesh2D.jl")
include("../DG2D/utils2D.jl")
include("../random/navier_stokes_structs.jl")
include("../DG2D/dg_poisson.jl")
include("../DG2D/dg_helmholtz.jl")
include("../DG2D/triangles.jl")

# List of notable points: computation of the symmetrized advection is potentially bad
# The previous way of solving for pressure was wrong
# Need something that handles neumann boundary conditions
# one way is to assume that the solution is mean zero
# need to include lift terms for the pressure gradient
# most of the incompressibility error comes from the jump in q_x and q_y
# this suggest that it is much better (for an incompressible model) to just solve
# stoke's equations directly

# define polynomial order, n=11 is about the right size
n = 9
neumann = false
wierd = false
plotting = false
timing = false

const debug = true
const ν = 1e-2

# load grids
#FileName = "pvortex4A01.neu"
FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName

# set up structs
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)
ι = ns_fields(mesh)

# construct boundary data
Nv, VX, VY, K, EtoV, bctype, bc_name = meshreader_gambit_bc_2D(filename)
mapT, vmapT, bc_label = build_bc_maps(mesh, bctype, bc_name)

# set time and time step and viscocity
t = 0.0


# evaluate analytic solution on grid
u_exact = eval_grid(u_analytic, mesh, t)
v_exact = eval_grid(v_analytic, mesh, t)
p_exact = eval_grid(p_analytic, mesh, t)

∂ˣu_exact = eval_grid(∂ˣu_analytic, mesh, t)
∂ˣv_exact = eval_grid(∂ˣv_analytic, mesh, t)
∂ˣp_exact = eval_grid(∂ˣp_analytic, mesh, t)

∂ʸu_exact = eval_grid(∂ʸu_analytic, mesh, t)
∂ʸv_exact = eval_grid(∂ʸv_analytic, mesh, t)
∂ʸp_exact = eval_grid(∂ʸp_analytic, mesh, t)

Δp_exact = eval_grid(Δp_analytic, mesh, t)


# bc_p, dbc_p = calculate_pearson_bc_p(mesh) #homogenous for pressure


# dirichlet
bc_p = (mesh.vmapB[1:end], mesh.mapB[1:end], p_exact[mesh.vmapB[1:end]])
dbc_p = ([],[], 0.0, 0.0)
# neumann, with normal
if neumann
    bc_p = ( [], [], 0.0 )
    dbc_p = (mesh.vmapB[1:end], mesh.mapB[1:end], ∂ˣp_exact[mesh.vmapB[1:end]], ∂ʸp_exact[mesh.vmapB[1:end]])
elseif wierd
    bc_p =  ([mesh.vmapB[1]], [mesh.mapB[1]], 0.0 )
    dbc_p = (mesh.vmapB[2:end], mesh.mapB[2:end], ∂ˣp_exact[mesh.vmapB[2:end]], ∂ʸp_exact[mesh.vmapB[2:end]])
end



# set up operators for u and v
τ = compute_τ(mesh)
# set up operators for p
params = [τ]
# set up matrix and affine component
#Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_wierd, bc_∇!, dbc_wierd)
Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_p, bc_∇!, dbc_p)
if neumann
    m,n = size(Δᵖ)
    nΔᵖ = spzeros(m+1,n+1)
    @. nΔᵖ[1:n, 1:m] = Δᵖ
    @. nΔᵖ[1:n,m+1] = 1.0
    @. nΔᵖ[n+1,1:m] = 1.0
    #nΔᵖ[1,m+1] = -1.0
    #nΔᵖ[n+1,1] = -1.0
    dropϵzeros!(nΔᵖ)
    maximum(abs.((nΔᵖ - nΔᵖ' ) ./ 2))
    nΔᵖ = (nΔᵖ + nΔᵖ' ) ./ 2
    dropϵzeros!(nΔᵖ)
    eigvals(Array(nΔᵖ))
    chol_Δᵖ = lu(-nΔᵖ)
else
    Δᵖ = (Δᵖ + Δᵖ' ) ./ 2
    dropϵzeros!(Δᵖ)
    chol_Δᵖ = cholesky(-Δᵖ)
end



# we are just going to check that we reproduce the correct solution from the pressure solve



println("=============================")

# check numerical incompressibility, (should be zero)
∇⨀𝐮 = similar(u_exact)
∇⨀!(∇⨀𝐮 , u_exact, v_exact, mesh)
max_div = maximum(abs.(∇⨀𝐮))
println("The maximum value of the numerically computed divergence is $(max_div)")

# check numerical advection
exact = u_exact .* ∂ˣu_exact + v_exact .* ∂ʸu_exact
u∂ˣu⨁v∂ʸu = similar(exact)
advec!(u∂ˣu⨁v∂ʸu , u_exact, v_exact, u_exact, mesh)
advection_error_u = rel_error(exact, u∂ˣu⨁v∂ʸu)
println("The relative error of the advection for u is $(advection_error_u )")

exact = u_exact .* ∂ˣv_exact + v_exact .* ∂ʸv_exact
u∂ˣv⨁v∂ʸv = similar(exact)
advec!(u∂ˣv⨁v∂ʸv , u_exact, v_exact, v_exact, mesh)
advection_error_v = rel_error(exact, u∂ˣv⨁v∂ʸv)
println("The relative error of the advection for v is $(advection_error_v )")

# the numerical error for the divergence of the nonlinear term is
exact = -Δp_exact
∇⨀ũ = similar(exact)
∇⨀!(∇⨀ũ , u∂ˣu⨁v∂ʸu, u∂ˣv⨁v∂ʸv, mesh)
Δ_error_p = rel_error(exact, ∇⨀ũ )
println("The relative error for the divergence of the nonlinear part is $(advection_error_v )")
# we now compute the symmetric advection term which should exactly cancel out the pressure gradient

#The numerical error of numerically computing the second derivative directly
exact = mesh.J .* ( mesh.M * ( Δp_exact  ) ) - bᵖ
Δp = reshape(Δᵖ * p_exact[:], size(exact))

Δ_error_p = rel_error(exact, Δp)
println("The relative error for the laplacian operator is $(Δ_error_p)")

println("The relative error in the pressure solve utilizing the nonlinear component is ")

if neumann
    m = length(bᵖ)+1;
    rhs_p = zeros(m)
    tmp = mesh.J .* ( mesh.M * ( ∇⨀ũ   ) ) + bᵖ
    @. rhs_p[1:(m-1)] = tmp[:]
    p̃ = reshape( (chol_Δᵖ \ rhs_p[:])[1:(m-1)], size(∇⨀ũ))
    gauge = sum(p_exact) / length(p_exact)
    @. p̃ += gauge
    error_p = rel_error(p_exact, p̃)
    println("The relative error for the pressure solve is $(error_p)")
    if timing
        @btime p̃ = reshape( (chol_Δᵖ \ rhs_p[:])[1:(m-1)], size(∇⨀ũ));
    end
elseif wierd
    rhs_p = mesh.J .* ( mesh.M * ( ∇⨀ũ   ) ) + bᵖ
    p̃ = reshape( (chol_Δᵖ \ rhs_p[:]) , size(∇⨀ũ) )
    gauge = sum(p_exact - p̃) / length(p_exact)
    @. p̃ += gauge
    error_p = rel_error(p_exact, p̃)
    println("The relative error for the pressure solve is $(error_p)")
    if timing
        @btime p̃ = reshape( (chol_Δᵖ \ rhs_p[:]) , size(∇⨀ũ) );
    end
else
    rhs_p = mesh.J .* ( mesh.M * ( ∇⨀ũ   ) ) + bᵖ
    p̃ = reshape( (chol_Δᵖ \ rhs_p[:]) , size(∇⨀ũ) )

    error_p = rel_error(p_exact, p̃)
    println("The relative error for the pressure solve is $(error_p)")
end



# now check to see if the pressure solve can eliminate the gradient of a potential
ϕ_analytic(x,y,t) = 1 / (2*π) * sin(2 * π * y ) * sin(2 * π * x)
ϕ_exact = eval_grid(ϕ_analytic, mesh, t)
# compute the gradient
# x-component
∂ˣϕ_analytic(x,y,t) = sin(2 * π * y ) * cos(2 * π * x)
∂ˣϕ_exact = eval_grid(∂ˣϕ_analytic, mesh, t)
# y-component
∂ʸϕ_analytic(x,y,t) = cos(2 * π * y ) * sin(2 * π * x)
∂ʸϕ_exact = eval_grid(∂ʸϕ_analytic, mesh, t)

# Δ
Δϕ_analytic(x,y,t) = - ( (2π) + (2π) ) * sin(2 * π * y ) * sin(2 * π * x)
Δϕ_exact = eval_grid(Δϕ_analytic, mesh, t)

# set new boundary conditions
# dirichlet
bc_p = (mesh.vmapB[1:end], mesh.mapB[1:end], ϕ_exact[mesh.vmapB[1:end]])
dbc_p = ([],[], 0.0, 0.0)
# neumann, with normal
if neumann
    bc_p = ( [], [], 0.0 )
    dbc_p = (mesh.vmapB[1:end], mesh.mapB[1:end], ∂ˣϕ_exact[mesh.vmapB[1:end]], ∂ʸϕ_exact[mesh.vmapB[1:end]])
elseif wierd
    bc_p =  ([mesh.vmapB[1]], [mesh.mapB[1]], 0.0 )
    dbc_p = (mesh.vmapB[2:end], mesh.mapB[2:end], ∂ˣϕ_exact[mesh.vmapB[2:end]], ∂ʸϕ_exact[mesh.vmapB[2:end]])
end

Δᵖ, bᵖ = poisson_setup_bc(field, params, mesh, bc!, bc_p, bc_∇!, dbc_p)
if neumann
    m,n = size(Δᵖ)
    nΔᵖ = spzeros(m+1,n+1)
    @. nΔᵖ[1:n, 1:m] = Δᵖ
    @. nΔᵖ[1:n,m+1] = 1.0
    @. nΔᵖ[n+1,1:m] = 1.0
    #nΔᵖ[1,m+1] = -1.0
    #nΔᵖ[n+1,1] = -1.0
    dropϵzeros!(nΔᵖ)
    maximum(abs.((nΔᵖ - nΔᵖ' ) ./ 2))
    nΔᵖ = (nΔᵖ + nΔᵖ' ) ./ 2
    dropϵzeros!(nΔᵖ)
    eigvals(Array(nΔᵖ))
    chol_Δᵖ = lu(-nΔᵖ)
else
    Δᵖ = (Δᵖ + Δᵖ' ) ./ 2
    dropϵzeros!(Δᵖ)
    chol_Δᵖ = cholesky(-Δᵖ)
end

# now compute the divergence of the test solution
fx = u_exact - ∂ˣϕ_exact
fy = v_exact - ∂ʸϕ_exact

∇⨀ũ = similar(fy)
∇⨀!(∇⨀ũ, fx, fy, mesh)
@. ∇⨀ũ *= 1.0
Δ_error = rel_error(-Δϕ_exact, ∇⨀ũ)


# now solve
println("Now we test the project part of the operator a little differently")
println("The error in computing the second derivative is $(Δ_error )")
if neumann
    m = length(bᵖ)+1;
    rhs_p = zeros(m)
    tmp = mesh.J .* ( mesh.M * ( ∇⨀ũ   ) ) + bᵖ
    @. rhs_p[1:(m-1)] = tmp[:]
    p̃ = reshape( (chol_Δᵖ \ rhs_p[:])[1:(m-1)], size(∇⨀ũ))
    gauge = sum(p_exact) / length(p_exact)
    @. p̃ += gauge
    error_p = rel_error(ϕ_exact, p̃)
    println("The relative error for the pressure solve is $(error_p)")
    if timing
        @btime p̃ = reshape( (chol_Δᵖ \ rhs_p[:])[1:(m-1)], size(∇⨀ũ));
    end
elseif wierd
    rhs_p = mesh.J .* ( mesh.M * ( ∇⨀ũ   ) ) + bᵖ
    p̃ = reshape( (chol_Δᵖ \ rhs_p[:]) , size(∇⨀ũ) )
    gauge = sum(p_exact - p̃) / length(p_exact)
    @. p̃ += gauge
    error_p = rel_error(ϕ_exact, p̃)
    println("The relative error for the pressure solve is $(error_p)")
    if timing
        @btime p̃ = reshape( (chol_Δᵖ \ rhs_p[:]) , size(∇⨀ũ) );
    end
else
    rhs_p = mesh.J .* ( mesh.M * ( ∇⨀ũ   ) ) + bᵖ
    p̃ = reshape( (chol_Δᵖ \ rhs_p[:]) , size(∇⨀ũ) )

    error_p = rel_error(ϕ_exact, p̃)
    println("The relative error for the pressure solve using Dirichlet bc is $(error_p)")
end

# check if new field is incompressible
∂ˣp̃ = similar(p̃)
∂ʸp̃ = similar(p̃)
∇!(∂ˣp̃, ∂ʸp̃,  p̃, mesh)

#include lift terms for DG consistency
xlift = @. mesh.nx[:] * (p̃[mesh.vmapP]-p̃[mesh.vmapM])
ylift = @. mesh.ny[:] * (p̃[mesh.vmapP]-p̃[mesh.vmapM])

xlift = reshape(xlift, mesh.nFaces * mesh.nfp, mesh.K)
ylift = reshape(ylift, mesh.nFaces * mesh.nfp, mesh.K)

fx = u_exact - ∂ˣϕ_exact + ∂ˣp̃ + mesh.lift * ( mesh.fscale .* xlift) * 0.5
fy = v_exact - ∂ʸϕ_exact + ∂ʸp̃ + mesh.lift * ( mesh.fscale .* ylift) * 0.5

∇⨀!(∇⨀ũ , fx, fy, mesh)
incomp = maximum(abs.(∇⨀ũ ))
println("The maximum incompressibility of the new solution is $(incomp)")

xlift = @. mesh.nx[:] * (fx[mesh.vmapP]-fx[mesh.vmapM])
ylift = @. mesh.ny[:] * (fy[mesh.vmapP]-fy[mesh.vmapM])
xlift = reshape(xlift, mesh.nFaces * mesh.nfp, mesh.K)
ylift = reshape(ylift, mesh.nFaces * mesh.nfp, mesh.K)

w∇⨀ũ = ∇⨀ũ + mesh.lift * ( mesh.fscale .* xlift) * 0.5 + mesh.lift * ( mesh.fscale .* ylift) * 0.5
incomp = maximum(abs.(w∇⨀ũ ))
println("The maximum incompressibility (in the weak sense) of the new solution is $(incomp)")


relx = rel_error(∂ˣϕ_exact , ∂ˣp̃ + mesh.lift * ( mesh.fscale .* xlift) * 0.5 )
println("The error in the x-derivative is $(relx)")

rely = rel_error(∂ʸϕ_exact , ∂ʸp̃ + mesh.lift * ( mesh.fscale .* ylift) * 0.5)
println("The error in the y-derivative is $(rely)")

discontinuity_error = maximum(abs.(p̃[mesh.vmapP] - p̃[mesh.vmapM]))
println("The maximum discontinuity in p is $(discontinuity_error)")

discontinuity_error = maximum(abs.(∂ˣp̃[mesh.vmapP] - ∂ˣp̃[mesh.vmapM]))
println("The maximum discontinuity in px is $(discontinuity_error)")

discontinuity_error = maximum(abs.(∂ʸp̃[mesh.vmapP] - ∂ʸp̃[mesh.vmapM]))
println("The maximum discontinuity in py is $(discontinuity_error)")

∂ˣq = ∂ˣp̃ + mesh.lift * ( mesh.fscale .* xlift) * 0.5
∂ʸq = ∂ʸp̃ + mesh.lift * ( mesh.fscale .* ylift) * 0.5

discontinuity_error = maximum(abs.(∂ˣq[mesh.vmapP] - ∂ˣq[mesh.vmapM]))
println("The maximum discontinuity in qx is $(discontinuity_error)")

discontinuity_error = maximum(abs.(∂ʸq[mesh.vmapP] - ∂ʸq[mesh.vmapM]))
println("The maximum discontinuity in qy is $(discontinuity_error)")

#=
thing = abs.(f)
thing = log.(abs.(- ∂ˣϕ_exact + ∂ˣp̃))
thing = log.(abs.( - ∂ʸϕ_exact + ∂ʸp̃ ))
p3 = surface(mesh.x[:],mesh.y[:], thing , camera = (0,90))
display(p3)
=#


###
