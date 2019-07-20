using BandedMatrices
using BenchmarkTools
using LinearAlgebra
using Plots

include("mesh2D.jl")
include("dg_advection.jl")
include("../DG2D/triangles.jl")
include("../DG2D/dg_poisson.jl")
include("../src/CuthillMckee.jl")


timings = false
plotting = true
check_correctness = true
# simulation parameters and grid
n = 6
FileName = "Maxwell0125.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)

# location of boundary grid points for dirichlet bc
bc = (mesh.vmapB, mesh.mapB)
# location of boundary grid points for neumann bc
dbc = ([],[])

#compute tau
τ = compute_τ(mesh)
params = [τ]
#homogenous dirichlet
function bc_u!(du, u, bc)
    @. du[bc[2]] = 2 * u[bc[1]]
    return nothing
end
#homogenous neumann
function bc_φ!(fˣ, fʸ, φˣ, φʸ, bc)
    @. fˣ[bc[2]] = 2 * φˣ[bc[1]]
    @. fʸ[bc[2]] = 2 * φʸ[bc[1]]
    return nothing
end
# define boundary conditions
# check that it doesn't crash
Δu = similar(field.u)
u = similar(field.u)
#dg_poisson!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)


# may take a while for larger matrices
∇² = poisson_setup(field, params, mesh, bc_u!, bc, bc_φ!, dbc)
# make sure its symmetric
∇² = (∇² + ∇²')/2

# output some matrix properties
println("The size of the matrix is $(size(∇²))")
i,j = findnz(∇²)
println("The bandwidth of the matrix is $(maximum(i-j)+1)")
println("The sparsity is $(length(nonzeros(∇²)) / length(∇²))")

p = symrcm(∇²)
cm∇² = sparse(∇²[p,p])
i,j = findnz(cm∇²)
println("The bandwidth of the reordered matrix is $(maximum(i-j)+1)")

if check_correctness
    # first create an exact solution
    exact(x,y,α,β) = cos(π/2 * x * α) * cos(π/2 * y * β)

    # then create a forcing function
    forcing(x,y,α,β) = - ( (α*π/2)^2 + (β*π/2)^2 ) * cos(π/2 * x * α) * cos(π/2 * y * β)

    #for convenience
    x = mesh.x
    y = mesh.y

    # evaluate at grid points with given values for α and β
    # odd for dirichlet, even for neumann
    α = 3
    β = 3
    frhs = [forcing(x[i,j],y[i,j],α,β) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]

    # adjust for J * mass matrix component
    frhs = mesh.J .* (mesh.M * frhs)

    fsol = [exact(x[i,j],y[i,j],α,β) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]

    # chech, W^{2,∞} error
    println("----------------")
    @. u = fsol
    dg_poisson!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)
    w2inf = maximum(abs.(Δu .- frhs)) / maximum(abs.(frhs))
    println("The relative error in computing the second derivative is $(w2inf)")
    println("This is a lower estimate since its on the grid points")

    # now to compute the solution
    #chol_∇² = cholesky(-∇²); #will need to multiply by -1
    chol_∇² = lu(-∇²); #will need to multiply by -1
    @. Δu = - frhs #due to cholesky nonsense
    tmpΔu = Δu[:]
    tmpu = u[:]
    #ldiv!(tmpΔu, chol_∇², tmpu)
    tmpu = chol_∇² \ tmpΔu #just using the fastest
    #modify for neumann
    tmpu = tmpu .- sum(tmpu)/length(tmpu) .+ sum(fsol)/length(fsol)
    #set values
    @. u[:] = tmpu
    w2inf = maximum(abs.(u .- fsol)) / maximum(abs.(u))
    println("The relative error in computing the solution is $(w2inf)")
    println("----------------")

    println("-------------")
    jump_max = maximum(abs.(u[mesh.vmapP] .- u[mesh.vmapM]))
    println("The maximum discontinuity across gridpoints is $(jump_max)")
    println("-------------")

end

if timings
    # create full matrix
    f∇² = Symmetric(Array(∇²));
    # create banded matrix
    mat_size = size(cm∇²)
    i,j = findnz(cm∇²)
    band = maximum(i-j) + 1
    b∇² = BandedMatrix(zeros(mat_size), (band,band))
    @. b∇² = cm∇²

    #for comparison
    println("------------")
    println("evaluating the second derivative takes")
    @btime dg_poisson!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)

    # now for timings,
    #these are somewhat irrelevant hence the commenting
    #=
    println("--------------")
    println("sparse")
    @btime ∇² \ u[:];
    println("full")
    @btime f∇² \ u[:];
    println("reordered")
    @btime cm∇² \ u[p];
    println("banded")
    @btime b∇² \ u[p];
    =#
    println("--------------")

    chol_f∇² = cholesky(-f∇²); #will need to multiply by -1
    chol_cm∇² = cholesky(-cm∇²); #will need to multiply by -1
    chol_∇² = cholesky(-∇²); #will need to multiply by -1

    println("cholesky sparse")
    @btime chol_∇² \ u[:];

    println("cholesky full")
    @btime chol_f∇² \ u[:];

    println("cholesky reordered")
    @btime chol_cm∇² \ u[p];

    println("cholesky banded is not an option")

    println("--------------")
    lu_f∇² = lu(f∇²);
    lu_∇² = lu(∇²) ;
    lu_cm∇² = lu(cm∇²) ;
    qr_b∇² = qr(b∇²) ;

    println("lu sparse")
    @btime lu_∇² \ u[:];

    println("lu full")
    @btime lu_f∇² \ u[:];

    println("lu reordered")
    @btime lu_cm∇² \ u[p];

    println("qr banded")
    @btime qr_b∇² \ u[p];
end

if plotting
    p1 = spy(∇²)
    p2 = spy(cm∇²)
    display(plot(p1,p2))
end
