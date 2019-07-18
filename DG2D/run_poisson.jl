using BandedMatrices
using BenchmarkTools
using LinearAlgebra
using Plots

include("mesh2D.jl")
include("dg_advection.jl")
include("../DG2D/triangles.jl")
include("../DG2D/dg_poisson.jl")
include("../src/CuthillMckee.jl")


timings = true
plotting = true
check_correctness = true
# simulation parameters and grid
n = 10
FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)

bc = (mesh.vmapB, mesh.mapB)
dbc = ([],[])
#compute tau
function compute_τ(mesh)
    matP = mesh.J[mesh.vmapP] ./ mesh.sJ[:]
    matM = mesh.J[mesh.vmapM] ./ mesh.sJ[:]
    hmin = zeros(length(matP))
    for i in 1:length(matP)
        matP[i] < matM[i] ? hmin[i] = 2 * matP[i] : hmin[i] = 2* matM[i]
    end
    np = (mesh.n + 1) * (mesh.n + 2) / 2
    return reshape(np ./ hmin, mesh.nfp * mesh.nFaces, mesh.K)
end
τ = compute_τ(mesh)
params = [τ]
#homogenous dirichlet
function bc_u!(du, u, bc)
    @. du[bc[2]] = 2 * u[bc[1]]
end
#homogenous neumann
function bc_φ!(fˣ, fʸ, φˣ, φʸ, bc)
    @. fˣ[bc[2]] = 2 * φˣ[bc[1]]
    @. fʸ[bc[2]] = 2 * φʸ[bc[1]]
end
# define boundary conditions
# check that it doesn't crash
Δu = similar(field.u)
u = similar(field.u)
#dg_poisson!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)


∇² = poisson_setup(field, params, mesh, bc_u!, bc, bc_φ!, dbc)
∇² = (∇² + ∇²')/2



println("The size of the matrix is $(size(∇²))")
i,j = findnz(∇²)
println("The bandwidth of the matrix is $(maximum(i-j)+1)")
println("The sparsity is $(length(nonzeros(∇²)) / length(∇²))")

p = symrcm(∇²)
cm∇² = sparse(∇²[p,p])
i,j = findnz(cm∇²)
println("The bandwidth of the reordered matrix is $(maximum(i-j)+1)")


if timings
    # create full matrix
    f∇² = Symmetric(Array(∇²));
    # create banded matrix
    mat_size = size(cm∇²)
    i,j = findnz(cm∇²)
    band = (maximum(i-j)+1)
    b∇² = BandedMatrix(zeros(mat_size), (band,band))
    @. b∇² = cm∇²

    #for comparison
    println("------------")
    println("evaluating the second derivative takes")
    @btime dg_poisson!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)

    # now for timings
    println("--------------")
    println("sparse")
    @btime ∇² \ u[:];
    println("full")
    @btime f∇² \ u[:];
    println("reordered")
    @btime cm∇² \ u[p];
    println("banded")
    @btime b∇² \ u[p];

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
