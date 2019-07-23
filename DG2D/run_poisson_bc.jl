# run_poisson_bc

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
plotting_matrix = true
check_correctness = true
plotting_solution = false
# simulation parameters and grid
n = 4
FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
mesh = garbage_triangle3(n, filename)
field = dg_garbage_triangle(mesh)
ι = field

# location of boundary grid points for dirichlet bc
bc = (mesh.vmapB, mesh.mapB)
bc = ([],[])
# location of boundary grid points for neumann bc
dbc = ([],[])
dbc = (mesh.vmapB, mesh.mapB)

#compute tau
τ = compute_τ(mesh)
params = [τ] 
#for the first poisson equation
#dirichlet
const shift = 00.0
function bc_u!(ι, mesh, bc)
    @. ι.fⁿ[bc[2]] = 1 * ι.u[bc[1]] - shift
    return nothing
end
#neumann
function bc_φ!(ι, mesh, bc)
    @. ι.fˣ[bc[2]] = ι.φˣ[bc[1]] - 0.0 - 2 * mesh.x[bc[1]]
    @. ι.fʸ[bc[2]] = ι.φʸ[bc[1]] - 0.0 - 2 * mesh.y[bc[1]]
    return nothing
end
# define boundary conditions
# check that it doesn't crash
Δu = similar(field.u)
u = similar(field.u)
#dg_poisson!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)

#( (mesh.x[bc[1]])^2 + (mesh.y[bc[1]])^2)*1

# may take a while for larger matrices
∇², b = poisson_setup_bc(field, params, mesh, bc_u!, bc, bc_φ!, dbc)
# make sure its numericall symmetric
symmetric_check = sum(abs.(∇² .- (∇² + ∇²')./2)) / length(∇²) / maximum(abs.(∇²))
if symmetric_check > eps(1.0)
    println("warning the matrix is not numerically symmetric")
    ∇² = (∇² + ∇²')/2
else
    ∇² = (∇² + ∇²')/2
end
i,j = findnz(∇²)
#manually drop the zeros
println("current number of nonzero entries is $(length(i))")
drop_criteria = eps(maximum(abs.(∇²)))
for loop in 1:length(i)
    if abs(∇²[i[loop],j[loop]]) < drop_criteria
        #println("dropping")
        #println(abs(∇²[i[loop],j[loop]]))
        #println((i[loop],j[loop]))
        ∇²[i[loop],j[loop]] = 0.0
    end
end
dropzeros!(∇²)
i,j = findnz(∇²)
println("now it is $(length(i))")

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

    exact(x,y,α,β) = cos(π/2 * x * α) * cos(π/2 * y * β) + shift
    #exact(x,y,α,β) = (x^2 -1 ) * (y^2 -1 )
    exact(x,y,α,β) = x^2 + y^2
    # then create a forcing function

    #forcing(x,y,α,β) = - ( (α*π/2)^2 + (β*π/2)^2 ) * cos(π/2 * x * α) * cos(π/2 * y * β)
    #forcing(x,y,α,β) = 2.0 * (y^2 - 1.0) + 2.0 * (x^2 - 1.0)
    forcing(x,y,α,β) = 4.0

    #for convenience
    x = mesh.x
    y = mesh.y

    # evaluate at grid points with given values for α and β
    # odd for dirichlet, even for neumann
    α = 2
    β = 2
    frhs = [forcing(x[i,j],y[i,j],α,β) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]
    #@. frhs[mesh.vmapB] = (x[mesh.vmapB])^2 + (y[mesh.vmapB])^2

    # adjust for J * mass matrix component
    frhs = mesh.J .* (mesh.M * frhs) - b
    frhs_with_affine = mesh.J .* (mesh.M * frhs)

    fsol = [exact(x[i,j],y[i,j],α,β) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]

    # chech, W^{2,∞} error
    println("----------------")
    @. u = fsol
    e1 = similar(u)
    @. e1 = 0.0
    e1[2] = 0.0
    # @. u = e1
    dg_poisson_bc!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)
    w2inf_Δ = maximum(abs.(Δu .- frhs_with_affine)) / maximum(abs.(frhs))
    cc = ∇² * u[:]
    ctmp = frhs[:]
    w2inf_Δ_2 = maximum(abs.( cc .- ctmp)) / maximum(abs.(frhs))
    println("The relative error in computing the second derivative is $(w2inf_Δ)")
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
    @btime dg_poisson_bc!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)

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

if plotting_matrix
    p1 = spy(∇²)
    p2 = spy(cm∇²)
    display(plot(p1,p2))
end

if plotting_solution
    gr()
    camera_top = 90 #this is a very hacky way to get a 2D contour plot
    camera_side = 0
    p1 = surface(x[:],y[:],u[:], camera = (camera_side,camera_top))
    p2 = surface(x[:],y[:],fsol[:], camera = (camera_side,camera_top))
    p3 = surface(x[:],y[:],fsol[:].-u[:], camera = (camera_side,camera_top))

    display(plot(p3))
end
