using BandedMatrices
using BenchmarkTools

include("mesh2D.jl")
include("dg_advection.jl")
include("../DG2D/triangles.jl")
include("../DG2D/dg_poisson.jl")


# simulation parameters and grid
n = 1
FileName = "Maxwell1.neu"
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
#dirichlet
function bc_u!(du, u, bc)
    @. du[bc[2]] = 2 * u[bc[1]]
end
#neumann
function bc_φ!(fˣ, fʸ,φˣ, φʸ, bc)
    @. fˣ[bc[2]] = 2 * φˣ[bc[1]]
    @. fʸ[bc[2]] = 2 * φʸ[bc[1]]
end
# define boundary conditions
# check that it doesn't crash
Δu = similar(field.u)
u = similar(field.u)
#dg_poisson!(Δu, u, field, params, mesh, bc_u!, bc, bc_φ!, dbc)


∇² = poisson_setup(field, params, mesh, bc_u!, bc, bc_φ!, dbc)
s∇² = (∇² + ∇²')/2
