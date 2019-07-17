# first define the stream function
include("grid2D.jl")
include("dg_advection.jl")
# choose the polynomial order
#3 seems to be pretty efficient
n = 1
timings = false
gradients_check = false
solve_ode = true
euler = false
upwind_check = false
#load file
FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
grid = garbage_triangle3(n, filename)
field = dg_garbage_triangle(grid)

#for convenience
x = grid.x
y = grid.y
#plot the total grid points
p1 = scatter(x,y,legend=false)
# plot boundary of triangles
scatter!(x[grid.vmapM] , y[grid.vmapM], color = "black", legend = false)
#plot boundary of domain
scatter!(x[grid.vmapB] , y[grid.vmapB], color = "yellow", legend = false)
display(plot(p1))

println("We have")
println(length(grid.x))
println("degrees of freedom")
#define stream function and components of velocity
ψ(x, y, γ) = exp(γ*(y-1)^2 ) * cos(π/2 * x) * cos(π/2 * y)
u1(x, y, γ) =  cos(π/2 * y) * cos(π/2 * x) * γ * 2 * (y-1) * exp(γ*(y-1)^2 )  - π / 2 * sin(π/2 * y) * exp(γ*(y-1)^2 ) * cos(π/2 * x)
u2(x, y, γ) = π / 2 * sin(π/2 * x) * exp(γ*(y-1)^2 ) * cos(π/2 * y)
u0(x, y, μ) = exp(-μ * x^2 - μ * (y+0.5)^2) * cos(π/2 * x) * cos(π/2 * y)

#simpler
#=
ψ(x, y, γ)  = y^2 + x
u1(x, y, γ) =  2 * y
u2(x, y, γ) = -1
=#
#u0(x, y, μ) = sin(x)*cos(y) + x
#u0(x, y, μ) =  1.0



#define initial conditions and velocity field
γ = -0.0
μ = 10.0
u⁰ = [u0(x[i,j],y[i,j],μ) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]
ψᵏ = [ψ(x[i,j],y[i,j],γ) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]
v¹ = [u1(x[i,j],y[i,j],γ) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]
v² = [u2(x[i,j],y[i,j],γ) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]

flux1 = v¹ .* u⁰
flux2 = v² .* u⁰

struct velocity_field{T}
    v1::T
    v2::T
    function velocity_field(w1, w2)
        v1 = w1
        v2 = w2
        return new{typeof(w1)}(v1, v2)
    end
end


#define params
tspan = (0.0, 8.0)
ι = field
ε = external
𝒢 = grid
#rhs! = dg_central_2D!
#rhs! = dg_rusonov_2D!
rhs! = dg_central_2D!
dt =  0.5 * (grid.r[2] - grid.r[1]) / grid.K / maximum([1, maximum(v¹)])
println("The time step size is ")
println(dt)
# find numerical velocity field
∇!(ι.φˣ, ι.φʸ, ψᵏ, 𝒢)
w¹ = copy(ι.φʸ)
w² = -copy(ι.φˣ)
external = velocity_field(v¹, v²) #use numerical instead of exact derivative



params = (grid, field, external)

@. field.u = u⁰
u = field.u
u̇ = field.u̇

if timings
    println("central")
    @btime dg_central_2D!(u̇, u, params, 0);
    println("upwind")
    @btime dg_upwind_2D!(u̇, u, params, 0);
    println("divergence")
    @btime ∇⨀!(u̇, ι.φˣ, ι.φʸ, 𝒢);
    println("lift")
    @btime lift = 𝒢.lift * (𝒢.fscale .* ι.fⁿ);
end

if gradients_check
    #check incompressibility 1
    ∇⨀!(u̇,v¹, v², 𝒢)
    println("The infinity norm of the divergence of velocity field 1 is ")
    infnorm = maximum(abs.(u̇))
    println(infnorm)
    #check incompressibility 2
    ∇⨀!(u̇, w¹, w², 𝒢)
    println("The infinity norm of the divergence of velocity field 2 is ")
    infnorm = maximum(abs.(u̇))
    println(infnorm)
    #check gradients
    ∇!(ι.φˣ, ι.φʸ, ψᵏ, 𝒢)
    println("Checking stream function vs velocity field")
    println("The infinity norm of the first component of velocity field is ")
    infnorm = maximum(abs.( ι.φʸ - v¹))
    println(infnorm)
    println("The infinity norm of the second component of velocity field is ")
    infnorm = maximum(abs.(  ι.φˣ + v²))
    println(infnorm)
    println("Checking v dot grad theta vs grad (u theta)")
    flux1 = v¹ .* u⁰
    flux2 = v² .* u⁰
    ∇⨀!(u̇,flux1, flux2, 𝒢)
    ∇!(ι.φˣ, ι.φʸ, u⁰, 𝒢)
    @. u = v¹ * ι.φˣ + v² * ι.φʸ
    println("The infinity norm commutation error is ")
    infnorm = maximum(abs.(  u̇ - u))
    println(infnorm)
    println("Checking w dot grad theta vs grad (u theta)")
    flux1 = w¹ .* u⁰
    flux2 = w² .* u⁰
    ∇⨀!(u̇,flux1, flux2, 𝒢)
    ∇!(ι.φˣ, ι.φʸ, u⁰, 𝒢)
    @. u =  w¹ * ι.φˣ + w² * ι.φʸ
    println("The infinity norm commutation error is ")
    infnorm = maximum(abs.(  u̇ - u))
    println(infnorm)
end


@. field.u = u⁰
u = copy(field.u)
u̇ = copy(field.u̇)

if solve_ode
    prob = ODEProblem(rhs!, u, tspan, params);
    sol  = solve(prob, RK4(), dt=dt, adaptive = false); # AB3(), RK4(), Tsit5()
end

#euler time-stepping for debugging


if euler
    @. field.u = u⁰
    u = copy(field.u)
    u̇ = copy(field.u̇)
    for i in 1:1
        dg_central_2D!(ι.u̇, ι.u, params, 0)
        @. ι.u += dt * ι.u̇
    end
end

if upwind_check
    @. field.u = u⁰
    u = copy(field.u)
    u̇ = copy(field.u̇)
    for i in 1:20
        dg_central_2D!(ι.u̇, ι.u, params, 0)
        @. ι.u += dt * ι.u̇
    end
    u = copy(ι.u)
    flux1 = v¹ .* u
    flux2 = v² .* u
    fxP  = flux1[𝒢.vmapP]
    fyP  = flux2[𝒢.vmapP]
    fxM  = flux1[𝒢.vmapM]
    fyM  = flux2[𝒢.vmapM]
    fnP  = @. fxP * 𝒢.nx[:] + fyP * 𝒢.ny[:]
    fnM  = @. fxM * 𝒢.nx[:] + fyM * 𝒢.ny[:]
    fnP = reshape( fnP, size(ι.fˣ) )
    fnM = reshape( fnM, size(ι.fˣ) )
    vn = reshape( ε.v1[𝒢.vmapM] .* 𝒢.nx[:] + ε.v2[𝒢.vmapM] .* 𝒢.ny[:] , size(ι.fˣ) )
    #now for the normal component along the faces, with upwind
    ujump = reshape( abs.(ε.v1[𝒢.vmapM] .* 𝒢.nx[:] + ε.v2[𝒢.vmapM] .* 𝒢.ny[:]) .* (u[𝒢.vmapM] - u[𝒢.vmapP]), size(ι.fˣ) )
    index1 = 5
    index2 = 25
    println("----------------------")
    println("the normal component of velocity on element 1")
    display(vn[index1,index2])
    println("the flux on the interior is")
    display(fnM[index1,index2])
    println("the flux on the exterior is")
    display(fnP[index1,index2])
    println("the chosen upwind flux is")
    upwind = @. (fnP + fnM)/2 - 0.5 * ujump
    display(upwind[index1,index2])
    whichflux = vn[index1,index2]<0 ? "interior" : "exterior"
    whichfluxval = vn[index1,index2]<0 ? fnM[index1,index2] : fnP[index1,index2]
    println("This should have chosen the one on the "*whichflux)
    println(whichfluxval)
    println(whichfluxval ≈ upwind[index1,index2])
    println("----------------------")
end


# [ max(x[i,j], y[i,j]) for i in 1:length(x[:,1]), j in 1:length(y[1,:]) ]


###
 p1 = scatter(x[:,1],y[:,1])
for i in 2:10
 scatter!(x[:,i],y[:,i], legend=false)
end
display(plot(p1))
###

###
gr()
endtime = length(sol.t)
steps = Int( floor(endtime / 80))

@gif for i in 1:steps:endtime
    println(i)
    u = copy(sol.u[i])
    p1 = surface(x[:],y[:],u[:], camera = (15,60))
    plot(p1)
end
###

###
gr()

u = copy(sol.u[1])
p1 = surface(x[:],y[:],u[:], camera = (15,60))
plot(p1)

###
