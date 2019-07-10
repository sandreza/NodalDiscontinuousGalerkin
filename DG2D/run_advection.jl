# first define the stream function

# choose the polynomial order
#3 seems to be pretty efficient
n = 3
timings = false
#load file
FileName = "Maxwell025.neu"
filepath = "./DG2D/grids/"
filename = filepath*FileName
grid = garbage_triangle3(n, filename)
field = dg_garbage_triangle(grid)

#plot the grid
p1 = scatter(grid.x,grid.y,legend=false)
display(plot(p1))
println("We have")
println(length(grid.x))
println("degrees of freedom")
#define stream function and components of velocity
ψ(x, y, γ) = exp(γ*(y-1)^2 ) * cos(π/2 * x) * cos(π/2 * y)
u1(x, y, γ) =  cos(π/2 * y) * cos(π/2 * x) * γ * 2 * y * exp(γ*(y-1)^2 )  - π / 2 * sin(π/2 * y) * exp(γ*(y-1)^2 ) * cos(π/2 * y)
u2(x, y, γ) = π / 2 * sin(π/2 * x) * exp(γ*(y-1)^2 ) * cos(π/2 * x) * cos(π/2 * y)

u0(x, y, μ) = exp(-μ * x^2 - μ * (y+0.5)^2)

#define initial conditions and velocity field
γ = -0.0
μ = 10.0
x = grid.x
y = grid.y
u⁰ = [u0(x[i,j],y[i,j],μ) for i in 1:length(x[:,1]), j in 1:length(y[1,:])]
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

external = velocity_field(v¹, v²)

#define params
tspace = (0, 2)
ι = field
ε = external
𝒢 = grid
rhs! = dg_central_2D!
dt =  (grid.r[2] - grid.r[1]) / grid.K / maximum([1, maximum(v¹)])

params = (grid, field, external)

@. field.u = u⁰
u = field.u
u̇ = field.u̇
if timings
    @btime dg_central_2D!(u̇, u, params, 0);
    @btime ∇⨀!(u̇, ι.φˣ, ι.φʸ, 𝒢);
    @btime lift = 𝒢.lift * (𝒢.fscale .* ι.fⁿ);
end

prob = ODEProblem(rhs!, u, tspan, params);
sol  = solve(prob, Tsit5(), dt=dt, adaptive = false); # AB3(), RK4(),

#=
#euler time-stepping for debugging
for i in 1:0
    dg_central_2D!(ι.u̇, ι.u, params, 0)
    @. ι.u += dt * ι.u̇
end
=#
