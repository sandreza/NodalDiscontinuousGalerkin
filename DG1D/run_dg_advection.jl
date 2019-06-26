include("dg1D.jl")
include("dg_advection.jl")

using Plots
using BenchmarkTools
using DifferentialEquations

K = 2^4 #number of elements
n = 2^2-1 #polynomial order,
# (for 2^8, 2^4 with 2^4-1 seems best)


println("The degrees of freedom are ")
println( (n+1)*K)
#domain parameters
L = 2π
xmin = 0.0
xmax = L
#speed of wave
v = 2π
α = 0.0 #1 is central flux, 0 is upwind

ι = dg(K, n, xmin, xmax)
struct extern_params_parametric{T,S}
    v::T
    α::S
end
ε = extern_params_parametric(v, α) # first is velocity, second is value for α
params = (ι, ε)
x = ι.x
u = ι.u
Δx =  minimum(x[2,:] - x[1,:])
CFL = 0.75
dt = CFL * Δx / v
dt *= 0.5 / 1
@. u = sin(ι.x) #initial condition for other problem
@. u = exp.(-4*(ι.x .- L/2).^2) #initial condition for periodic problem
rhsu = ι.rhsu
#=
t = 0
dg_upwind!(rhsu, u, params,
dt = 0.1
@. u += dt * rhsu
t += dt
dg_upwind!(rhsu, u, params, t)
=#
make_periodic1D!(ι.vmapP, ι.u)
rhs! = dg_upwind_p!
#run code
tspan = (0.0,2)

prob = ODEProblem(rhs!, u, tspan, params);
#AB3(), RK4(), Tsit5(), Heun()
sol = solve(prob, Tsit5(), dt=dt, adaptive = false);
# @code_warntype dg_upwind!(ι.rhsu, ι.u, params, 0)
# @btime dg_upwind!(ι.rhsu, ι.u, params, 0)
#@btime sol = solve(prob, Tsit5(), dt=dt, adaptive = false);


#plotting
theme(:juno)
L = 2π
plt = plot(x,real.(sol.u[1]), xlims=(0,L), ylims = (-1.1,1.1))
nt = length(sol.t)
num = 20
ind = Int(floor(nt/num)) * collect(1:num)
ind[end] = length(sol.t)
for i in 1:num
    plt = plot(x, sol.u[ind[i]] ,xlims=(0,L), ylims = (-1.1,1.1), leg = false, marker = 3)
    plot!(x, sol.u[1] ,xlims=(0,L), ylims = (-1.1,1.1), color= "red", leg = false)
    display(plt)
    sleep(0.25)
end
#get_color_palette(:auto, plot_color(:white), 18)
#can set color = 1, to an index
plt = plot(x, sol.u[ind[end]] ,xlims=(0,L), ylims = (-1.1,1.1), leg = false, marker = 3)
plot!(x, sol.u[1] ,xlims=(0,L), ylims = (-1.1,1.1), color= "red", leg = false)
display(plt)

println("The error in the solution is ")
wrongness = norm(sol.u[1]-sol.u[end]) / norm(sol.u[1])
println(wrongness)
println("To find where the computational bottleneck is")
println("Evaluating the right hand side takes")
@btime dg_upwind!(ι.rhsu, ι.u, params, 0)
println("Performing a matrix multiplication")
@btime mul!(ι.rhsu, ι.D, ι.u)
