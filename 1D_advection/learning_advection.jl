# ∂_t θ + ∂_x (u θ) = 0
#we define different discretizations for the rhs
using Plots
using FFTW
using LinearAlgebra
using DifferentialEquations
using BenchmarkTools
using Profile
using Revise

#number of discretization points
n = 2^6
#domain length
L = 2 * pi
#collocation points
x = collect(0:(n-1))/n * L
#time-step size
dt = 0.0001*pi
#CFL
println("The CFL number is")
println((x[2]-x[1])/dt)
#define the rhs!
#include("1D_spectral.jl")
include("1D_central.jl")
#include("1D_fourth_order.jl")
#include("1D_sixth_order.jl")
#include("1D_upwind.jl")

t = 1.0;
θ₀ = exp.(-4*(x.-L/2).^2) #.+ 0 * im  #needs this for spectral
θ₁ = similar(θ₀) #just for consistencey
tspan = (0.0,2*pi)
u = [1.0 ] #parameters, velocity speed
#@btime rhs!(θ₁, θ₀, u, t)
prob = ODEProblem(rhs!, θ₀, tspan, u)

#AB3(), RK4(), Tsit5(), Heun()
@time sol = solve(prob, RK4(), dt=dt, adaptive = false)

#plotting
plt = plot(x,real.(sol.u[1]),xlims=(0,L), ylims = (0,1))
nt = length(sol.t)
ind = Int(floor(nt/20)) * collect(1:20)
ind[end] = length(sol.t)
for i in 1:20
    plt = plot(x,real.(sol.u[ind[i]]),xlims=(0,L), ylims = (0,1), label = "evolution", marker = 3)
    plot!(x,real.(sol.u[1]),xlims=(0,L), ylims = (0,1), color= "red", label = "starting")
    display(plt)
    #gui()
end

println("The error in the solution is ")
wrongness = norm(sol.u[1]-sol.u[end]) / norm(sol.u[1])
println(wrongness)
