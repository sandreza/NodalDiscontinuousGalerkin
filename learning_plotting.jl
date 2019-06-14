using Plots
using Revise
using BenchmarkTools

include("learning_functions.jl")
mutable struct exp_down
    dt; x; f; t;
end
decay = exp_down((dt = 0.2, x = 4.0, f = damped, t = 0.0)...)
function step!(ee::exp_down)
    ee.x = ee.x + ee.f(ee.x) * ee.dt
    ee.t = ee.t + ee.dt
end

plt = plot(1, marker = 3, color=:red)#basically creates empty plot array
plot!(xlims=(0,5), ylims=(0,5))
push!(plt,decay.t,decay.x)
display(plt)
for i in 1:20
    step!(decay)
    push!(plt,decay.t,decay.x)
    display(plt)
end
#=
include("trial_functions.jl")
x = collect(1:10)
y = cos.(x)

plot(x,y,title="hi",color="red")

testing(3,4)
=#
# define the Lorenz attractor

#=
mutable struct Lorenz
    dt; σ; ρ; β; x; y; z
end

function step!(l::Lorenz)
    dx = l.σ*(l.y - l.x)       ; l.x += l.dt * dx
    dy = l.x*(l.ρ - l.z) - l.y ; l.y += l.dt * dy
    dz = l.x*l.y - l.β*l.z     ; l.z += l.dt * dz
end

attractor = Lorenz((dt = 0.02, σ = 10., ρ = 28., β = 8//3, x = 1., y = 1., z = 1.)...)


# initialize a 3D plot with 1 empty series
plt = plot3d(1, xlim=(-25,25), ylim=(-25,25), zlim=(0,50),
                title = "Lorenz Attractor", marker = 3)
=#
# build an animated gif by pushing new points to the plot, saving every 10th frame
#=
@gif for i=1:1500
    step!(attractor)
    push!(plt, attractor.x, attractor.y, attractor.z)
    plot!()
end every 10
=#
#=
for i=1:25
    step!(attractor)
    push!(plt, attractor.x, attractor.y, attractor.z)
    display(plt)
end
=#
