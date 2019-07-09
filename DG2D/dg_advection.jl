# first define the stream function
n = 3

#define stream function and components of velocity
ψ(x, y, γ) = exp(γ*(y-1)^2 ) * cos(π/2 * x) * cos(π/2 * y)
u1(x, y, γ) = - π / 2 * sin(π/2 * y) * cos(π/2 * x) * γ * 2 * y * exp(γ*(y-1)^2 )
u2(x, y, γ) = π / 2 * sin(π/2 * x) * exp(γ*(y-1)^2 ) * cos(π/2 * x) * cos(π/2 * y)

u0(x, y, μ) = exp(-μ * x^2 - μ * (y-0.5)^2)

#=
γ = 1.0
μ = 10.0
pyplot()
p1 = contourf(x[:], y[:], (x, y) -> ψ(x,y,γ), title = "Stream Function", xlabel = "x", ylabel = "y")

p2 = contourf(x[:], y[:], (x, y) -> u0(x,y,μ), title = "Initial Condition", xlabel = "x", ylabel = "y")
display(plot(p1,p2))
=#
println("The number of degrees of freedom are")
println(length(x))


function dg_central_2D!(u̇, u, params, t)
    # unpack params
    𝒢 = params[1] # grid parameters
    ι = params[2] # internal parameters
    ε = params[3] # external parameters
    periodic = params[4]

    # calculate fluxes, assigns memory
    flux1 = ε.v1 .* u
    flux2 = ε.v2 .* u

    # Form field differences at faces
    diffs = reshape( flux1[𝒢.vmapM] - flux1[𝒢.vmapP], size(ι.flux))
    @. ι.flux = 0.5 * diffs *  𝒢.normals
    # now for the other velocity
    diffs = reshape( flux2[𝒢.vmapM] - flux2[𝒢.vmapP], size(ι.flux) )
    @. ι.flux += 0.5 * diffs * 𝒢.normals

    # now for the boundary conditions
    # neumann boundary conditions (reflecting)
    @. ι.flux[mapB] = 2*u[vmapB]

    # rhs of the semi-discerte PDE, ∂ᵗu = -∂ˣ(v1*u) - ∂ʸ(v2*u)
    # compute divergence
    du = ∇⨀(flux1, flux2, 𝒢)
    @. u̇ = - du
    lift = 𝒢.lift * (𝒢.fscale .* ι.flux )
    @. u̇ += lift
    return nothing
end
