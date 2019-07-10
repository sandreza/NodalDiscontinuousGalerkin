
println("The number of degrees of freedom are")
println(length(x))
include("utils2D.jl")

function dg_central_2D!(u̇, u, params, t)
    # unpack params
    𝒢 = params[1] # grid parameters
    ι = params[2] # internal parameters
    ε = params[3] # external parameters

    # calculate fluxes, assigns memory
    @. ι.φˣ = ε.v1 * u
    @. ι.φʸ = ε.v2 * u

    # Form field differences at faces, computing central flux
    @. ι.fˣ[:] = (ι.φˣ[𝒢.vmapM] - ι.φˣ[𝒢.vmapP])/2
    @. ι.fʸ[:] = (ι.φʸ[𝒢.vmapM] - ι.φʸ[𝒢.vmapP])/2
    #now for the normal component along the faces
    @. ι.fⁿ = ι.fˣ * 𝒢.nx + ι.fʸ * 𝒢.ny

    # now for the boundary conditions
    # neumann boundary conditions (reflecting)
    @. ι.fⁿ[𝒢.mapB] = 2*u[𝒢.vmapB]

    # rhs of the semi-discrete PDE, ∂ᵗu = ∂ˣ(v1*u) + ∂ʸ(v2*u)
    # compute divergence
    ∇⨀!(u̇, ι.φˣ, ι.φʸ, 𝒢)
    lift = 𝒢.lift * (𝒢.fscale .* ι.fⁿ) #inefficient part
    @. u̇ += lift
    return nothing
end
