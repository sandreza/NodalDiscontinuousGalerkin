include("utils2D.jl")

function stommel!(ϕ̇, ϕ, params, t)
    # unpack params
    𝒢 = params[1] # grid parameters
    ε = params[2]  # external parameters
    𝓊 = params[3] # internal parameters
    𝓋 = params[4] # internal parameters
    𝒽 = params[4] # internal parameters

    u = ϕ[1]
    v = ϕ[2]
    η = ϕ[3]


    # calculate fluxes for 𝓊
    @. 𝓊.φˣ = η
    @. 𝓊.φʸ = 0.0

    # calculate fluxes for 𝓋
    @. 𝓋.φˣ = 0.0
    @. 𝓋.φʸ = η

    # calculate fluxes for 𝒽
    @. 𝒽.φˣ = u
    @. 𝒽.φʸ = v


    # now for the boundary conditions
    @. 𝓊.fⁿ[𝒢.mapB] =  -u[𝒢.vmapB]
    @. 𝓋.fⁿ[𝒢.mapB] =  -v[𝒢.vmapB]

    # Form field differences at faces, computing central flux
    @. 𝓊.fˣ[:] = (𝓊.φˣ[𝒢.vmapM] - 𝓊.φˣ[𝒢.vmapP])/2
    @. 𝓊.fʸ[:] = (𝓊.φʸ[𝒢.vmapM] - 𝓊.φʸ[𝒢.vmapP])/2
    #now for the normal component along the faces
    @. 𝓊.fⁿ = 𝓊.fˣ * 𝒢.nx + 𝓊.fʸ * 𝒢.ny

    # Form field differences at faces, computing central flux
    @. 𝓋.fˣ[:] = (𝓋.φˣ[𝒢.vmapM] - 𝓋.φˣ[𝒢.vmapP])/2
    @. 𝓋.fʸ[:] = (𝓋.φʸ[𝒢.vmapM] - 𝓋.φʸ[𝒢.vmapP])/2
    #now for the normal component along the faces
    @. 𝓋.fⁿ = 𝓋.fˣ * 𝒢.nx + 𝓋.fʸ * 𝒢.ny

    # Form field differences at faces, computing central flux
    @. 𝒽.fˣ[:] = (𝒽.φˣ[𝒢.vmapM] - 𝒽.φˣ[𝒢.vmapP])/2
    @. 𝒽.fʸ[:] = (𝒽.φʸ[𝒢.vmapM] - 𝒽.φʸ[𝒢.vmapP])/2
    #now for the normal component along the faces
    @. 𝒽.fⁿ = 𝒽.fˣ * 𝒢.nx + 𝒽.fʸ * 𝒢.ny


    # rhs of the semi-discrete PDE, ∂ᵗu = -∂ˣ(v1*u) - ∂ʸ(v2*u)
    # compute divergence
    ∇⨀!(u̇, ι.φˣ, ι.φʸ, 𝒢)
    @. u̇ *= -1.0
    lift = 𝒢.lift * (𝒢.fscale .* ι.fⁿ)
    @. u̇ +=  lift #inefficient part, has to be done pointwise
        # now hack in zeroness on boundary
    return nothing
end
