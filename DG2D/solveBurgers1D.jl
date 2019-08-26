include("field2D.jl")
include("utils2D.jl")

"""
solveBurgers1D!(fields, params)

# Description

    numerical solution to Chorin Navier Stokes equation
    in vector form:
    ∂ᵗu = -∇(u²/2) + ν∇²u
    written out component wise for DG formulation:
    ∂ᵗu = -∂ˣ(u²/2 - νuˣ) - ∂ʸ(u²/2 - νuʸ)
    we are setting the flux in the y-direction to be zero for the 1D case


# Arguments

-   `fields = (u)`: velocity in each dimension
-   `auxils = (uˣ, uʸ, u²)`: auxiliary fields for computation
-   `params = (𝒢, ν, α, β)`: grid struct, viscosity, nonlinear switch, and 2D switch
-   `t`: time to compute BC at

"""
function solveBurgers1D!(fields, auxil, params, t)
    # unpack params
    𝒢 = params[1] # grid parameters
    ν = params[2]
    α = params[3]
    β = params[4]

    # unpack fields
    u  = fields[1]

    # auxiliary fields
    u² = auxil[1]
    uˣ = auxil[2]
    uʸ = auxil[3]

    for Ωᵏ in 𝒢.Ω
        # get volume nodes
        iⱽ = Ωᵏ.iⱽ

        # compute volume contribution to uˣ and uʸ
        ∇!(u.φˣ, u.φʸ, u.ϕ, Ωᵏ)
        @. uˣ.ϕ[iⱽ] = sqrt(ν) * u.φˣ[iⱽ]
        @. uʸ.ϕ[iⱽ] = sqrt(ν) * u.φʸ[iⱽ]

        # define physical fluxes for uˣ and uʸ
        @. uˣ.φˣ[iⱽ] = sqrt(ν) * u.ϕ[iⱽ]
        @. uʸ.φʸ[iⱽ] = sqrt(ν) * u.ϕ[iⱽ]

        # compute surface contributions to uˣ, uʸ
        for f in Ωᵏ.faces
            # get face nodes
            i⁻ = f.i⁻
            i⁺ = f.i⁺

            computeCentralFluxes!(uˣ, f)
            computeCentralFluxes!(uʸ, f)

            # impose BC
            if f.isBoundary[1]
                uᴮ = [u⁰(𝒢.x[i,1],t) for i in f.i⁻]
                @. uˣ.fˣ[f.i⁻] = sqrt(ν) * uᴮ
                @. uʸ.fʸ[f.i⁻] = sqrt(ν) * uᴮ
            end

            computeSurfaceTerms!(uˣ, Ωᵏ, f)
            computeSurfaceTerms!(uʸ, Ωᵏ, f)
        end

        # compute u²
        @. u².ϕ[iⱽ] = u.ϕ[iⱽ]^2

        # define physical fluxes
        @. u.φˣ[iⱽ] = 0.5 * α * u².ϕ[iⱽ] - sqrt(ν) * uˣ.ϕ[iⱽ]
        @. u.φʸ[iⱽ] = 0.5 * β * (α * u².ϕ[iⱽ] - sqrt(ν) * uʸ.ϕ[iⱽ])

        # compute volume contributions
        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ωᵏ)
        @. u.ϕ̇[iⱽ] = -u.𝚽[iⱽ]

        # compute surface contributions to tendency
        for f in Ωᵏ.faces
            computeCentralDifference!(uˣ, f)
            computeCentralDifference!(uʸ, f)
            computeCentralDifference!(u², f)

            # impose BC on uˣ, uʸ, and u²
            if f.isBoundary[1]
                uᴮ = [u⁰(𝒢.x[i,1],t) for i in f.i⁻]
                @. uˣ.ϕ°[f.i⁻] = uˣ.ϕ[f.i⁻]
                @. uʸ.ϕ°[f.i⁻] = uʸ.ϕ[f.i⁻]
                @. u².ϕ°[f.i⁻] = uᴮ^2
            end

            # evaluate numerical flux for u
            C = maximum(abs.(u.ϕ[f.i⁻]))
            @. u.fˣ[f.i⁻] = 0.5 * α * u².ϕ°[f.i⁻] - sqrt(ν) * uˣ.ϕ°[f.i⁻]
            @. u.fʸ[f.i⁻] = 0.5 * β * (α * u².ϕ°[f.i⁻] - sqrt(ν) * uʸ.ϕ°[f.i⁻])
            computeLaxFriedrichsFluxes!(u, f, C)

            computeSurfaceTerms!(u, Ωᵏ, f)
        end
    end

    return nothing
end
