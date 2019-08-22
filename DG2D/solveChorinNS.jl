include("field2D.jl")
include("utils2D.jl")

"""
solveSalmonCNS!(fields, params)

# Description

    numerical solution to Chorin Navier Stokes equation
    in vector form:
    ∂ᵗũ = -∇(ũ⨂ũ) + c²∇(∇⋅ũ) + ν∇²ũ
    written out component wise for DG formulation:
    ∂ᵗu = -∂ˣ(uu - (ν+c²)uˣ - c²vʸ) - ∂ʸ(uv - νuʸ)
    ∂ᵗv = -∂ʸ(vv - (ν+c²)vʸ - c²uˣ) - ∂ˣ(vu - νvˣ)


# Arguments

-   `fields = (u, v)`: velocity in each dimension
-   `params = (𝒢, ν, c²)`: grid struct, viscosity, and speed of sound
-   `BCᵈ = (Dᵘ, Dᵛ)`: dirichlet boundary conditions for each field
-   `BCⁿ = (Nᵘ, Nᵛ)`:   neumann boundary conditions for each field

"""
function solveChorinNS!(fields, params, time; BCᵈ = [nothing, nothing, nothing], BCⁿ = [nothing, nothing, nothing])
    # unpack parameters
    𝒢  = params[1]
    ν  = params[2]
    c² = params[3]

    # main velocity fields
    u  = fields[1]
    v  = fields[2]

    # utility fields for first derivatives
    uˣ = fields[3]
    uʸ = fields[4]
    vˣ = fields[5]
    vʸ = fields[6]

    # utility fields for second order terms
    uu = fields[7]
    uv = fields[8]
    vu = fields[9]
    vv = fields[10]

    # for convenience
    nonlinear   = [uu, uv, vu, vv]
    derivatives = [uˣ, uʸ, vˣ, vʸ]
    auxiliary   = nonlinear + derivatives

    for Ωᵏ in 𝒢.Ω
        # get volume nodes
        iⱽ = Ωᵏ.iⱽ

        # compute volume contributions to first derivatives
        ∇!(u.φˣ, u.φʸ, u.ϕ, Ωᵏ)
        @. uˣ.ϕ[iⱽ] = u.φˣ[iⱽ]
        @. uʸ.ϕ[iⱽ] = u.φʸ[iⱽ]

        ∇!(v.φˣ, v.φʸ, v.ϕ, Ωᵏ)
        @. vˣ.ϕ[iⱽ] = v.φˣ[iⱽ]
        @. vʸ.ϕ[iⱽ] = v.φʸ[iⱽ]

        # define physical fluxes for first derivatives
        @. uˣ.φˣ[iⱽ] = u.ϕ[iⱽ]
        @. uʸ.φʸ[iⱽ] = u.ϕ[iⱽ]

        @. vˣ.φˣ[iⱽ] = v.ϕ[iⱽ]
        @. vʸ.φʸ[iⱽ] = v.ϕ[iⱽ]

        # compute surface contributions to first derivatives
        for f in Ωᵏ.faces
            for 𝑓 in derivatives
                computeCentralFluxes!(𝑓, f)
            end

            # impose BC
            if f.isBoundary[1]
                uᴮ = [u⁰(𝒢.x[i],t) for i in f.i⁻]
                @. uˣ.fˣ[f.i⁻] = uᴮ
                @. uʸ.fʸ[f.i⁻] = uᴮ

                vᴮ = [v⁰(𝒢.x[i],t) for i in f.i⁻]
                @. vˣ.fˣ[f.i⁻] = vᴮ
                @. vʸ.fʸ[f.i⁻] = vᴮ
            end

            for 𝑓 in derivatives
                computeSurfaceTerms!(𝑓, Ωᵏ, f)
            end
        end

        # compute non-linear terms
        @. uu.ϕ[iⱽ] = u.ϕ[iⱽ] * u.ϕ[iⱽ]
        @. uv.ϕ[iⱽ] = u.ϕ[iⱽ] * v.ϕ[iⱽ]
        @. vu.ϕ[iⱽ] = v.ϕ[iⱽ] * u.ϕ[iⱽ]
        @. vv.ϕ[iⱽ] = v.ϕ[iⱽ] * v.ϕ[iⱽ]

        # define physical fluxes for u and v
        @. u.φˣ[iⱽ] = uu.ϕ[iⱽ] - (ν+c²) * uˣ.ϕ[iⱽ] - c² * vʸ.ϕ[iⱽ]
        @. u.φʸ[iⱽ] = uv.ϕ[iⱽ] - ν * uʸ.ϕ[iⱽ]

        @. v.φˣ[iⱽ] = vu.ϕ[iⱽ] - ν * vˣ.ϕ[iⱽ]
        @. v.φʸ[iⱽ] = vv.ϕ[iⱽ] - (ν+c²) * vʸ.ϕ[iⱽ] - c² * uˣ.ϕ[iⱽ]

        # compute volume contributions to the tendecies
        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ωᵏ)
        @. u.ϕ̇[iⱽ] = -u.𝚽[iⱽ]

        ∇⨀!(v.𝚽, v.φˣ, v.φʸ, Ωᵏ)
        @. v.ϕ̇[iⱽ] = -v.𝚽[iⱽ]

        # compute surface contributions to tendency
        for f in Ωᵏ.faces
            for 𝑓 in auxiliary
                computeCentralDifference!(𝑓, f)
            end

            # impose BC on auxiliary fields
            if f.isBoundary[1]
                uᴮ = [u⁰(𝒢.x[i],t) for i in f.i⁻]
                vᴮ = [v⁰(𝒢.x[i],t) for i in f.i⁻]
                @. uu.ϕ°[f.i⁻] = uᴮ * uᴮ
                @. uv.ϕ°[f.i⁻] = uᴮ * vᴮ
                @. vu.ϕ°[f.i⁻] = vᴮ * uᴮ
                @. vv.ϕ°[f.i⁻] = vᴮ * vᴮ

                @. uˣ.ϕ°[f.i⁻] = uˣ.ϕ[f.i⁻]
                @. uʸ.ϕ°[f.i⁻] = uʸ.ϕ[f.i⁻]
                @. vˣ.ϕ°[f.i⁻] = vˣ.ϕ[f.i⁻]
                @. vʸ.ϕ°[f.i⁻] = vʸ.ϕ[f.i⁻]
            end

            Cᵘ = []
            @. u.fˣ[f.i⁻] = uu.ϕ°[f.i⁻] - (ν+c²) * uˣ.ϕ°[f.i⁻] - c² * vʸ.ϕ°[f.i⁻]
            @. u.fʸ[f.i⁻] = uv.ϕ°[f.i⁻] - ν * uʸ.ϕ°[f.i⁻]
            computeLaxFriedrichsFluxes!(u, f, Cᵘ)
            computeSurfaceTerms!(u, Ωᵏ, f)

            Cᵛ = []
            @. v.fˣ[f.i⁻] = vu.ϕ°[f.i⁻] - ν * vˣ.ϕ°[f.i⁻]
            @. v.fʸ[f.i⁻] = vv.ϕ°[f.i⁻] - (ν+c²) * vʸ.ϕ°[f.i⁻] - c² * uˣ.ϕ°[f.i⁻]
            computeLaxFriedrichsFluxes!(v, f, Cᵛ)
            computeSurfaceTerms!(v, Ωᵏ, f)
        end
    end

    return nothing
end
