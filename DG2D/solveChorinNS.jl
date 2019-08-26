include("field2D.jl")
include("utils2D.jl")

"""
solveChorinNS!(fields, auxils, params, time)

# Description

    numerical solution to Chorin Navier Stokes equation
    in vector form:
    ∂ᵗũ = -∇(ũ⨂ũ) + c²∇(∇⋅ũ) + ν∇²ũ
    written out component wise for DG formulation:
    ∂ᵗu = -∂ˣ(uu - (ν+c²)uˣ - c²vʸ) - ∂ʸ(uv - νuʸ)
    ∂ᵗv = -∂ʸ(vv - (ν+c²)vʸ - c²uˣ) - ∂ˣ(vu - νvˣ)


# Arguments

-   `fields = (u, v)`: velocity in each dimension
-   `auxils = (uˣ, uʸ, vˣ, vʸ, uu, uv, vu, vv)`: auxiliary fields for computation
-   `params = (𝒢, ν, c²)`: grid struct, viscosity, speed of sound, and nonlinear switch
-   `t`: time to compute BC at

"""
function solveChorinNS!(fields, auxils, params, t)
    # unpack parameters
    𝒢  = params[1]
    ν  = params[2]
    c² = params[3]
    α  = params[4]

    # main velocity fields
    u  = fields[1]
    v  = fields[2]

    # utility fields for first derivatives
    uˣ = auxils[1]
    uʸ = auxils[2]
    vˣ = auxils[3]
    vʸ = auxils[4]

    # utility fields for second order terms
    uu = auxils[5]
    uv = auxils[6]
    vu = auxils[7]
    vv = auxils[8]

    # for convenience
    nonlinear   = [uu, uv, vu, vv]
    derivatives = [uˣ, uʸ, vˣ, vʸ]

    # compute volume contributions to first derivatives
    for Ωᵏ in 𝒢.Ω
        # get volume nodes
        iⱽ = Ωᵏ.iⱽ

        # define physical fluxes for first derivatives
        @. uˣ.φˣ[iⱽ] = u.ϕ[iⱽ]
        @. uʸ.φʸ[iⱽ] = u.ϕ[iⱽ]

        @. vˣ.φˣ[iⱽ] = v.ϕ[iⱽ]
        @. vʸ.φʸ[iⱽ] = v.ϕ[iⱽ]

        ∇!(u.φˣ, u.φʸ, u.ϕ, Ωᵏ)
        @. uˣ.ϕ[iⱽ] = u.φˣ[iⱽ]
        @. uʸ.ϕ[iⱽ] = u.φʸ[iⱽ]

        ∇!(v.φˣ, v.φʸ, v.ϕ, Ωᵏ)
        @. vˣ.ϕ[iⱽ] = v.φˣ[iⱽ]
        @. vʸ.ϕ[iⱽ] = v.φʸ[iⱽ]
    end

    # compute surface contributions to first derivatives
    for Ωᵏ in 𝒢.Ω
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
    end

    # compute volume contributions to the tendecies
    for Ωᵏ in 𝒢.Ω
        # get volume nodes
        iⱽ = Ωᵏ.iⱽ
        # compute non-linear terms
        @. uu.ϕ[iⱽ] = u.ϕ[iⱽ] * u.ϕ[iⱽ]
        @. uv.ϕ[iⱽ] = u.ϕ[iⱽ] * v.ϕ[iⱽ]
        @. vu.ϕ[iⱽ] = v.ϕ[iⱽ] * u.ϕ[iⱽ]
        @. vv.ϕ[iⱽ] = v.ϕ[iⱽ] * v.ϕ[iⱽ]

        # define physical fluxes for u and v
        @. u.φˣ[iⱽ] = α * uu.ϕ[iⱽ] - (ν+c²) * uˣ.ϕ[iⱽ] - c² * vʸ.ϕ[iⱽ]
        @. u.φʸ[iⱽ] = α * uv.ϕ[iⱽ] - ν * uʸ.ϕ[iⱽ]

        @. v.φˣ[iⱽ] = α * vu.ϕ[iⱽ] - ν * vˣ.ϕ[iⱽ]
        @. v.φʸ[iⱽ] = α * vv.ϕ[iⱽ] - (ν+c²) * vʸ.ϕ[iⱽ] - c² * uˣ.ϕ[iⱽ]

        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ωᵏ)
        @. u.ϕ̇[iⱽ] = -u.𝚽[iⱽ]

        ∇⨀!(v.𝚽, v.φˣ, v.φʸ, Ωᵏ)
        @. v.ϕ̇[iⱽ] = -v.𝚽[iⱽ]
    end

    # compute surface contributions to tendency
    for Ωᵏ in 𝒢.Ω
        for f in Ωᵏ.faces
            for 𝑓 in auxils
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

            ṽ⁻ = @. abs(f.nˣ * u.ϕ[f.i⁻] + f.nʸ * v.ϕ[f.i⁻])
            ṽ⁺ = @. abs(f.nˣ * u.ϕ[f.i⁺] + f.nʸ * v.ϕ[f.i⁺])
            C = maximum([ṽ⁻, ṽ⁺])

            @. u.fˣ[f.i⁻] = α * uu.ϕ°[f.i⁻] - (ν+c²) * uˣ.ϕ°[f.i⁻] - c² * vʸ.ϕ°[f.i⁻]
            @. u.fʸ[f.i⁻] = α * uv.ϕ°[f.i⁻] - ν * uʸ.ϕ°[f.i⁻]
            computeLaxFriedrichsFluxes!(u, f, C)
            computeSurfaceTerms!(u, Ωᵏ, f)

            @. v.fˣ[f.i⁻] = α * vu.ϕ°[f.i⁻] - ν * vˣ.ϕ°[f.i⁻]
            @. v.fʸ[f.i⁻] = α * vv.ϕ°[f.i⁻] - (ν+c²) * vʸ.ϕ°[f.i⁻] - c² * uˣ.ϕ°[f.i⁻]
            computeLaxFriedrichsFluxes!(v, f, C)
            computeSurfaceTerms!(v, Ωᵏ, f)
        end
    end

    return nothing
end
