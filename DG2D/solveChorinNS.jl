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
function solveChorinNS!(fields, fluxes, auxils, params, t)
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

    # fluxes
    φᵘ  = fluxes[1]
    φᵛ  = fluxes[2]
    φˣᵤ = fluxes[3]
    φʸᵤ = fluxes[4]
    φˣᵥ = fluxes[5]
    φʸᵥ = fluxes[6]

    # for convenience
    nonlinear   = [uu, uv, vu, vv]
    derivatives = [uˣ, uʸ, vˣ, vʸ]

    # compute volume contributions to first derivatives
    for Ω in 𝒢.Ω
        # define physical fluxes for first derivatives
        computePhysicalFlux!(uˣ.φˣ, φᵘ, Ω)
        computePhysicalFlux!(uʸ.φʸ, φᵘ, Ω)

        computePhysicalFlux!(vˣ.φˣ, φᵛ, Ω)
        computePhysicalFlux!(vʸ.φʸ, φᵛ, Ω)

        # volume contribs
        ∇!(u.φˣ, u.φʸ, u.ϕ, Ω)
        @. uˣ.ϕ[Ω.iⱽ] = u.φˣ[Ω.iⱽ]
        @. uʸ.ϕ[Ω.iⱽ] = u.φʸ[Ω.iⱽ]

        # volume contribs
        ∇!(v.φˣ, v.φʸ, v.ϕ, Ω)
        @. vˣ.ϕ[Ω.iⱽ] = v.φˣ[Ω.iⱽ]
        @. vʸ.ϕ[Ω.iⱽ] = v.φʸ[Ω.iⱽ]
    end

    # compute surface contributions to first derivatives
    for Ω in 𝒢.Ω
        for f in Ω.faces
            computeCentralDifference!(u, f)
            computeCentralDifference!(v, f)

            # impose BC
            if f.isBoundary[1]
                uᴮ = [u⁰(𝒢.x[i],t) for i in f.i⁻]
                vᴮ = [v⁰(𝒢.x[i],t) for i in f.i⁻]

                @. u.ϕ°[f.i⁻] = uᴮ
                @. v.ϕ°[f.i⁻] = vᴮ
            end

            computeNumericalFlux!(uˣ.fˣ, φᵘ, f)
            computeNumericalFlux!(uʸ.fʸ, φᵘ, f)

            computeNumericalFlux!(vˣ.fˣ, φᵛ, f)
            computeNumericalFlux!(vʸ.fʸ, φᵛ, f)

            for 𝑓 in derivatives
                computeSurfaceTerms!(𝑓.ϕ, 𝑓, Ω, f)
            end
        end
    end

    # compute non-linear terms
    @. uu.ϕ = u.ϕ * u.ϕ
    @. uv.ϕ = u.ϕ * v.ϕ
    @. vu.ϕ = v.ϕ * u.ϕ
    @. vv.ϕ = v.ϕ * v.ϕ

    # compute volume contributions to the tendecies
    for Ω in 𝒢.Ω
        computePhysicalFlux!(u.φˣ, φˣᵤ, Ω)
        computePhysicalFlux!(u.φʸ, φʸᵤ, Ω)

        computePhysicalFlux!(v.φˣ, φˣᵥ, Ω)
        computePhysicalFlux!(v.φʸ, φʸᵥ, Ω)

        ∇⨀!(u.𝚽, u.φˣ, u.φʸ, Ω)
        ∇⨀!(v.𝚽, v.φˣ, v.φʸ, Ω)

        @. u.ϕ̇[Ω.iⱽ] = u.𝚽[Ω.iⱽ]
        @. v.ϕ̇[Ω.iⱽ] = v.𝚽[Ω.iⱽ]
    end

    # compute surface contributions to tendency
    for Ω in 𝒢.Ω
        for f in Ω.faces
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

            computeNumericalFlux!(u.fˣ, φˣᵤ, f)
            computeNumericalFlux!(u.fʸ, φʸᵤ, f)
            computeNumericalFlux!(v.fˣ, φˣᵥ, f)
            computeNumericalFlux!(v.fʸ, φʸᵥ, f)

            ṽ⁻ = @. abs(f.nˣ * u.ϕ[f.i⁻] + f.nʸ * v.ϕ[f.i⁻])
            ṽ⁺ = @. abs(f.nˣ * u.ϕ[f.i⁺] + f.nʸ * v.ϕ[f.i⁺])
            C = -maximum([ṽ⁻, ṽ⁺])

            computeLaxFriedrichsFluxes!(u, f, C)
            computeLaxFriedrichsFluxes!(v, f, C)

            computeSurfaceTerms!(u.ϕ̇, u, Ω, f)
            computeSurfaceTerms!(v.ϕ̇, v, Ω, f)
        end
    end

    return nothing
end
