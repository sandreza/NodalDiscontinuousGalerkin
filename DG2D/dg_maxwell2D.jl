include("dg2D.jl")

"""
dg_maxwell!(u̇, u, params, t)

# Description

    numerical solution to 1D maxwell's equation

# Arguments

-   `u̇ = (Eʰ, Hʰ)`: container for numerical solutions to fields
-   `u  = (E , H )`: container for starting field values
-   `params = (𝒢, E, H, ext)`: mesh, E sol, H sol, and material parameters
-   `t`: time to evaluate at

"""
function dg_maxwell2D!(u̇, u, params, t)
    # unpack params
    𝒢  = params[1] # grid parameters
    Hˣ = params[2] # internal parameters for E
    Hʸ = params[3] # internal parameters for H
    Eᶻ = params[4] # external parameters

    # define field differences at faces
    dHˣ = similar(Hˣ.flux)
    @. dHˣ[:] = Hˣ.u[𝒢.vmap⁻] - Hˣ.u[𝒢.vmap⁺]
    dHʸ = similar(Hʸ.flux)
    @. dHʸ[:] = Hʸ.u[𝒢.vmap⁻] - Hʸ.u[𝒢.vmap⁺]
    dEᶻ = similar(Eᶻ.flux)
    @. dEᶻ[:] = Eᶻ.u[𝒢.vmap⁻] - Eᶻ.u[𝒢.vmap⁺]

    # impose reflective BC
    dHˣ[𝒢.mapᴮ] = @. 0
    dHʸ[𝒢.mapᴮ] = @. 0
    dEᶻ[𝒢.mapᴮ] = 2*Eᶻ.u[𝒢.vmapᴮ]

    # perform calculations over elements
    for Ω in 𝒢.Ω
        # evaluate upwind fluxes
        α = 1
        n̂ˣ = Ω.n̂[:,1]
        n̂ʸ = Ω.n̂[:,2]
        n̂⨂dH = n̂ˣ * dHˣ + n̂ʸ * dHʸ
        @. Hˣ.flux =  n̂ʸ * dEᶻ + α * (n̂ˣ * n̂⨂dH - dHˣ)
        @. Hʸ.flux = -n̂ˣ * dEᶻ + α * (n̂ʸ * n̂⨂dH - dHʸ)
        @. Eᶻ.flux = -n̂ˣ * dHʸ + n̂ʸ * dHˣ - α * dEᶻ

        # local derivatives of the fields
        dˣEᶻ,dʸEᶻ = ∇(Eᶻ, Ω)
        ∇⨂H = ∇⨂(Hˣ, Hʸ, Ω)

        # compute RHS of PDE's
        Hˣ.u̇ += -dʸEᶻ + 1//2 * Ω.lift * (Ω.volume .* Hˣ.flux)
        Hʸ.u̇ +=  dˣEᶻ + 1//2 * Ω.lift * (Ω.volume .* Hʸ.flux)
        Eᶻ.u̇ +=  ∇⨂H + 1//2 * Ω.lift * (Ω.volume .* Eᶻ.flux)
    end
    
    return nothing
end
