"""
∇(u, Ωᵏ)

# Description

    Compute gradient of u wrt physical grid

# Arguments

-   `u`: scalar to take gradient of
-   `Ωᵏ`: element to compute in

# Return Values

-   `uˣ`: first component of the gradient
-   `uʸ`: second component of the gradient

"""
function ∇(u, Ωᵏ)
    # compute partial derivatives on ideal grid
    uʳ = Ωᵏ.Dʳ * u
    uˢ = Ωᵏ.Dˢ * u

    # compute partial derivatives on physical grid
    uˣ = @. Ωᵏ.rˣ * uʳ + Ωᵏ.sˣ * uˢ
    uʸ = @. Ωᵏ.rʸ * uʳ + Ωᵏ.sʸ * uˢ

    return uˣ,uʸ
end

"""
∇⨀(x, y, Ωᵏ)

# Description

    Compute the divergence of u=(x,y) wrt physical grid

# Arguments

-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ωᵏ`: element to compute in

# Return Values

-   `∇⨀u`: the divergence of u

"""
function ∇⨀(x, y, Ωᵏ)
    # compute partial derivatives on ideal grid
    xʳ = Ωᵏ.Dʳ * x
    xˢ = Ωᵏ.Dˢ * x
    yʳ = Ωᵏ.Dʳ * y
    yˢ = Ωᵏ.Dˢ * y

    # compute gradient on physical grid
    ∇⨀u = @. Ωᵏ.rˣ * xʳ + Ωᵏ.sˣ * xˢ + Ωᵏ.rʸ * yʳ + Ωᵏ.sʸ * yˢ

    return ∇⨀u
end

"""
∇⨂(x, y, Ωᵏ)

# Description

    Compute the curl of u=(x,y) wrt physical grid

# Arguments

-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ωᵏ`: element to compute in

# Return Values

-   `∇⨂u`: the curl of u

"""
function ∇⨂(x, y, Ωᵏ)
    # compute partial derivatives on ideal grid
    xʳ = Ωᵏ.Dʳ * x
    xˢ = Ωᵏ.Dˢ * x
    yʳ = Ωᵏ.Dʳ * y
    yˢ = Ωᵏ.Dˢ * y

    # compute gradient on physical grid
    ∇⨂u = @. Ωᵏ.rˣ * yʳ + Ωᵏ.sˣ * yˢ - Ωᵏ.rʸ * xʳ - Ωᵏ.sʸ * xˢ

    return ∇⨂u
end
