
"""
∇⨀(x, y, Ω)

# Description

    Compute the divergence of u=(x,y) wrt physical grid

# Arguments

-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ω`: element to compute in

# Return Values

-   `∇⨀u`: the divergence of u

"""
function ∇⨀(x, y, Ω)
    # compute partial derivatives on ideal grid
    xʳ = Ω.Dʳ * x
    xˢ = Ω.Dˢ * x
    yʳ = Ω.Dʳ * y
    yˢ = Ω.Dˢ * y

    # compute gradient on physical grid
    ∇⨀u = @. Ω.rˣ * xʳ + Ω.sˣ * xˢ + Ω.rʸ * yʳ + Ω.sʸ * yˢ

    return ∇⨀u
end

"""
∇⨂(x, y, Ω)

# Description

    Compute the curl of u=(x,y) wrt physical grid

# Arguments

-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ω`: element to compute in

# Return Values

-   `∇⨂u`: the curl of u

"""
function ∇⨂(x, y, Ω)
    # compute partial derivatives on ideal grid
    xʳ = Ω.Dʳ * x
    xˢ = Ω.Dˢ * x
    yʳ = Ω.Dʳ * y
    yˢ = Ω.Dˢ * y

    # compute gradient on physical grid
    ∇⨂u = @. Ω.rˣ * yʳ + Ω.sˣ * yˢ - Ω.rʸ * xʳ - Ω.sʸ * xˢ

    return ∇⨂u
end




"""
∇⨀!(∇⨀u, fx, fy, Ω)

# Description

    Compute the divergence of u=(fx,fy) wrt physical grid

# Arguments
-   `∇⨀u`: allocated memory for result
-   `x`: first component of vector u
-   `y`: second component of vector u
-   `Ω`: element to compute in

# Return Values

-   `∇⨀u`: the divergence of u

"""
function ∇⨀!(∇⨀u, x, y, Ω)
    # compute partial derivatives on ideal grid
    xʳ = Ω.Dʳ * x
    xˢ = Ω.Dˢ * x
    yʳ = Ω.Dʳ * y
    yˢ = Ω.Dˢ * y

    # compute gradient on physical grid
    @. ∇⨀u = Ω.rx * xʳ + Ω.sx * xˢ + Ω.ry * yʳ + Ω.sy * yˢ
    return nothing
end
