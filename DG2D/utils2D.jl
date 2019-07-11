using Plots

"""
∇!(uˣ, uʸ, u, Ω)

# Description

    Compute gradient of u wrt physical grid

# Arguments

-   `uˣ`: first component of the gradient, overwitten
-   `uʸ`: second component of the gradient, overwritten
-   `u`: scalar to take gradient of
-   `Ω`: element to compute in

# Return Values



"""
function ∇!(uˣ, uʸ, u, Ω)
    # compute partial derivatives on ideal grid
    uʳ = Ω.Dʳ * u
    uˢ = Ω.Dˢ * u

    # compute partial derivatives on physical grid
    @. uˣ =  Ω.rx * uʳ + Ω.sx * uˢ
    @. uʸ =  Ω.ry * uʳ + Ω.sy * uˢ

    return nothing
end

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

"""
plotgrid2D(𝒢::Grid2D)

# Description

    Plot the GL points, element boundaries, and domain boundaries of a grid

# Arguments

-   `𝒢`: grid to plot

# Return Values

    Displays a plot

"""
function plotgrid2D(𝒢::Grid2D)
    # short cuts
    x = 𝒢.x[:, 1]
    y = 𝒢.x[:, 2]

    # plot the total grid points
    grid = scatter(x, y, legend = false)

    # plot boundary of the elements
    scatter!(x[𝒢.vmap⁻] , y[𝒢.vmap⁻], color = "black", legend = false)

    # plot boundary of domain
    scatter!(x[𝒢.vmapᴮ] , y[𝒢.vmapᴮ], color = "yellow", legend = false)

    # display
    display(plot(grid))
end
