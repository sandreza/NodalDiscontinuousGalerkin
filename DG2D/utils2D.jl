using Plots

"""
partials(rˣ)

# Description

    Convert array of jacobian matrices to four arrays of individual partial derivatives

# Arguments

-   `rˣ`: array of matrices to convert

# Return Values

-   `rˣ`: array of [1,1] entries
-   `sˣ`: array of [2,1] entries
-   `rʸ`: array of [1,2] entries
-   `sʸ`: array of [2,2] entries

"""
function partials(rˣ)
    # pull partials out from Jacobian
    rˣ = rˣ[:,1,1]
    sˣ = rˣ[:,2,1]
    rʸ = rˣ[:,1,2]
    sʸ = rˣ[:,2,2]

    return rˣ,sˣ,rʸ,sʸ
end


"""
∇!(uˣ, uʸ, u, Ω)

# Description

    Compute gradient of u wrt physical grid

# Arguments

-   `u`: scalar to take gradient of
-   `Ω`: element to compute in

# Return Values

-   `uˣ`: first component of the gradient
-   `uʸ`: second component of the gradient

"""
function ∇(u, Ω)
    # compute partial derivatives on ideal grid
    uʳ = Ω.D[1] * u
    uˢ = Ω.D[2] * u

    # pull partials out from Jacobian
    rˣ,sˣ,rʸ,sʸ = partials(Ω.rˣ)

    # compute partial derivatives on physical grid
    @. uˣ = rˣ * uʳ + sˣ * uˢ
    @. uʸ = rʸ * uʳ + sʸ * uˢ

    return uˣ,uʸ
end

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
    xʳ = Ω.D[1] * x
    xˢ = Ω.D[2] * x
    yʳ = Ω.D[1] * y
    yˢ = Ω.D[2] * y

    # pull partials out from Jacobian
    rˣ,sˣ,rʸ,sʸ = partials(Ω.rˣ)

    # compute gradient on physical grid
    ∇⨀u = @. rˣ * xʳ + sˣ * xˢ + rʸ * yʳ + sʸ * yˢ

    return ∇⨀u
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
    xʳ = Ω.D[1] * x
    xˢ = Ω.D[2] * x
    yʳ = Ω.D[1] * y
    yˢ = Ω.D[2] * y

    # pull partials out from Jacobian
    rˣ,sˣ,rʸ,sʸ = partials(Ω.rˣ)

    # compute gradient on physical grid
    ∇⨂u = @. rˣ * yʳ + sˣ * yˢ - rʸ * xʳ - sʸ * xˢ

    return ∇⨂u
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
    # whole grid plotting
    x = 𝒢.x[:, 1]
    y = 𝒢.x[:, 2]

    # initial grid (mainly for the to make for loop simpler)
    grid = scatter(x, y, legend = false)

    # plot GL points elementwise
    for Ω in 𝒢.Ω
        r = Ω.x[:, 1]
        s = Ω.x[:, 2]

        scatter!(r, s, legend = false)
    end

    # plot boundary of the elements
    scatter!(x[𝒢.vmap⁻] , y[𝒢.vmap⁻], color = "black", legend = false)

    # plot boundary of domain
    scatter!(x[𝒢.vmapᴮ] , y[𝒢.vmapᴮ], color = "yellow", legend = false)

    # display
    display(plot(grid))
end
