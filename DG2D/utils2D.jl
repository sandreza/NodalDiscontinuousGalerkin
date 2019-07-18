using Plots

"""
partials(r̃ˣ)

# Description

    Convert array of jacobian matrices to four arrays of individual partial derivatives

# Arguments

-   `r̃ˣ`: array of matrices to convert

# Return Values

-   `rˣ`: array of [1,1] entries
-   `sˣ`: array of [2,1] entries
-   `rʸ`: array of [1,2] entries
-   `sʸ`: array of [2,2] entries

"""
function partials(r̃ˣ)
    # pull partials out from Jacobian
    rˣ = r̃ˣ[:,1,1]
    sˣ = r̃ˣ[:,2,1]
    rʸ = r̃ˣ[:,1,2]
    sʸ = r̃ˣ[:,2,2]

    return rˣ,sˣ,rʸ,sʸ
end


"""
∇!(uˣ, uʸ, u, Ω)

# Description

    Compute gradient of u wrt physical grid

# Arguments

-   `uˣ`: where to store first component of the gradient
-   `uʸ`: where to store second component of the gradient
-   `u`: scalar to take gradient of
-   `Ω`: element to compute in

# Return Values

"""
function ∇!(uˣ,uʸ, u, Ω)
    # compute partial derivatives on ideal grid
    uʳ = Ω.D[1] * u
    uˢ = Ω.D[2] * u

    # pull partials out from Jacobian
    rˣ,sˣ,rʸ,sʸ = partials(Ω.rˣ)

    # compute partial derivatives on physical grid
    @. uˣ = rˣ * uʳ + sˣ * uˢ
    @. uʸ = rʸ * uʳ + sʸ * uˢ

    return nothing
end

"""
∇⨀!(∇⨀u, uˣ, uʸ, Ω)

# Description

    Compute the divergence of u=(uˣ,uʸ) wrt physical grid

# Arguments

-   `∇⨀u`: place to store the divergence of u
-   `uˣ`: first component of vector u
-   `uʸ`: second component of vector u
-   `Ω`: element to compute in

# Return Values

"""
function ∇⨀!(∇⨀u, uˣ, uʸ, Ω)
    # compute partial derivatives on ideal grid
    xʳ = Ω.D[1] * uˣ
    xˢ = Ω.D[2] * uˣ
    yʳ = Ω.D[1] * uʸ
    yˢ = Ω.D[2] * uʸ

    # pull partials out from Jacobian
    rˣ,sˣ,rʸ,sʸ = partials(Ω.rˣ)

    # compute gradient on physical grid
    @. ∇⨀u = rˣ * xʳ + sˣ * xˢ + rʸ * yʳ + sʸ * yˢ

    return nothing
end


"""
∇⨂!(∇⨂u, uˣ, uʸ, Ω)

# Description

    Compute the curl of u=(uˣ,uʸ) wrt physical grid

# Arguments

-   `∇⨂u`: place to store the curl of u
-   `uˣ`: first component of vector u
-   `uʸ`: second component of vector u
-   `Ω`: element to compute in

# Return Values

"""
function ∇⨂!(∇⨂u, uˣ, uʸ, Ω)
    # compute partial derivatives on ideal grid
    xʳ = Ω.D[1] * uˣ
    xˢ = Ω.D[2] * uˣ
    yʳ = Ω.D[1] * uʸ
    yˢ = Ω.D[2] * uʸ

    # pull partials out from Jacobian
    rˣ,sˣ,rʸ,sʸ = partials(Ω.rˣ)

    # compute gradient on physical grid
    @. ∇⨂u = rˣ * yʳ + sˣ * yˢ - rʸ * xʳ - sʸ * xˢ

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
    scatter!(x[𝒢.nodes⁻] , y[𝒢.nodes⁻], color = "black", legend = false)

    # plot boundary of domain
    scatter!(x[𝒢.nodesᴮ] , y[𝒢.nodesᴮ], color = "yellow", legend = false)

    # display
    display(plot(grid))
end


"""
plotfield2D(times, solutions, x, y)

# Description

    Plots the fields as a function of time

# Arguments

-   `times`: time steps to plot
-   `solutions`: fields to plot
-   `x`: x coordinates of the GL points
-   `y`: y coordinates of the GL points

# Return Values

    Displays a plot

"""
function plotfield2D(times, solutions, x, y)
    @animate for t in times
        plots = []
        for (i,sol) in enumerate(solutions)
            ploti = surface(x[:],y[:],sol[t][:], camera = (0,90))# (15,60))
            push!(plots, ploti)
        end
        display(plot(plots...))
    end
end
