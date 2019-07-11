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
make_periodic2D(Ω)

# Description

- Takes a rectangular grid and modifies vmapP so that the domain becomes periodic

# Arguments

- `Ω` : the mesh struct

# Return : nothing


"""
function make_periodic2D(Ω)
    boundary_index = findall(grid.vmapM - grid.vmapP .≈ 0.0)
    ax = minimum(grid.x)
    bx = maximum(grid.x)
    ay = minimum(grid.y)
    by = maximum(grid.y)
    leftface_index   =  findall(grid.x[:] .== ax)
    rightface_index  =  findall(grid.x[:] .== bx)
    bottomface_index =  findall(grid.y[:] .== ay)
    topface_index    =  findall(grid.y[:] .== by)

    for j in boundary_index
        global_index = grid.vmapM[j]
        x = grid.x[j]
        y = grid.y[j]
        if j in leftface_index
            findall(grid.y[rightface] .== y)
        end
    end


end
