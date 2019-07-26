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
∇⨂∇⨂(u, v, Ω)

# Description

    Compute the curl of curl of u

# Arguments

-   `ux`: first component of vector
-   `uy`: second component of vector
-   `Ω`: element to compute in

# Return Values

-   `(∇⨂∇⨂u)ˣ`: first component the curl of curl u
-   `(∇⨂∇⨂u)ʸ`: second component the curl of curl u

"""
function ∇⨂∇⨂(ux, uy, Ω)
    ∂ˣux = Ω.rx .* ( Ω.Dʳ * ux  + Ω.Dʳ * ux )
    ∂ʸux = Ω.ry .* ( Ω.Dʳ * ux  + Ω.Dʳ * ux )
    ∂ˣuy = Ω.rx .* ( Ω.Dʳ * uy  + Ω.Dʳ * uy )
    ∂ʸuy = Ω.ry .* ( Ω.Dʳ * uy  + Ω.Dʳ * uy )

    ∂ˣ∂ʸux = Ω.ry .* ( Ω.Dʳ * ∂ˣux  + Ω.Dʳ * ∂ˣux )
    ∂ˣ∂ʸuy = Ω.ry .* ( Ω.Dʳ * ∂ˣuy  + Ω.Dʳ * ∂ˣuy )

    ∂ˣ∂ˣuy = Ω.rx .* ( Ω.Dʳ * ∂ˣuy  + Ω.Dʳ * ∂ˣuy )
    ∂ʸ∂ʸux = Ω.ry .* ( Ω.Dʳ * ∂ʸux  + Ω.Dʳ * ∂ʸux )

    tmpˣ = ∂ˣ∂ʸuy - ∂ʸ∂ʸux
    tmpʸ = ∂ˣ∂ʸux - ∂ˣ∂ˣuy
    return tmpˣ , tmpʸ
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
advec(∇⨀u, fx, fy, Ω)

# Description

    Compute the advection of a scalar θ by flow field (vx,vy)

# Arguments
-   `u⨀∇θ`: allocated memory for result
-   `vx`: first component of vector u
-   `vy`: second component of vector u
-   `θ`: the scalar
-   `Ω`: element to compute in

# Return Values

-   `∇⨀u`: the divergence of u

"""
function advec!(u⨀∇θ, vx, vy, θ, Ω)
    # compute gradient on physical grid
    tmpˣ =  Ω.rx .* ( vx .* (Ω.Dʳ * θ) )
    tmpˣ += Ω.sx .* ( vx .* (Ω.Dˢ * θ) )
    tmpʸ =  Ω.ry .* ( vy .* (Ω.Dʳ * θ) )
    tmpʸ += Ω.sy .* ( vy .* (Ω.Dˢ * θ) )

    @. u⨀∇θ = (tmpˣ + tmpʸ)

    return nothing
end


"""
sym_advec(∇⨀u, fx, fy, Ω)

# Description

-    Compute the advection of a scalar θ by flow field (vx,vy), symmetrized advection

# Arguments
-   `u⨀∇θ`: allocated memory for result
-   `vx`: first component of vector u
-   `vy`: second component of vector u
-   `θ`: the scalar
-   `Ω`: mesh to compute in

# Return Values

-   `u⨀∇θ`: symmetric advective component

"""
function sym_advec!(u⨀∇θ, vx, vy, θ, Ω)

    # compute gradient on physical grid
    tmpˣ = Ω.rx .* ( Ω.Dʳ * ( vx .* θ )  + vx .* (Ω.Dʳ * θ) )
    tmpˣ += Ω.sx .* ( Ω.Dˢ * ( vx .* θ )  + vx .* (Ω.Dˢ * θ) )
    tmpʸ = Ω.ry .* ( Ω.Dʳ * ( vy .* θ )  + vy .* (Ω.Dʳ * θ) )
    tmpʸ += Ω.sy .* ( Ω.Dˢ * ( vy .* θ )  + vy .* (Ω.Dˢ * θ) )

    @. u⨀∇θ = (tmpˣ + tmpʸ) * 0.5
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
