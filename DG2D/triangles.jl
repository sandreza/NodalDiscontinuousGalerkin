include("../utils.jl")

"""
rstoab(r,s)

# Description

Changes from (r,s) coordinate to (a,b) coordinates.

# Inputs

- `r`: a vector
- `s`: a vector

# Return

- `a`: a vector: first coordinate pair
- `b`: a vector; second coordinate pair

# Example

a, b = rstoab([-0.5, 0.5], [0.5, -0.5])


"""
function rstoab(r,s)
    np = length(r);
    a = zeros(np)
    for i in 1:np
        if (s[i] ≈ 1)
            a[i] = -1
        else
            a[i] = 2*(1+r[i]) / (1-s[i]) -1
        end
    end
    b = s
    return a, b
end

function warp_factor(n, rout)
    LGLr = jacobiGL(0, 0, n)
    req = collect(range(-1,1, length = n+1))
    Veq = ones(n+1,n+1)
    jacobi!(Veq, req, 0, 0)
    nr = length(rout)
    pmat = zeros(n+1,nr)
    for i in 1:(n+1)
        pmat[i,:] = jacobi(rout, 0, 0, i-1)
    end
    return
end

"""
nodes2D(n)

Description:

    Returns the interpolation nodes for the 2D GL points

"""
function nodes2D(n)
    α° = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258];
    if n < 16
        α = α°[n]
    elseif
        α = 5/3
    end
    np = (n+1) * (n+2) / 2;
    L1 = zeros(np)
    L2 = zeros(np)
    L3 = zeros(np)
    sk = 1;
    for j ∈ 1:(n+1)
        for m ∈ 1:(n+2-j)
            L1[sk] = (j-1)/n
            L3[sk] = (m-1)/n
            sk += 1
        end
    end
    @. L2 = 1.0 - L1 - L3
    x = -L2 + L3
    y = (-L2 - L3 + 2 .* L1 ) / sqrt(3)
    blend1 = 4 * L2 .* L3
    blend2 = 4 * L1 .* L3
    blend3 = 4 * L1 .* L3
    warpf1 = warp_factor(n, L3 - L2)
    warpf2 = warp_factor(n, L1 - L3)
    warpf3 = warp_factor(n, L2 - L1)
    warp1 = @. blend1 * warpf1 * (1 + (α * L1 )^2 )
    warp2 = @. blend2 * warpf2 * (1 + (α * L2 )^2 )
    warp3 = @. blend3 * warpf3 * (1 + (α * L3 )^2 )
    x = x + 1*warp1 + cos(2π/3) * warp2 + cos(4π/3) * warp3
    y = y + 0*warp1 + sin(2π/3) * warp2 + sin(4π/3) * warp3
    return x, y
end


"""
xy2rs(n)

Description:

    Returns the rs coordinates

"""
function xy2rs(x,y)
    return
end

"""
vandermonde2D(n,r,s)

Description:

    Matrix to convert from modal to nodal basis

"""
function vandermonde2D(n,r,s)
    return
end
