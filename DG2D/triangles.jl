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
    veq = vandermonde(req, 0, 0, n)
    nr = length(rout)
    pmat = vandermonde(rout, 0, 0, n)'
    lmat = veq' \ pmat
    warp = lmat' * (LGLr - req)
    zerof = (abs.(rout) .< (1.0 -  1.0e-10) )
    sf = @. 1.0 - (zerof *rout )^2
    warp = warp ./ sf + warp .* (zerof .- 1)
    return warp
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
    else
        α = 5/3
    end
    np = Int((n+1) * (n+2) / 2);
    L1 = zeros(np)
    L2 = zeros(np)
    L3 = zeros(np)
    let sk = 1
        for j ∈ 1:(n+1)
            for m ∈ 1:(n+2-j)
                L1[sk] = (j-1)/n
                L3[sk] = (m-1)/n
                sk += 1
            end
        end
    end
    @. L2 = 1.0 - L1 - L3
    x = -L2 + L3
    y = (-L2 - L3 + 2 .* L1 ) / sqrt(3)
    blend1 = 4 * L2 .* L3
    blend2 = 4 * L3 .* L1
    blend3 = 4 * L1 .* L2
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
xytors(n)

# Description

-    maps the (x,y) coordinates of an equilateral triangle to the standard right triangle

# Arguments

-    `x`: x-coordinate values in the equilateral triangle, an array
-    `y`: y-coordinate values in the equilateral triangle, an array

# Outputs: r,s

-    `r`: r-coordinate values in the standard right triangle, an array
-    `s`: s-coordinate values in the standard right triangle, an array

"""
function xytors(x,y)
    L1 = @. (sqrt(3.0) * y + 1) / 3.0
    L2 = @. (-3.0 * x - sqrt(3.0) * y + 2.0) / 6.0
    L3 = @. ( 3.0 * x - sqrt(3.0) * y + 2.0) / 6.0
    r = -L2 + L3 - L1
    s = -L2 -L3 + L1
    return r, s
end

"""
simplex2DP

# Description

- Evaluate the Jacobi polynomials at nodal values on the a,b grid

# Inputs

- `a`: first coordinate
- `b`: second coordinate
- `i`: jacobi polynomial parameter, Legendre hard coded
- `j`: jacobi polynomial parameter, Legendre hard coded

# Outputs

- `p`: value of the jacobi polynomial

"""
function simplex2DP(a,b, i::Int, j::Int)
    h1 = @. jacobi(a, 0, 0, i)
    h2 = @. jacobi(b, 2*i + 1, 0, j)
    p = @. sqrt(2.0) * h1 * h2 * (1-b)^i
    return p
end


"""
vandermonde2D(n,r,s)

# Description

    Matrix to convert from modal to nodal basis. Depends on Simplex2DP

# Arguments

- `n`:
- `r`:
- `s`:



"""
function vandermonde2D(n,r,s)
    np = Int((n+1) * (n+2) / 2);
    V2D = zeros(length(r), np)
    a,b = rstoab(r,s)
    let sk = 1
        for i ∈ 0:n
            for j ∈ 0:(n - i)
                V2D[:, sk] = simplex2DP(a,b,i,j)
                sk += 1
            end
        end
    end
    return V2D
end




"""
∇()

Description:

    gradient

"""
function ∇(u)
    return nothing
end


"""
∇o()

Description:

    Divergenence

"""
function ∇o(u)
    return nothing
end


"""
∇x()

# Description

- The curl operator in 2D. unicode is (nabla x)

    Curl

"""
function ∇x(u)
    return nothing
end
