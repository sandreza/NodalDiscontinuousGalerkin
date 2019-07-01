include("../DG1D/dg_utils.jl")

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
    return
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
