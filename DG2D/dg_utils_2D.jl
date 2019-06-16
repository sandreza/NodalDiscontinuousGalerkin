include("dg_utils.jl")

"""
rstoab(r,s)

Description:

Changes from (r,s) coordinate to (a,b) coordinates

Inputs: r,s

    r: a vector
    s: a vector

Return: a,b

    a: a vector: first coordinate pair
    b: a vector; second coordinate pair
"""
function rs2ab(r,s)
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
    return a,b
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

function ∇()
    return
end


"""
∇∘()

Description:

    Divergenence

"""

function ∇∘()
    return
end


"""
∇×()

Description:

    Curl

"""

function ∇×()
    return
end
