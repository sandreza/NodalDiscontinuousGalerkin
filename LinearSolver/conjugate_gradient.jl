using LinearAlgebra

"""
conjugate_gradient!(L, x⁰, b; P = x->x, tolerance = eps(1.0), maximum_iterations = false, track_residual = false)

# Description

- An implementation of conjugate gradient. Solves Lx = b

# Arguments
- 'L': (function), a linear operator
- 'x⁰': (vector), initial condition, [OVERWRITTEN]
- 'b': (vector), right-hand side

# Keyword Arguments

- 'P': (function). This is the preconditioner. The default is the identity function x->x.
- 'tolerance': (number). relative error tolerance. default = eps(1.0)
- 'maximum_iterations': (integer). default = length(x⁰), Maximum iterations for conjugate gradient
- 'track residual': (boolean). keeps track of norm of residual and returns it

# Return
- Nothing if track_residual is false. the norm of the residual if track_residual is true

"""
function conjugate_gradient!(L, x⁰, b; P = x->x, tolerance = eps(1.0), maximum_iterations = length(x⁰), track_residual = false)
    # calculate the residual and auxillary field
    r⁰ = b - L(x⁰)
    z⁰ = P(r⁰)
    p⁰ = copy(z⁰)
    # check to see if the guess was fantastic
    if tolerance_boolean(r⁰, b, tolerance)
        return nothing
    end
    if track_residual
        r_tracked = copy(r⁰)
    end

    # start searching
    for j in 1:maximum_iterations
        # create search step size
        Lp = L(p⁰)
        α = (r⁰' * z⁰) / (p⁰' * Lp)
        # update along preconditioned direction
        @. x⁰ += α .* p⁰
        # form new residual
        r¹ = r⁰ - α .* Lp
        # track it
        if track_residual
            r_tracked[j] = norm(r¹)
        end
        # check to see if the update was reasonable

        if tolerance_boolean(r¹, b, tolerance)
            return nothing
        end
        # update p⁰
        z¹ = P(r¹)
        β  = (z¹' * r¹) / (z⁰' * r⁰)
        @. p⁰ = z¹ + β .* p⁰
        # rinse repeate
        @. z⁰ = z¹
    end
    if track_residual
        return r_tracked
    end
    return nothing
end


"""
tolerance_boolean(r, b, tolerance)

# Description
- Calculates tolerance

# Arguments
- 'r': (vector), the residual
- 'b': (vector), takes into account relative size
- 'tolerance': checks the tolerance

# Return
- boolean value: true or false
"""
function tolerance_boolean(r, b, tolerance)
    return norm(r)/norm(b) < tolerance
end
