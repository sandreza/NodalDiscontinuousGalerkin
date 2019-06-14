dx = x[2] - x[1]
#1/48, -(9/16), 9/16, -(1/48)
weights = [1/12; -(2/3); 2/3; -(1/12)] / (1 * dx )

let dx = x[2] - x[1]
global central4!
function central4!(θ₁, θ₀, u, t)
    nx = length(θ₀)
    for i in 3:(nx-2)
        θ₁[i]  = weights[1] * θ₀[i - 2] * u[1]
        θ₁[i] += weights[2]* θ₀[i - 1] * u[1]
        θ₁[i] += weights[3]* θ₀[i + 1] * u[1]
        θ₁[i] += weights[4]* θ₀[i + 2] * u[1]
    end
    θ₁[1] = (weights[1] * θ₀[end-1] + weights[2] * θ₀[end])* u[1]
    θ₁[1] +=  (weights[3] * θ₀[1 + 1] +  weights[4] * θ₀[1 + 2] ) * u[1]
    θ₁[2] = (weights[1] * θ₀[end] + weights[2] * θ₀[1])* u[1]
    θ₁[2] +=  (weights[3] * θ₀[2 + 1] +  weights[4] * θ₀[2 + 2] ) * u[1]
    θ₁[end] = (weights[1] * θ₀[end-2] + weights[2] * θ₀[end-1])* u[1]
    θ₁[end] +=  (weights[3] * θ₀[0 + 1] +  weights[4] * θ₀[0 + 2] ) * u[1]
    θ₁[end-1] = (weights[1] * θ₀[end-3] + weights[2] * θ₀[end-2])* u[1]
    θ₁[end-1] +=  (weights[3] * θ₀[end] +  weights[4] * θ₀[1] ) * u[1]
end
end

rhs! = central4!
