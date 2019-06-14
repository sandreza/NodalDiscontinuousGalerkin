dx = x[2] - x[1]

function upwind!(θ₁, θ₀, u, t)
    nx = length(θ₀)
    @inbounds for i in 2:(nx-1)
        θ₁[i] = ( θ₀[i + 1] - θ₀[i - 0] ) / (dx) * u[1]
    end
    θ₁[1] = (θ₀[1 + 1] - θ₀[1]) / (dx) * u[1]
    θ₁[end] = (θ₀[1] - θ₀[end]) / (dx) * u[1]
end

function upwind2!(θ₁, θ₀, u, t)
    nx = length(θ₀)
    @inbounds for i in 2:(nx-1)
        θ₁[i] = ( θ₀[i] - θ₀[i - 1] ) / (dx) * u[1]
    end
    θ₁[1] = (θ₀[1] - θ₀[end]) / (dx) * u[1]
    θ₁[end] = (θ₀[end] - θ₀[end-1]) / (dx) * u[1]
end

rhs! = upwind!
