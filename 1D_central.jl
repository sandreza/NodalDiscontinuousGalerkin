let dx = x[2] - x[1]
global central!
function central!(θ₁, θ₀, u, t)
    nx = length(θ₀)
    @inbounds for i in 2:(nx-1)
        θ₁[i] = ( θ₀[i + 1] - θ₀[i - 1] ) / (2*dx) * u[1]
    end
    θ₁[1] = (θ₀[1 + 1] - θ₀[end]) / (2*dx) * u[1]
    θ₁[end] = (θ₀[1] - θ₀[end-1]) / (2*dx) * u[1]
end
end


rhs! = central!
