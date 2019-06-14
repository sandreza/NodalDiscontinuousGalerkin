dx = x[2] - x[1]
#1/48, -(9/16), 9/16, -(1/48)
dx = x[2] - x[1]
weights = [-(1/60); 3/20; -(3/4); 0; 3/4; -(3/20); 1/60] / (1 * dx )
function central6!(θ₁, θ₀, u, t)
    nx = length(θ₀)
    for i in eachindex(θ₀)
        if mod(i - 3,nx)==0
            ii = nx
        else
            ii = mod(i - 3,nx)
        end
        θ₁[i]  = weights[1] * θ₀[ii] * u[1]
        for j in 2:7
            if mod(i + j - 4,nx)==0
                ii = nx
                else
                ii = mod(i + j - 4,nx)
                end
            θ₁[i] += weights[j]* θ₀[ii] * u[1]
        end
    end
end
