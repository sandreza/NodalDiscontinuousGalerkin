#θ₀ should be defined

 ik = similar(x * im)
for i in 1:(Int(n/2)+1)
    ik[i] = (i-1)* (2* pi /L) * im
end
for i in 0:(Int(n/2)-2)
    ik[n-i] = - (i+1) * (2* pi /L) * im
end

θ₀ = exp.(-4*(x.-L/2).^2) .+ im*0
#θ₀ = cos.(x) .+ im*0
θ̂₀ = similar(θ₀)
θ₁ = similar(θ₀)

P = FFTW.plan_fft(x*im)
iP = FFTW.plan_ifft(x*im)
mul!(θ̂₀,P,θ₀ )

#basically takes the derivative
function spectral!(θ₁, θ₀, u, t)
    mul!(θ₁,P,θ₀)
    tmp = ik .* θ₁ * u[1]
    mul!(θ₁,iP,tmp)
end

function spectral2!(θ₁, θ₀, u, t)
    θ₁ .= θ₀
    fft!(θ₁)
    θ₁ .= ik .* θ₁ * u[1]
    ifft!(θ₁)
    return nothing
end

function spectral!(θ₁, θ₀, u, t, ik)
    mul!(θ₁,P,θ₀)
    tmp = ik .* θ₁ * u[1]
    mul!(θ₁,iP,tmp)
end

rθ = similar(x)
rθ .= real(θ₀)
plt = plot(x, rθ, marker=4,title="initial condition");
display(plt)
plt = plot(abs.(θ̂₀),marker=4, title="fourier amplitudes")
display(plt)

spectral!(θ₁,θ₀,[1.0],1.0)
plt = plot(x, real.(θ₁), marker=4,title="derivative of initial condition");
display(plt)

rhs! = spectral!
