include("field2D.jl")

abstract type AbstractFlux2D end

"""
Flux2D()

# Description

# Arguments

"""
struct Flux2D{S, T} <: AbstractFlux2D
    # inputs
    fields::S
    weights::T

    function Flux2D(fields, weights)
        if length(fields) == length(weights)
            return new{typeof(fields),typeof(weights)}(fields, weights)
        else
            println("Please pass equal number of weights and fields.")
            return nothing
        end
    end
end

function computePhysicalFlux!(φ, flux::Flux2D, Ω::Element2D)
    @. φ[Ω.iⱽ] = 0.0
    for (𝑓, 𝑤) in zip(flux.fields, flux.weights)
        @. φ[Ω.iⱽ] += 𝑤 * 𝑓.ϕ[Ω.iⱽ]
    end
end

function computeNumericalFlux!(φ, flux::Flux2D, f::Face2D)
    @. φ[f.i⁻] = 0.0
    for (𝑓, 𝑤) in zip(flux.fields, flux.weights)
        @. φ[f.i⁻] += 𝑤 * 𝑓.ϕ°[f.i⁻]
    end
end
