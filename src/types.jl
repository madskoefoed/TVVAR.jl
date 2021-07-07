mutable struct Priors{T<:AbstractFloat}
    m::Array{T, 2}
    P::Array{T, 2}
    S::Array{T, 2}
    δ::T
    β::T

    function Priors(m::Array{T, 2},
                    P::Array{T, 2},
                    S::Array{T, 2},
                    δ::T,
                    β::T) where T <: AbstractFloat

        str = "m is $(size(m, 1)) x $(size(m, 2))," *
              "\nP is $(size(P, 1)) x $(size(P, 2))," *
              "\nS is $(size(S, 1)) x $(size(S, 2))."

        !(size(m, 1) == size(P, 1) == size(P, 2)) && throw(DimensionMismatch(str))
        !(size(y, 2) == size(m, 2) == size(S, 1) == size(S, 2)) && throw(DimensionMismatch(str))

        any(diag(S) .<= 0) && throw(ArgumentError("The diagonal elements of S must be strictly positive."))
        any(diag(P) .<= 0) && throw(ArgumentError("The diagonal elements of P must be strictly positive."))

        (δ >   0 && δ ≤ 1) || throw(ArgumentError("0 < δ ≤ 1 required (currently $δ)."))
        (β ≥ 2/3 && β < 1) || throw(ArgumentError("2/3 ≤ β < 1 required (currently $β)."))

        new{T}(m, P, S, δ, β)
    end
end

mutable struct Model{T<:AbstractFloat}
    y::Array{T, 2}
    m::Array{T, 3}
    P::Array{T, 3}
    S::Array{T, 3}
    δ::T
    β::T
    ν::Array{T, 1}
    μ::Array{T, 2}
    Σ::Array{T, 3}
end