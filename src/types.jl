using Base: AbstractFloat
using Parameters: @unpack
using Distributions: MultivariateNormal, Wishart, InverseWishart

mutable struct TVVAR{T<:AbstractFloat}
    y::Matrix{<:T}
    Φ::MultivariateNormal
    Σ⁻¹::Wishart

    function TVVAR(y::Matrix{<:T},
        Φ::MultivariateNormal
    Σ⁻¹::Wishart
                   ) where T <: AbstractFloat

        str = "y is $(length(y))," *
              "\nm is $(size(x, 1)) x $(size(x, 2))," *
              "\nm is $(size(m, 1)) x $(size(m, 2))," *
              "\nP is $(size(P, 1)) x $(size(P, 2))," *
              "\nS is $(size(S, 1)) x $(size(S, 2))."

        !(size(m, 1) == size(P, 1) == size(P, 2)) && throw(DimensionMismatch(str))
        !(size(m, 2) == size(S, 1) == size(S, 2)) && throw(DimensionMismatch(str))

        any(diag(S) .<= 0) && throw(ArgumentError("The diagonal elements of S must be strictly positive."))
        any(diag(P) .<= 0) && throw(ArgumentError("The diagonal elements of P must be strictly positive."))

        (δ > 0 && δ ≤ 1) || throw(ArgumentError("0 < δ ≤ 1 required (currently $δ)."))
        ν ≥ 2 || throw(ArgumentError("ν ≥ 2 required (currently $ν)."))

        return new{T}(y, x, m, P, S, δ, ν)
    end
end