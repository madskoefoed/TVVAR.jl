using Base: AbstractFloat

abstract type StateSpaceModel end

mutable struct TVVAR{T<:AbstractFloat} <: StateSpaceModel
    m::Matrix{<:T}
    P::Matrix{<:T}
    S::Matrix{<:T}
    δ::T
    ν::T
    μ::Vector{<:T}
    Σ::Matrix{<:T}

    function TVVAR(m::Matrix{<:T}, P::Matrix{<:T}, S::Matrix{<:T}, δ::T, ν::T) where T <: AbstractFloat
        str = "nm is $(size(m, 1)) x $(size(m, 2)),\nP is $(size(P, 1)) x $(size(P, 2)), \nS is $(size(S, 1)) x $(size(S, 2))."

        !(size(m, 1) == size(P, 1) == size(P, 2)) && throw(DimensionMismatch(str))
        !(size(m, 2) == size(S, 1) == size(S, 2)) && throw(DimensionMismatch(str))

        any(diag(S) .<= 0) && throw(ArgumentError("The diagonal elements of S must be strictly positive."))
        any(diag(P) .<= 0) && throw(ArgumentError("The diagonal elements of P must be strictly positive."))

        (δ > 0 && δ ≤ 1) || throw(ArgumentError("0 < δ ≤ 1 required (currently $δ)."))
        ν ≥ 2 || throw(ArgumentError("2 ≤ ν required (currently $ν)."))

        return new{T}(m, P, S, δ, ν)
    end
end

mutable struct VAR{T<:AbstractFloat} <: StateSpaceModel
    m::Matrix{<:T}
    P::Matrix{<:T}
    S::Matrix{<:T}
    δ::T
    ν::T
    μ::Vector{<:T}
    Σ::Matrix{<:T}

    function VAR(m::Matrix{<:T}, P::Matrix{<:T}, S::Matrix{<:T}, δ::T, ν::T) where T <: AbstractFloat
        str = "nm is $(size(m, 1)) x $(size(m, 2)),\nP is $(size(P, 1)) x $(size(P, 2)), \nS is $(size(S, 1)) x $(size(S, 2))."

        !(size(m, 1) == size(P, 1) == size(P, 2)) && throw(DimensionMismatch(str))
        !(size(m, 2) == size(S, 1) == size(S, 2)) && throw(DimensionMismatch(str))

        any(diag(S) .<= 0) && throw(ArgumentError("The diagonal elements of S must be strictly positive."))
        any(diag(P) .<= 0) && throw(ArgumentError("The diagonal elements of P must be strictly positive."))

        (δ > 0 && δ ≤ 1) || throw(ArgumentError("0 < δ ≤ 1 required (currently $δ)."))
        ν > 3 || throw(ArgumentError("ν ≤ 3 required (currently $ν)."))

        return new{T}(m, P, S, δ, ν)
    end
end