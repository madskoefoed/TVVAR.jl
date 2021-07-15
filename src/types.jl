##############
### Hypers ###
##############

abstract type Hypers end

mutable struct StochasticVol{T<:AbstractFloat} <: Hypers
    δ::T
    β::T
    n::T
    ν::T
    function StochasticVol(δ::T, β::T) where T <: AbstractFloat
        (δ >   0 && δ ≤ 1) || throw(ArgumentError("0 < δ ≤ 1 required (currently $δ)."))
        (β ≥ 2/3 && β < 1) || throw(ArgumentError("2/3 ≤ β < 1 required (currently $β)."))
        n = 1/(1 - β)
        ν = n * β
        new{T}(δ, β, n, ν)
    end
end

mutable struct ConstantVol{T<:AbstractFloat} <: Hypers
    δ::T
    β::T
    n::T
    ν::T
    function ConstantVol(δ::T, n::T) where T <: AbstractFloat
        (δ > 0 && δ ≤ 1) || throw(ArgumentError("0 < δ ≤ 1 required (currently $δ)."))
        n ≥ 3 || throw(ArgumentError("n ≥ 3 required (currently $n)."))
        β = 1.0
        ν = n * β
        new{T}(δ, β, n, ν)
    end
end

##############
### Priors ###
##############

mutable struct Priors{T<:AbstractFloat}
    m::Array{T, 2}
    P::Array{T, 2}
    S::Array{T, 2}

    function Priors(m::Array{T, 2},
                    P::Array{T, 2},
                    S::Array{T, 2}) where T <: AbstractFloat#,

        str = "m is $(size(m, 1)) x $(size(m, 2))," *
              "\nP is $(size(P, 1)) x $(size(P, 2))," *
              "\nS is $(size(S, 1)) x $(size(S, 2))."

        !(size(m, 1) == size(P, 1) == size(P, 2)) && throw(DimensionMismatch(str))
        !(size(m, 2) == size(S, 1) == size(S, 2)) && throw(DimensionMismatch(str))

        any(diag(S) .<= 0) && throw(ArgumentError("The diagonal elements of S must be strictly positive."))
        any(diag(P) .<= 0) && throw(ArgumentError("The diagonal elements of P must be strictly positive."))

        new{T}(m, P, S)
    end
end

# Outer constructor for univariate case (p = 1)
function Priors(m::Array{T, 1}, P::Array{T, 2}, S::T) where T <: AbstractFloat
    m = reshape(m, 1, 1)
    S = reshape([S], 1, 1)
    return (Priors(m, P, S))
end

# Outer constructor for univariate stochastic volatility case (p = 1, d = 0)
function Priors(m::T, P::T, S::T) where T <: AbstractFloat
    m = reshape([m], 1, 1)
    P = reshape([P], 1, 1)
    S = reshape([S], 1, 1)
    return (Priors(m, P, S))
end

# Outer constructor for multivariate stochastic volatility case (d = 0)
function Priors(m::T, P::T, S::Array{T, 2}) where T <: AbstractFloat
    m = reshape([m], 1, 1)
    P = reshape([p], 1, 1)
    return (Priors(m, P, S))
end

##################
### Simulation ###
##################

mutable struct Simulation{T<:AbstractFloat}
    m::Array{T, 3}
    Σ::Array{T, 3}

    function Simulation(m::Array{T, 3}, Σ::Array{T, 3}) where T <: AbstractFloat

        str = "m is $(size(m, 1)) x $(size(m, 2)) x $(size(m, 3))," *
              "\nΣ is $(size(Σ, 1)) x $(size(Σ, 2)) x $(size(Σ, 3))."

        !(size(m, 1) == size(Σ, 1)) && throw(DimensionMismatch(str))
        !(size(m, 3) == size(Σ, 2) == size(Σ, 3)) && throw(DimensionMismatch(str))

        #any(diag(Σ) .<= 0) && throw(ArgumentError("The diagonal elements of Σ must be strictly positive."))

        new{T}(m, Σ)
    end
end

# Outer constructor
function Simulation(m::Array{T, 2}, Σ::Array{T, 2}, N::Integer) where T <: AbstractFloat
    Simulation(expandArray(m, N), expandArray(Σ, N))
end

#############
### Model ###
#############

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
    e::Array{T, 2}
end