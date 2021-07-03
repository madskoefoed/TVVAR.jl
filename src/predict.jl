using Base: AbstractFloat

"""
estimate(model::MSV)

Estimate a stochastic volatility model. Volatility is modelled as a random walk.

The hyperparameter 2/3 < β < 1 is a discount factor controlling the shocks to Σ.
"""

# The kalman filter update
# m = m
# P = F' * P * F + 1.0

function predict!(prior::TVVAR, y::Vector{<:T}, x::Matrix{<:T}) where T <: AbstractFloat
    β = get_β(prior.ν)
    k = get_k(length(y), β)
    F = get_F(x)
    prior.P = prior.P / prior.δ
    Q = F' * prior.P * F + 1.0
    prior.S = prior.S / k
    prior.μ = prior.m' * F
    prior.Σ = Q * prior.S
end

function predict!(prior::VAR, y::Vector{<:T}, x::Matrix{<:T}) where T <: AbstractFloat
    F = get_F(x)
    prior.ν = prior.ν + 1.0
    prior.P = prior.P / prior.δ
    Q = prior.F' * prior.P * prior.F + 1.0
    prior.μ = prior.m' * prior.F
    prior.Σ = Q * (1 - β) / (3 - 2) * prior.S
end

function update!(prior::TVVAR, y::Vector{<:T}, x::Matrix{<:T}) where T <: AbstractFloat
    β = get_β(prior.ν)
    k = get_k(length(y), β)
    F = get_F(x)
    prior.P = prior.P / prior.δ
    Q = F' * prior.P * F + 1.0
    prior.S = prior.S / k
    prior.μ = prior.m' * F
    prior.Σ = Q * prior.S
end