"""
update!(model::TVVAR)

Estimate a stochastic volatility model. Volatility is modelled as a random walk.

The hyperparameter 2/3 < β < 1 is a discount factor controlling the shocks to Σ.
"""

function update!(model::TVVAR, y::Vector{<:T}, F::Vector{<:T}) where T <: AbstractFloat
    @unpack m, P, S, δ, β, p, d, n, k, ν, Q, μ, Σ = model
    # Update step
    e = y - μ
    S = S + e * e' / Q
    K = P * F * Q
    m = m + K * e'
    P = P - K * K' / Q
    # Update model container
    model.m = m
    model.P = P
    model.S = S
end