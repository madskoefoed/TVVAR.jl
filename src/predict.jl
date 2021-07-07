"""
estimate(model::TVVAR)

Estimate a stochastic volatility model. Volatility is modelled as a random walk.

The hyperparameter 2/3 < β < 1 is a discount factor controlling the shocks to Σ.
"""

function predict!(model::TVVAR, F::Vector{<:T}) where T <: AbstractFloat
    @unpack m, P, S, δ, β, p, d, n, k, ν, Q = model
    # Predict step
    P = P / δ
    S = S / k
    Q = F' * P * F + 1.0
    μ = m' * F
    Σ = Q * (1 - β) / (3 - 2) * S
    # Update model container
    model.P = P
    model.S = S
    model.Q = Q
    model.μ = μ
    model.Σ = Σ
end