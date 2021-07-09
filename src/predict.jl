"""
estimate(model::TVVAR)

Estimate a stochastic volatility model. Volatility is modelled as a random walk.

The hyperparameter 2/3 < β < 1 is a discount factor controlling the shocks to Σ.
"""

function predict(hypers::StochasticVol, m, P, S, k, F)
    P = P / hypers.δ
    Q = F' * P * F + 1.0
    S = S / k
    μ = m' * F
    Σ = Q * (1 - hypers.β)/(3*hypers.β - 2) * S
    return (m, P, S, Q, μ, Σ)
end

function predict(hypers::ConstantVol, m, P, S, k, F)
    P = P / hypers.δ
    Q = F' * P * F + 1.0
    μ = m' * F
    Σ = Q * S
    return (m, P, S, Q, μ, Σ)
end