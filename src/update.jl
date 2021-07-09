"""
update!(model::TVVAR)

Estimate a stochastic volatility model. Volatility is modelled as a random walk.

The hyperparameter 2/3 < β < 1 is a discount factor controlling the shocks to Σ.
"""

function update!(hypers::StochasticVol, m, P, S, Q, F, e)
    K = P * F / Q
    m = m + K * e'
    P = P - K * K' * Q
    S = S + e * e' / Q
    return (m, P, S)
end

function update!(hypers::ConstantVol, m, P, S, Q, F, e)
    K = P * F / Q
    m = m + K * e'
    P = P - K * K' * Q
    S = (hypers.ν * S + e * e' / Q) / (hypers.ν + 1)
    hypers.ν += 1
    return (m, P, S)
end