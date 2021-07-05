using Base: AbstractFloat

"""
estimate(model::MSV)

Estimate a stochastic volatility model. Volatility is modelled as a random walk.

The hyperparameter 2/3 < β < 1 is a discount factor controlling the shocks to Σ.
"""

function predict!(model::TVVAR)
    @unpack y, x, m, P, S, δ, ν = prior

    F = get_F(x)
    β = get_β(ν)
    k = get_k(length(y), β)
    P = P / δ
    Q = F' * P * F + 1.0
    μ = m' * F
    Σ = Q * (1 - β) / (3 - 2) * S
    e = y .- μ

    # Update
    model.μ = μ
    model.Σ = Σ
    model.e = e
    model.Q = Q
end

function update(model::TVVAR)

end

function timestep(prior::Prior)
    # Unpack
    m, P, S, δ, ν = posterior

    F = get_F(x)
    β = get_β(tvvar.ν)
    k = get_k(length(y), β)
    m = tvvar.m
    P = tvvar.P / tvvar.δ
    Q = F' * P * F + 1.0
    μ = m' * F
    Σ = Q * (1 - β) / (3 - 2) * S
    e = y .- μ
    S = S / k + e * e' / Q
    K = P * F * Q
    m = tvvar.m + K * e
    tvvar = TVVAR(y, x, m, P, S, δ, ν, μ, Σ)
end