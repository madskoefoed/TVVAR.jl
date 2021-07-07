"""
batch(model::TVVAR)

Estimate a stochastic volatility model. Volatility is modelled as a random walk.

The hyperparameter 2/3 < β < 1 is a discount factor controlling the shocks to Σ.
"""

function estimation(priors::Priors, y::Matrix{T}) where T <: AbstractFloat
    @unpack m, P, S, δ, β = priors

    @assert size(y, 2) == size(m, 2) "The second dimensiones of y and m do not match."

    N, p = size(y)
    n = 1/(1 - β)
    ν = n * β
    k = (β - p*β + p)/(2β - p*β + p - 1)
    d = size(m, 1)
    d = convert(Integer, (d - 1)/p)

    # Prepare model output
    model = Model(y,                 # y
                  expandArray(m, N), # m
                  expandArray(P, N), # P
                  expandArray(S, N), # S
                  δ,                 # δ
                  β,                 # β
                  repeat([ν], N),    # ν
                  zeros(T, N, p),    # μ
                  zeros(T, N, p, p)) # Σ

    for t in (d + 1):N

        if d == 0
            F = getF()
        else
            F = getF(y[t - d:t - 1, :])
        end

        # Predict
        m, P, S, Q, μ, Σ = predict(m, P, S, δ, β, k, F)

        # Prediction error
        e = y[t, :] - μ

        # Storage
        model.m[t, :, :] = m
        model.P[t, :, :] = P
        model.S[t, :, :] = S
        model.μ[t, :, :] = μ
        model.Σ[t, :, :] = Σ

        # Update
        m, P, S = update(m, P, S, Q, β, k, F, e)
        
    end
    return model
end

function predict(m, P, S, δ, β, k, F)
    P = P / δ
    Q = F' * P * F + 1.0
    S = S / k
    μ = m' * F
    Σ = Q * (1 - β)/(3β - 2) * S
    return (m, P, S, Q, μ, Σ)
end

function update(m, P, S, Q, β, k, F, e)
    K = P * F / Q
    m = m + K * e'
    P = P - K * K' * Q
    S = S + e * e' / Q
    return (m, P, S)
end