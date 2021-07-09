"""
estimate(mod)

Estimate a multivariate volatility model.

If δ < 1, the parameters are time-varying. If δ = 1, the model is a standard VAR model.
If β < 1, volatility is modelled as a random walk. If β = 1, the model assumes constant volatility.
"""

function estimate(priors::Priors, hypers::Hypers, y::Matrix{T}) where T <: AbstractFloat
    m, P, S = priors.m, priors.P, priors.S

    @assert size(y, 2) == size(m, 2) "The second dimensiones of y and m do not match."

    N, p = size(y)
    k = getk(p, hypers.β)
    d = size(m, 1)
    d = convert(Integer, (d - 1)/p)

    # Prepare model output
    model = Model(y[d+1:end, :],       # y
                  expandArray(m, N-d), # m
                  expandArray(P, N-d), # P
                  expandArray(S, N-d), # S
                  hypers.δ,            # δ
                  hypers.β,            # β
                  zeros(T, N-d),       # ν
                  zeros(T, N-d, p),    # μ
                  zeros(T, N-d, p, p), # Σ
                  zeros(T, N-d, p))    # e

    for t in (d + 1):N
        if d == 0
            F = getF()
        else
            F = getF(y[t - d:t - 1, :])
        end

        # Predict
        m, P, S, Q, μ, Σ = predict(hypers, m, P, S, k, F)

        # Prediction error
        e = y[t, :] - μ

        # Storage
        model.m[t-d, :, :] = m
        model.P[t-d, :, :] = P
        model.S[t-d, :, :] = S
        model.μ[t-d, :]    = μ
        model.Σ[t-d, :, :] = Σ
        model.e[t-d, :]    = e
        model.ν[t-d]       = hypers.ν
        
        # Update
        m, P, S = update!(hypers, m, P, S, Q, F, e)
        
    end
    return model
end