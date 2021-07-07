"""
step(model::TVVAR)

Estimate a stochastic volatility model. Volatility is modelled as a random walk.

The hyperparameter 2/3 < β < 1 is a discount factor controlling the shocks to Σ.
"""

function step(model::TVVAR, y::Vector{<:T}, x::Matrix{<:T}) where T <: AbstractFloat
    # Vectorize x
    F = get_F(x)
    # Predict step
    predict!(model, F)
    # Update step
    update!(model, y, F)
end