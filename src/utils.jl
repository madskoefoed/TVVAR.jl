num2mat(x::Real) = fill(x, 1, 1)
vec2mat(x::Vector) = repeat(x, 1, 1)

expandArray(x, N) = repeat(reshape(x, 1, size(x, 1), size(x, 2)), N, 1, 1)

function getF()
    F = [1.0]
    return F
end

function getF(x::Matrix{<:AbstractFloat})
    p, d = size(x)
    F = [1.0]
    append!(F, vec(x))
    return F
end