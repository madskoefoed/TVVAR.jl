using Base: AbstractFloat

num2mat(x::Real) = fill(x, 1, 1)
vec2mat(x::Vector) = repeat(x, 1, 1)

get_ν(β) = β/(1 - β)
get_β(ν) = ν/(1 + ν)
get_k(p, β) = (β - p*β + p)/(2β - p*β + p - 1)

function get_F()
    F = [1.0]
    return F
end
function get_F(x::Matrix{<:AbstractFloat})
    p, d = size(x)
    F = [1.0]
    append!(F, vec(x))
    return F
end