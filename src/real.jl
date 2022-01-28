# making FUSE work with Reals that are non-Floats requires defining some functions

# =============== #
#  Measurements.jl
# =============== #

import Ratios

function Int(x::Measurement)
    return Int(x.val)
end

function Base.convert(::Type{Measurement{T}}, x::Ratios.SimpleRatio{S}) where {T<:AbstractFloat, S}
    return x.num/x.den ± 0.0
end

function Base.unsafe_trunc(::Type{Int64}, x::Measurement{Float64})
    return Int(x.val)
end
