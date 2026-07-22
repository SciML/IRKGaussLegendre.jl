struct CompensatedSum{T <: Number}
    hi::T
    lo::T
end

function compensated_sum(x::Number, y::Number)
    x, y = promote(x, y)
    hi = x + y
    z = hi - x
    return CompensatedSum(hi, (x - (hi - z)) + (y - z))
end

function Base.:+(x::CompensatedSum, y::Number)
    sum = compensated_sum(x.hi, y)
    return compensated_sum(sum.hi, sum.lo + x.lo)
end

Base.:+(x::Number, y::CompensatedSum) = y + x

function Base.:+(x::CompensatedSum, y::CompensatedSum)
    sum = compensated_sum(x.hi, y.hi)
    return compensated_sum(sum.hi, sum.lo + x.lo + y.lo)
end

Base.:-(x::CompensatedSum, y::Number) = x + (-y)

function Base.:*(x::Number, y::CompensatedSum)
    return compensated_sum(x * y.hi, x * y.lo)
end

Base.:*(x::CompensatedSum, y::Number) = y * x

normalize_verbose(verbose::AbstractVerbosityPreset) = DEVerbosity(verbose)
normalize_verbose(verbose::DEVerbosity) = verbose
