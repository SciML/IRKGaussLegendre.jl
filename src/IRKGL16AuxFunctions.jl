#
#   Auxiliarfunctions.jl file:
#       ErrorEst:
#       MyNorm
#       Rdigits

function ErrorEst(U, F, dt, beta, abstol, reltol)
    uiType = eltype(U[1])

    (s,) = size(F)
    D = length(U[1])

    est = zero(uiType)

    @inbounds begin for k in eachindex(U[1])
        sum = zero(uiType)
        maxU = zero(uiType)
        for is in 1:s
            sum += beta[is] * F[is][k]
            maxU = max(maxU, abs(U[is][k]))
        end

        est += (abs(dt * sum))^2 / (abstol + maxU^2 * reltol)
    end end

    return (est / D)
end

function MyNorm(u, abstol, reltol)
    uiType = eltype(u)
    norm = zero(uiType)

    @inbounds begin for k in eachindex(u)
        aux = u[k] / (abstol + abs(u[k]) * reltol)
        norm += aux * aux
    end end

    norm = sqrt(norm / (length(u)))

    return (norm)
end

function Rdigits(x::Real, r::Integer)
    mx = r * x
    mxx = mx + x
    return (mxx - mx)
end
