"""
    IRKGLInterpolation

Dense output of an [`IRKGL16`](@ref) solution: on each step `[t[i], t[i+1]]` the
solution is evaluated with the collocation polynomial of that step,
`u[i] + sum_j kappa_j(θ) * L[i][j, :]`, where `θ = (t - t[i]) / h[i]` and `h[i]`
is the signed length of the step — the same polynomial the `saveat` path uses,
except for the position block of `second_order_ode` solutions described below.
At the saved time points the saved values are returned.

For `second_order_ode` solutions (state `[q; v]`, `q'' = f(u)`) the steppers
only converge the velocity half of `L`, so the position block is instead
evaluated with the Nyström interpolant
`q(t) = q[i] + (t - t[i]) * v[i] + h[i] * sum_j κ̃_j(θ) * L[i][j, v-block]`,
where `κ̃_j` interpolates the double-integral coefficients `eta` on the nodes
`[0, c, 1]` — including the endpoint value `κ̃_j(1) = 1 - c_j`, which makes the
interpolant agree with the saved endpoint values. `κ̃` is not the exact double
integral of the acceleration polynomial (they differ in the highest-degree
term), so `Val{1}` differentiates `κ̃` itself:
`q'(t) = v[i] + sum_j κ̃'_j(θ) * L[i][j, v-block]`.

`Val{0}` interpolates the solution and `Val{1}` its first derivative; higher
derivative orders throw an `ArgumentError`. A query time carrying AD
partials (e.g. `ForwardDiff.Dual`) propagates through the polynomial, and at a
saved time the step is selected on the primal value with the same `continuity`
rule as a `Float64` query.
"""
struct IRKGLInterpolation{tType, uType, LType} <: SciMLBase.AbstractDiffEqInterpolation
    t::Vector{tType}
    u::Vector{uType}
    h::Vector{tType}
    L::Vector{LType}
    X2::Vector{tType}
    Y2::Matrix{tType}
    X3::Vector{tType}
    Y3::Matrix{tType}
    XD::Vector{tType}
    YD::Matrix{tType}
    lenq::Int
end

SciMLBase.interp_summary(::IRKGLInterpolation) = "Gauss-Legendre collocation polynomial"

function _store_step_L!(Ls, L, indices, s, use_simd)
    M = Matrix{eltype(eltype(Ls))}(undef, s, length(indices))
    for (n, k) in enumerate(indices)
        if use_simd
            Lk = L[k]
            for j in 1:s
                M[j, n] = Lk[j]
            end
        else
            for j in 1:s
                M[j, n] = L[j][k]
            end
        end
    end
    push!(Ls, M)
    return nothing
end

function _irkgl_interp_point(
        id::IRKGLInterpolation, tval, idxs, ::Type{Val{D}}, continuity::Symbol
    ) where {D}
    if D > 1
        throw(
            ArgumentError(
                "IRKGL16 dense output only supports interpolating the solution " *
                    "(Val{0}) and its first derivative (Val{1}), got Val{$D}."
            )
        )
    end
    t = id.t
    tdir = sign(t[end] - t[1])
    # Bounds checks, interval lookup, and knot equality run on the primal value
    # of the query time: dual-numbered queries compare on their partials too,
    # which would sort a knot query past the knot and bypass the `continuity`
    # branch. The polynomial itself still sees `tval`, so derivatives propagate.
    tv = SciMLBase.unitfulvalue(tval)
    if tdir * tv > tdir * t[end] || tdir * tv < tdir * t[1]
        error("Solution interpolation cannot extrapolate outside of the time span [$(t[1]), $(t[end])].")
    end
    i = searchsortedfirst(t, tv; rev = tdir < 0)
    if D == 0 && tval isa AbstractFloat && t[i] == tv
        u = id.u[i]
        return idxs === nothing ? copy(u) : u[idxs]
    end
    step = if t[i] == tv
        i == firstindex(t) ? i : (continuity === :right && i < lastindex(t) ? i : i - 1)
    else
        i - 1
    end
    u0 = id.u[step]
    Lstep = id.L[step]
    h = id.h[step]
    θ = (tval - t[step]) / h
    s = size(Lstep, 1)
    lenq = id.lenq
    u = similar(u0, promote_type(eltype(u0), typeof(θ)))
    if D == 0
        κv = vec(PolInterp(id.X2, id.Y2, [θ]))
        if lenq > 0
            κq = vec(PolInterp(id.X3, id.Y3, [θ]))
            for n in 1:lenq
                nv = n + lenq
                acc = κq[1] * Lstep[1, nv]
                for j in 2:s
                    acc = muladd(κq[j], Lstep[j, nv], acc)
                end
                u[n] = u0[n] + (tval - t[step]) * u0[nv] + h * acc
            end
        end
        for (n, k) in enumerate(eachindex(u0))
            lenq > 0 && n <= lenq && continue
            acc = κv[1] * Lstep[1, n]
            for j in 2:s
                acc = muladd(κv[j], Lstep[j, n], acc)
            end
            u[k] = u0[k] + acc
        end
    else
        κd = vec(PolInterp(id.XD, id.YD, [θ]))
        if lenq > 0
            κqd = vec(PolInterpDer(id.X3, id.Y3, [θ]))
            for n in 1:lenq
                nv = n + lenq
                acc = κqd[1] * Lstep[1, nv]
                for j in 2:s
                    acc = muladd(κqd[j], Lstep[j, nv], acc)
                end
                u[n] = u0[nv] + acc
            end
        end
        for (n, k) in enumerate(eachindex(u0))
            lenq > 0 && n <= lenq && continue
            acc = κd[1] * Lstep[1, n]
            for j in 2:s
                acc = muladd(κd[j], Lstep[j, n], acc)
            end
            u[k] = acc / h
        end
    end
    return idxs === nothing ? u : u[idxs]
end

_irkgl_interp(id, tval::Number, idxs, deriv, continuity) =
    _irkgl_interp_point(id, tval, idxs, deriv, continuity)
function _irkgl_interp(id, tvals, idxs, deriv, continuity)
    return DiffEqArray(
        [_irkgl_interp_point(id, tv, idxs, deriv, continuity) for tv in tvals],
        tvals
    )
end

function _irkgl_interp!(val, id, tval::Number, idxs, deriv, continuity)
    val .= _irkgl_interp_point(id, tval, idxs, deriv, continuity)
    return val
end
function _irkgl_interp!(vals, id, tvals, idxs, deriv, continuity)
    for (j, tv) in enumerate(tvals)
        vals[j] = _irkgl_interp_point(id, tv, idxs, deriv, continuity)
    end
    return vals
end

function (id::IRKGLInterpolation)(
        tvals, idxs, deriv::D, p, continuity::Symbol = :left
    ) where {D}
    return _irkgl_interp(id, tvals, idxs, deriv, continuity)
end

function (id::IRKGLInterpolation)(
        val, tvals, idxs, deriv::D, p, continuity::Symbol = :left
    ) where {D}
    return _irkgl_interp!(val, id, tvals, idxs, deriv, continuity)
end
