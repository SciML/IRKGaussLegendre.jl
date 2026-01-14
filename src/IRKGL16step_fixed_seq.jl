#
#  IRKstep_fixed!
#  IRKNGLstep_fixed_2nd!

function IRKstep_fixed!(
        ttj::Array{tType, 1},
        uj::uType,
        ej::uType,
        dts::Array{tType, 1},
        stats::SciMLBase.DEStats,
        coeffs::tcoeffs{tType},
        cache::tcache{uType, realuType, tType, fT, pT}
    ) where {
        uType, realuType, tType, fT, pT,
    }
    @unpack mu, c, b, nu, alpha = coeffs
    @unpack p, U, U_, L, L_, F, Dmin, tf, lambdas = cache

    f = cache.odef
    initial_extrapolation = cache.initial_extrap
    step_number = cache.step_number[]
    maxiters = (step_number == 1 ? 10 + cache.maxiters : cache.maxiters)

    s = length(b)
    tj = ttj[1]
    te = ttj[2]
    indices = eachindex(uj)
    #indices = firstindex(uj):lastindex(uj)
    dtmax = abs((tf - ttj[1]) - ttj[2])

    uiType = eltype(uj)
    realuiType = real(uiType)

    dt = dts[1]
    dtprev = dts[2]
    signdt = dts[3]
    sdt = signdt * dt

    #if (dt != dtprev)
    #    ExtrapolationCoefficients!(nu, mu, c,  sdt, signdt*dtprev, realuiType)
    #end

    step_retcode = true

    if initial_extrapolation && step_number > 1
        @inbounds begin
            for is in 1:s
                for k in indices
                    dUik = muladd(nu[is, 1], L[1][k], ej[k])
                    for js in 2:s
                        dUik = muladd(nu[is, js], L[js][k], dUik)
                    end
                    U[is][k] = uj[k] + dUik
                end
            end
        end
    else
        @inbounds begin
            for is in 1:s
                for k in indices
                    U[is][k] = uj[k] + ej[k]
                end
            end
        end
    end

    Dmin .= Inf

    iter = true
    plusIt = true
    diffU = false
    j_iter = 0
    nf = 0

    while (j_iter < maxiters && iter)
        iter = false
        j_iter += 1

        nf += s
        @inbounds begin
            for is in 1:s
                f(F[is], U[is], p, muladd(sdt, c[is], tj))
                for k in indices
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
            end
        end

        for is in 1:s
            for k in indices
                U_[is][k] = U[is][k]
                dUik = muladd(mu[is, 1], L[1][k], ej[k])
                for js in 2:s
                    dUik = muladd(mu[is, js], L[js][k], dUik)
                end
                U[is][k] = uj[k] + dUik
            end
        end

        diffU = false

        for k in indices
            DY = abs(U[1][k] - U_[1][k])
            for is in 2:s
                DY = max(abs(U[is][k] - U_[is][k]), DY)
            end

            if DY > 0
                diffU = true
                if DY < Dmin[k]
                    Dmin[k] = DY
                    iter = true
                end
            end
        end

        if (!iter && diffU && plusIt)
            iter = true
            plusIt = false
        else
            plusIt = true
        end
    end  # while

    if iter  # iter=true implies that j_iter==maxiters
        @warn "Interrupted. Reached maximum number of iterations (maxiters=$maxiters). The value dt=$dt may be too large."
        step_retcode = false
    end

    if step_retcode
        @inbounds if diffU
            j_iter += 1

            nf += s
            @inbounds begin
                for is in 1:s
                    f(F[is], U[is], p, muladd(sdt, c[is], tj))
                    for k in indices
                        L[is][k] = sdt * (b[is] * F[is][k])
                    end
                end
            end
        end

        #Equivalent to compensated summation

        @inbounds begin
            for k in indices
                L_sum = L[1][k]
                for is in 2:s
                    L_sum += L[is][k]
                end
                res = Base.TwicePrecision(uj[k], ej[k]) + L_sum

                uj[k] = res.hi
                ej[k] = res.lo
            end
        end

        if abs(sdt) >= dtmax
            ttj[1] = tf
            ttj[2] = 0
        else
            res = Base.TwicePrecision(tj, te) + sdt
            ttj[1] = res.hi
            ttj[2] = res.lo
        end

        dtmax = abs((tf - ttj[1]) - ttj[2])
        dts[1] = min(dt, dtmax)
        dts[2] = dt

        stats.nfpiter += j_iter
        stats.nf += nf
    end

    return step_retcode
end

function IRKNGLstep_fixed_2nd!(
        ttj::Array{tType, 1},
        uj::uType,
        ej::uType,
        dts::Array{tType, 1},
        stats::SciMLBase.DEStats,
        coeffs::tcoeffs{tType},
        cache::tcache{uType, realuType, tType, fT, pT}
    ) where {
        uType, realuType, tType, fT, pT,
    }
    @unpack mu, c, b, nu, alpha, eta = coeffs
    @unpack p, U, U_, L, L_, F, Dmin, tf, lambdas = cache

    f = cache.odef
    initial_extrap = cache.initial_extrap
    step_number = cache.step_number[]
    maxiters = (step_number == 1 ? 10 + cache.maxiters : cache.maxiters)

    len = length(uj)
    #lenq = cache.length_q
    lenq = div(len, 2)
    tf = cache.tf

    s = length(b)
    tj = ttj[1]
    te = ttj[2]
    dtmax = abs((tf - ttj[1]) - ttj[2])

    indices = eachindex(uj)
    indices1 = 1:lenq
    indices2 = (lenq + 1):len

    dt = dts[1]
    dtprev = dts[2]
    signdt = dts[3]
    sdt = dt * signdt

    j_iter = 0  # counter of fixed_point iterations
    nf = 0
    nf2 = 0

    step_retcode = true

    if initial_extrap && step_number > 1
        for is in 1:s
            for k in indices2
                dUik = muladd(nu[is, 1], L[1][k], ej[k])
                for js in 2:s
                    dUik = muladd(nu[is, js], L[js][k], dUik)
                end
                U[is][k] = uj[k] + dUik
            end
        end
    else
        for is in 1:s
            for k in indices2
                U[is][k] = uj[k] + ej[k]
            end
        end
    end

    for is in 1:s
        #nf += 1
        for k in indices1
            L[is][k] = sdt * (b[is] * U[is][k + lenq])
        end
    end

    for is in 1:s
        for k in indices1
            dUik = muladd(mu[is, 1], L[1][k], ej[k])
            for js in 2:s
                dUik = muladd(mu[is, js], L[js][k], dUik)
            end
            U[is][k] = uj[k] + dUik
        end
    end

    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt = true
    diffU = false

    @inbounds while (j_iter < maxiters && iter)
        iter = false
        j_iter += 1

        nf2 += s
        for is in 1:s
            f(F[is], U[is], p, tj + sdt * c[is])
            for k in indices2
                L[is][k] = sdt * (b[is] * F[is][k])
            end
        end

        for is in 1:s
            for k in indices1
                kv = k + lenq
                dUik = eta[is, 1] * L[1][kv]
                for js in 2:s
                    dUik = muladd(eta[is, js], L[js][kv], dUik)
                end
                U_[is][k] = U[is][k]
                U[is][k] = uj[k] + (sdt * (c[is] * uj[kv]) + ej[k]) + sdt * dUik
            end
        end

        diffU = false

        for k in indices1
            DY = abs(U[1][k] - U_[1][k])
            for is in 2:s
                DY = max(abs(U[is][k] - U_[is][k]), DY)
            end

            if DY > 0
                diffU = true
                if DY < Dmin[k]
                    Dmin[k] = DY
                    iter = true
                end
            end
        end

        if (!iter && diffU && plusIt)
            iter = true
            plusIt = false
        else
            plusIt = true
        end
    end # while

    if iter  # iter=true implies that j_iter==maxiters
        @warn "Interrupted. Reached maximum number of iterations (maxiters=$maxiters). The value dt=$dt may be too large."
        step_retcode = false
    end

    if step_retcode
        @inbounds if diffU
            j_iter += 1

            for is in 1:s
                nf2 += 1
                f(F[is], U[is], p, tj + sdt * c[is])
                for k in indices2
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
            end

            for is in 1:s
                for k in indices1
                    kv = k + lenq
                    dUik = eta[is, 1] * L[1][kv]
                    for js in 2:s
                        dUik = muladd(eta[is, js], L[js][kv], dUik)
                    end
                    U_[is][k] = U[is][k]
                    U[is][k] = uj[k] + (sdt * (c[is] * uj[kv]) + ej[k]) + sdt * dUik
                end
            end
        end

        @inbounds for k in indices1    #Equivalent to compensated summation
            kv = k + lenq
            Sum = sdt * (1 - c[1]) * L[1][kv]
            for is in 2:s
                Sum = muladd(sdt * (1 - c[is]), L[is][kv], Sum)
            end
            res = Base.TwicePrecision(uj[k], ej[k]) +
                sdt * Base.TwicePrecision(uj[kv], ej[kv]) + Sum
            uj[k] = res.hi
            ej[k] = res.lo
        end

        @inbounds for k in indices2    #Equivalent to compensated summation
            L_sum = L[1][k]
            for is in 2:s
                L_sum += L[is][k]
            end
            res = Base.TwicePrecision(uj[k], ej[k]) + L_sum
            uj[k] = res.hi
            ej[k] = res.lo
        end

        if abs(sdt) >= dtmax
            ttj[1] = tf
            ttj[2] = 0
        else
            res = Base.TwicePrecision(tj, te) + sdt
            ttj[1] = res.hi
            ttj[2] = res.lo
        end

        dtmax = abs((tf - ttj[1]) - ttj[2])
        dts[1] = min(abs(sdt), dtmax)
        dts[2] = dt

        stats.nfpiter += j_iter
        stats.nf += nf
        stats.nf2 += nf2
    end

    return step_retcode
end
