
#
#  IRKstep_SIMD_adaptive!
#  IRKNGLstep_SIMD_adaptive_simpl!

function IRKstep_SIMD_adaptive!(ttj::Array{tType, 1},
        uj::uType,
        ej::uType,
        dts::Array{tType, 1},
        stats::SciMLBase.DEStats,
        coeffs::tcoeffs_SIMD{floatT},
        cache::IRKGL_SIMD_Cache{realuType, floatT, fType, pType, s_, dim_}) where {
        uType, tType, realuType, floatT, fType, pType, s_, dim_}
    @unpack mu, c, b, nu, alpha, X, Y, Z, nu_ = coeffs
    @unpack p, abstol, reltol, U, U_, L, L_, F, Dmin, tf, lambdas = cache

    f = cache.odef
    initial_extrap = cache.initial_extrap
    step_number = cache.step_number[]
    maxiters = (step_number == 1 ? 10 + cache.maxiters : cache.maxiters)
    maxtrials = (step_number == 1 ? 4 * cache.maxtrials : cache.maxtrials)

    uiType = eltype(uj)
    realuiType = real(uiType)

    lambda = lambdas[1]
    lambdaprev = lambdas[2]

    dt = dts[1]
    dtprev = dts[2]
    signdt = dts[3]
    sdt = signdt * dt

    s = length(b)
    pow = realuiType(1 / (2 * s))

    tj = ttj[1]
    te = ttj[2]
    len = length(uj)
    indices = eachindex(uj)
    dtmax = abs((tf - ttj[1]) - ttj[2])

    accept = false
    estimate = zero(eltype(uj))

    j_iter = 0
    ntrials = 0

    L_.data .= L.data

    step_retcode = true

    nf = 0
    nreject = 0
    diffU = false

    while (!accept && ntrials < maxtrials)
        if initial_extrap && step_number > 1
            if (dt != dtprev)
                #ExtrapolationCoefficients!(nu, mu, c,  sdt, signdt*dtprev, realuiType)
                gamma = sdt / (signdt * dtprev)
                #@. Z = gamma*c
                for is in 1:s
                    Z[is] = gamma * sdt * c[is]
                end
                nu_ .= -(PolInterp(sdt * X, Y, Z))'
                nu = VecArray{s, realuiType, 2}(nu_)
            end

            for k in indices
                Lk = getindex_(L_, k)
                dUk = muladd(nu[1], Lk[1], ej[k])
                for js in 2:s
                    dUk = muladd(nu[js], Lk[js], dUk)
                end
                setindex_!(U, uj[k] + dUk, k)
            end

        else
            for k in indices   # 2024-07-29 behin-behinekoa expolaziorik gabe
                uej = uj[k] + ej[k]
                setindex_!(U, uej, k)
            end
        end

        iter = true
        plusIt = true
        diffU = false
        #j_iter = 0

        Dmin .= Inf

        while (j_iter < maxiters && iter)
            iter = false
            j_iter += 1

            U_.data .= U.data
            nf += s
            f(F, U, p, tj + sdt * c)

            for k in indices
                Fk = getindex_(F, k)
                Lk = sdt * (b * Fk)
                dUk = muladd(mu[1], Lk[1], ej[k])
                for is in 2:s
                    dUk = muladd(mu[is], Lk[is], dUk)
                end
                Uk = uj[k] + dUk
                setindex_!(U, Uk, k)
                setindex_!(L, Lk, k)
                Uk_ = getindex_(U_, k)
                DY = maximum(abs(Uk - Uk_))

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
        end # while iter

        ntrials += 1

        #estimate = ErrorEst(U, F, dt, alpha, abstol, reltol)
        estimate = ErrorEst_SIMD(U, F, len, indices, dt, alpha, abstol, reltol)

        lambda = (estimate)^pow

        if (estimate < 2 && j_iter < maxiters)
            accept = true
        else
            nreject += 1
            #dt = dt / lambda
            dt = min(dt / lambda, dtmax)
            sdt = signdt * dt
        end
    end # while accept

    if (!accept && ntrials == maxtrials)
        @warn("Failure (adaptive step): maximum number of trials=", maxtrials,
            " at step=", step_number,
            " dt=", dts[1])

        step_retcode = false
    end

    if step_retcode
        @inbounds if diffU
            j_iter += 1
            nf += s
            f(F, U, p, tj + sdt * c)

            for k in indices
                Fk = getindex_(F, k)
                Lk = sdt * (b * Fk)
                setindex_!(L, Lk, k)
            end
        end

        #Equivalent to compensated summation

        @inbounds for k in indices    #Equivalent to compensated summation
            Lk = getindex_(L, k)
            L_sum = sum(Lk)
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

        if (step_number == 1)
            dts[1] = min(max(dt / 2, min(2 * dt, dt / lambda)), dtmax)
        else
            hath1 = dt / lambda
            hath2 = dtprev / lambdaprev
            tildeh = hath1 * (hath1 / hath2)^(lambda / lambdaprev)
            barlamb1 = (dt + tildeh) / (hath1 + tildeh)
            barlamb2 = (dtprev + dt) / (hath2 + hath1)
            barh = hath1 * (hath1 / hath2)^(barlamb1 / barlamb2)
            dts[1] = min(max(dt / 2, min(2 * dt, barh)), dtmax)
        end

        dts[2] = dt
        lambdas[2] = lambda

        stats.nfpiter += j_iter
        stats.nf += nf
        stats.nreject += nreject
    end

    return step_retcode
end

function IRKNGLstep_SIMD_adaptive_simpl!(ttj::Array{tType, 1},
        uj::uType,
        ej::uType,
        dts::Array{tType, 1},
        stats::SciMLBase.DEStats,
        coeffs::tcoeffs_SIMD{floatT},
        cache::IRKGL_SIMD_Cache{realuType, floatT, fType, pType, s_, dim_}) where {
        uType, tType, realuType, floatT, fType, pType, s_, dim_}
    @unpack mu, c, b, nu, alpha, X, Y, Z, nu_ = coeffs
    @unpack p, abstol, reltol, U, U_, L, L_, F, Dmin, tf, lambdas = cache

    f = cache.odef
    initial_extrap = cache.initial_extrap
    step_number = cache.step_number[]
    maxiters = (step_number == 1 ? 10 + cache.maxiters : cache.maxiters)
    maxtrials = (step_number == 1 ? 4 * cache.maxtrials : cache.maxtrials)

    len = length(uj)
    #lenq = cache.length_q
    lenq = div(len, 2)

    uiType = eltype(uj)
    realuiType = real(uiType)

    lambda = lambdas[1]
    lambdaprev = lambdas[2]

    dt = dts[1]
    dtprev = dts[2]
    signdt = dts[3]
    sdt = signdt * dt

    s = length(b)
    pow = realuiType(1 / (2 * s))

    tj = ttj[1]
    te = ttj[2]
    dtmax = abs((tf - ttj[1]) - ttj[2])

    len = length(uj)
    indices = eachindex(uj)
    indices1 = 1:lenq
    indices2 = (lenq + 1):len

    accept = false
    estimate = zero(eltype(uj))

    j_iter = 0
    ntrials = 0

    L_.data .= L.data

    step_retcode = true
    nf = 0
    nf2 = 0
    nreject = 0
    diffU = false

    while (!accept && ntrials < maxtrials)
        if (dt != dtprev)
            #ExtrapolationCoefficients!(nu, mu, c,  sdt, signdt*dtprev, realuiType)
            if dtprev == 0
                nu_ .= zeros(realuiType, s, s)
            else
                gamma = sdt / (signdt * dtprev)
                #@. Z = gamma*c
                for is in 1:s
                    Z[is] = gamma * sdt * c[is]
                end
                nu_ .= -(PolInterp(sdt * X, Y, Z))'
            end

            nu = VecArray{s, realuiType, 2}(nu_)
        end

        for k in indices2
            Lk = getindex_(L_, k)
            dUk = muladd(nu[1], Lk[1], ej[k])
            for is in 2:s
                dUk = muladd(nu[is], Lk[is], dUk)
            end
            setindex_!(U, uj[k] + dUk, k)
        end

        nf += s
        for k in indices1
            Uk = getindex_(U, k + lenq)
            setindex_!(F, Uk, k)
            Fk = getindex_(F, k)
            Lk = sdt * (b * Fk)
            setindex_!(L, Lk, k)
            dUk = muladd(mu[1], Lk[1], ej[k])
            for is in 2:s
                dUk = muladd(mu[is], Lk[is], dUk)
            end
            setindex_!(U, uj[k] + dUk, k)
        end

        Dmin .= Inf

        iter = true
        plusIt = true
        diffU = false

        while (j_iter < maxiters && iter)
            iter = false
            j_iter += 1

            U_.data .= U.data

            nf2 += s
            f(F, U, p, tj + sdt * c)

            for k in indices2
                Fk = getindex_(F, k)
                Lk = sdt * (b * Fk)
                setindex_!(L, Lk, k)
                dUk = muladd(mu[1], Lk[1], ej[k])
                for is in 2:s
                    dUk = muladd(mu[is], Lk[is], dUk)
                end
                setindex_!(U, uj[k] + dUk, k)
            end

            nf += s
            for k in indices1
                Uk = getindex_(U, k + lenq)
                setindex_!(F, Uk, k)
                Fk = getindex_(F, k)
                Lk = sdt * (b * Fk)
                setindex_!(L, Lk, k)
                dUk = muladd(mu[1], Lk[1], ej[k])
                for is in 2:s
                    dUk = muladd(mu[is], Lk[is], dUk)
                end
                setindex_!(U, uj[k] + dUk, k)
            end

            diffU = false

            for k in indices
                Uk = getindex_(U, k)
                Uk_ = getindex_(U_, k)
                DY = maximum(abs(Uk - Uk_))

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
        end # while iter

        ntrials += 1

        #estimate = ErrorEst(U, F, dt, alpha, abstol, reltol)
        estimate = ErrorEst_SIMD(U, F, len, indices, dt, alpha, abstol, reltol)

        lambda = (estimate)^pow

        if (estimate < 2 && j_iter < maxiters)
            accept = true
        else
            nreject += 1
            #dt = dt / lambda
            dt = min(dt / lambda, dtmax)
            sdt = signdt * dt
        end
    end # while accept

    if (!accept && ntrials == maxtrials)
        @warn("Failure (adaptive step): maximum number of trials=", maxtrials,
            " at step=", step_number,
            " dt=", dts[1])
        step_retcode = false
    end

    if step_retcode
        @inbounds if diffU
            j_iter += 1

            nf2 += s
            f(F, U, p, tj + sdt * c)

            for k in indices2
                Fk = getindex_(F, k)
                Lk = sdt * (b * Fk)
                dUk = muladd(mu[1], Lk[1], ej[k])
                for is in 2:s
                    dUk = muladd(mu[is], Lk[is], dUk)
                end
                setindex_!(U, uj[k] + dUk, k)
                setindex_!(L, Lk, k)
            end

            nf += s
            for k in indices1
                Uk = getindex_(U, k + lenq)
                setindex_!(F, Uk, k)
                Fk = getindex_(F, k)
                Lk = sdt * (b * Fk)
                setindex_!(L, Lk, k)
            end
        end

        #Equivalent to compensated summation

        @inbounds for k in indices
            Lk = getindex_(L, k)
            L_sum = sum(Lk)
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

        if (step_number == 1)
            dts[1] = min(max(dt / 2, min(2 * dt, dt / lambda)), dtmax)
        else
            hath1 = dt / lambda
            hath2 = dtprev / lambdaprev
            tildeh = hath1 * (hath1 / hath2)^(lambda / lambdaprev)
            barlamb1 = (dt + tildeh) / (hath1 + tildeh)
            barlamb2 = (dtprev + dt) / (hath2 + hath1)
            barh = hath1 * (hath1 / hath2)^(barlamb1 / barlamb2)
            dts[1] = min(max(dt / 2, min(2 * dt, barh)), dtmax)
        end
        dts[2] = dt
        lambdas[2] = lambda

        stats.nfpiter += j_iter
        stats.nf += nf
        stats.nf2 += nf2
        stats.nreject += nreject
    end

    return step_retcode
end
