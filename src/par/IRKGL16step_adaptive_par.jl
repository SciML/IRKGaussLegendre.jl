#
#  IRKstep_par_adaptive!
#  IRKstepDynODE_par_adaptive!
#  IRKNGLstep_adaptive_par_simpl!

function IRKstep_par_adaptive!(
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
    @unpack mu, c, b, nu, alpha, X, Y, Z = coeffs
    @unpack p, abstol, reltol, U, U_, L, L_, F, Dmin, tf, lambdas = cache

    f = cache.odef
    initial_extrapolation = cache.initial_extrap
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
    indices = eachindex(uj)
    dtmax = abs((tf - ttj[1]) - ttj[2])

    accept = false
    estimate = zero(eltype(uj))

    j_iter = 0
    ntrials = 0

    for is in 1:s
        L_[is] .= L[is]
    end

    step_retcode = true

    nf = 0
    nreject = 0
    diffU = false

    while (!accept && ntrials < maxtrials)
        if initial_extrapolation && step_number > 1
            if (dt != dtprev)
                #ExtrapolationCoefficients!(nu, mu, c,  sdt, signdt*dtprev, realuiType)
                gamma = sdt / (signdt * dtprev)
                for is in 1:s
                    Z[is] = gamma * sdt * c[is]
                end
                nu .= -(PolInterp(sdt * X, Y, Z))'
            end

            @inbounds begin
                for is in 1:s
                    for k in indices
                        dUik = muladd(nu[is, 1], L_[1][k], ej[k])
                        for js in 2:s
                            dUik = muladd(nu[is, js], L_[js][k], dUik)
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

        iter = true
        plusIt = true
        diffU = false
        #j_iter = 0

        Dmin .= Inf

        while (j_iter < maxiters && iter)
            iter = false
            j_iter += 1

            nf += s
            @inbounds begin
                Threads.@threads for is in 1:s
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
        end # while iter

        ntrials += 1

        estimate = ErrorEst(U, F, dt, alpha, abstol, reltol)
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
        @warn(
            "Failure (adaptive step): maximum number of trials=", maxtrials,
            " at step=", step_number,
            " dt=", dts[1]
        )

        step_retcode = false
    end

    if step_retcode
        @inbounds if diffU
            j_iter += 1

            nf += s
            @inbounds begin
                Threads.@threads for is in 1:s
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

function IRKstepDynODE_par_adaptive!(
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
    @unpack mu, c, b, nu, alpha, X, Y, Z = coeffs
    @unpack p, abstol, reltol, U, U_, L, L_, F, Dmin, tf, lambdas = cache

    f = cache.odef
    initial_extrap = cache.initial_extrap
    step_number = cache.step_number[]
    maxiters = (step_number == 1 ? 10 + cache.maxiters : cache.maxiters)
    maxtrials = (step_number == 1 ? 4 * cache.maxtrials : cache.maxtrials)
    f1 = f.f1
    f2 = f.f2

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

    indices = eachindex(uj)
    indices1 = eachindex(uj.x[1])
    indices2 = indices1[end] .+ eachindex(uj.x[2])

    accept = false
    estimate = zero(eltype(uj))

    j_iter = 0
    ntrials = 0

    for is in 1:s
        L_[is] .= L[is]
    end

    step_retcode = true
    nf = 0
    nf2 = 0
    nreject = 0
    diffU = false

    while (!accept && ntrials < maxtrials)
        if (dt != dtprev)
            # ExtrapolationCoefficients!(nu, mu, c,  sdt, signdt*dtprev, realuiType)
            if dtprev == 0
                nu .= zeros(realuiType, s, s)
            else
                gamma = sdt / (signdt * dtprev)
                #@. Z = gamma*c
                for is in 1:s
                    Z[is] = gamma * sdt * c[is]
                end
                nu .= -(PolInterp(sdt * X, Y, Z))'
            end
        end

        for is in 1:s
            for k in indices2
                dUik = muladd(nu[is, 1], L_[1][k], ej[k])
                for js in 2:s
                    dUik = muladd(nu[is, js], L_[js][k], dUik)
                end
                U[is][k] = uj[k] + dUik
            end
        end

        nf += s
        @inbounds begin
            for is in 1:s
                f1(F[is].x[1], U[is].x[1], U[is].x[2], p, muladd(sdt, c[is], tj))
                for k in indices1
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
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

        iter = true
        plusIt = true
        diffU = false

        Dmin .= Inf

        while (j_iter < maxiters && iter)
            iter = false
            j_iter += 1

            nf2 += s
            Threads.@threads for is in 1:s
                f2(F[is].x[2], U[is].x[1], U[is].x[2], p, muladd(sdt, c[is], tj))
                for k in indices2
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
            end

            for is in 1:s
                for k in indices2
                    dUik = muladd(mu[is, 1], L[1][k], ej[k])
                    for js in 2:s
                        dUik = muladd(mu[is, js], L[js][k], dUik)
                    end
                    U_[is][k] = U[is][k]
                    U[is][k] = uj[k] + dUik
                end
            end

            nf += s
            for is in 1:s
                f1(F[is].x[1], U[is].x[1], U[is].x[2], p, muladd(sdt, c[is], tj))
                for k in indices1
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
            end

            for is in 1:s
                for k in indices1
                    dUik = muladd(mu[is, 1], L[1][k], ej[k])
                    for js in 2:s
                        dUik = muladd(mu[is, js], L[js][k], dUik)
                    end
                    U_[is][k] = U[is][k]
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
        end # while iter

        ntrials += 1

        estimate = ErrorEst(U, F, dt, alpha, abstol, reltol)
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
        @warn(
            "Failure (adaptive step): maximum number of trials=", maxtrials,
            " at step=", step_number,
            " dt=", dts[1]
        )
        step_retcode = false
    end

    if step_retcode
        @inbounds if diffU
            j_iter += 1

            nf2 += s
            Threads.@threads for is in 1:s
                f2(F[is].x[2], U[is].x[1], U[is].x[2], p, muladd(sdt, c[is], tj))
                for k in indices2
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
            end

            for is in 1:s
                for k in indices2
                    dUik = muladd(mu[is, 1], L[1][k], ej[k])
                    for js in 2:s
                        dUik = muladd(mu[is, js], L[js][k], dUik)
                    end
                    U[is][k] = uj[k] + dUik
                end
            end

            nf += s
            for is in 1:s
                f1(F[is].x[1], U[is].x[1], U[is].x[2], p, muladd(sdt, c[is], tj))
                for k in indices1
                    L[is][k] = sdt * (b[is] * F[is][k])
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

function IRKNGLstep_par_adaptive_simpl!(
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
    @unpack mu, c, b, nu, alpha, X, Y, Z = coeffs
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

    indices = eachindex(uj)
    indices1 = 1:lenq
    indices2 = (lenq + 1):len

    accept = false
    estimate = zero(eltype(uj))

    j_iter = 0
    ntrials = 0

    for is in 1:s
        L_[is] .= L[is]
    end

    step_retcode = true
    nf = 0
    nf2 = 0
    nreject = 0
    diffU = false

    while (!accept && ntrials < maxtrials)
        if (dt != dtprev)
            # ExtrapolationCoefficients!(nu, mu, c,  sdt, signdt*dtprev, realuiType)
            if dtprev == 0
                nu .= zeros(realuiType, s, s)
            else
                gamma = sdt / (signdt * dtprev)
                #@. Z = gamma*c
                for is in 1:s
                    Z[is] = gamma * sdt * c[is]
                end
                nu .= -(PolInterp(sdt * X, Y, Z))'
            end
        end

        for is in 1:s
            for k in indices2
                dUik = muladd(nu[is, 1], L_[1][k], ej[k])
                for js in 2:s
                    dUik = muladd(nu[is, js], L_[js][k], dUik)
                end
                U[is][k] = uj[k] + dUik
            end
        end

        nf += s
        @inbounds begin
            for is in 1:s
                for k in indices1
                    F[is][k] = U[is][k + lenq]
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
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

        iter = true
        plusIt = true
        diffU = false

        Dmin .= Inf

        while (j_iter < maxiters && iter)
            iter = false
            j_iter += 1

            nf2 += s
            Threads.@threads for is in 1:s
                f(F[is], U[is], p, muladd(sdt, c[is], tj))
                for k in indices2
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
            end

            for is in 1:s
                for k in indices2
                    dUik = muladd(mu[is, 1], L[1][k], ej[k])
                    for js in 2:s
                        dUik = muladd(mu[is, js], L[js][k], dUik)
                    end
                    U_[is][k] = U[is][k]
                    U[is][k] = uj[k] + dUik
                end
            end

            nf += s
            for is in 1:s
                for k in indices1
                    F[is][k] = U[is][k + lenq]
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
            end

            for is in 1:s
                for k in indices1
                    dUik = muladd(mu[is, 1], L[1][k], ej[k])
                    for js in 2:s
                        dUik = muladd(mu[is, js], L[js][k], dUik)
                    end
                    U_[is][k] = U[is][k]
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
        end # while iter

        ntrials += 1

        estimate = ErrorEst(U, F, dt, alpha, abstol, reltol)
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
        @warn(
            "Failure (adaptive step): maximum number of trials=", maxtrials,
            " at step=", step_number,
            " dt=", dts[1]
        )
        step_retcode = false
    end

    if step_retcode
        @inbounds if diffU
            j_iter += 1

            nf2 += s
            Threads.@threads for is in 1:s
                f(F[is], U[is], p, muladd(sdt, c[is], tj))
                for k in indices2
                    L[is][k] = sdt * (b[is] * F[is][k])
                end
            end

            for is in 1:s
                for k in indices2
                    dUik = muladd(mu[is, 1], L[1][k], ej[k])
                    for js in 2:s
                        dUik = muladd(mu[is, js], L[js][k], dUik)
                    end
                    U[is][k] = uj[k] + dUik
                end
            end

            nf += s
            for is in 1:s
                for k in indices1
                    F[is][k] = U[is][k + lenq]
                    L[is][k] = sdt * (b[is] * F[is][k])
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
