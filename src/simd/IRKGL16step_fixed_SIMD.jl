#
#  IRKGL16Step fixed SIMD functions

#      IRKGLstep_SIMD_fixed!
#      IRKNGLstep_SIMD_fixed_simpl!

function IRKGLstep_SIMD_fixed!(
        ttj::Array{tType, 1},
        uj::uType,
        ej::uType,
        dts::Array{tType, 1},
        stats::SciMLBase.DEStats,
        coeffs::tcoeffs_SIMD{floatT},
        cache::IRKGL_SIMD_Cache{realuType, floatT, fType, pType, s_, dim_}
    ) where {
        uType, tType, realuType, floatT, fType, pType, s_, dim_,
    }
    @unpack mu, c, b, nu = coeffs
    @unpack p, U, U_, L, F, Dmin, tf = cache

    f = cache.odef
    initial_extrap = cache.initial_extrap
    step_number = cache.step_number[]
    maxiters = (step_number == 1 ? 10 + cache.maxiters : cache.maxiters)

    s = length(b)
    tj = ttj[1]
    te = ttj[2]
    indices = eachindex(uj)
    dtmax = abs((tf - ttj[1]) - ttj[2])

    dt = dts[1]
    dtprev = dts[2]
    signdt = dts[3]
    sdt = dt * signdt

    step_retcode = true

    if initial_extrap && step_number > 1
        for k in indices
            Lk = L[k]
            dUk = muladd(nu[1], Lk[1], ej[k])
            for js in 2:s
                dUk = muladd(nu[js], Lk[js], dUk)
            end
            U[k] = uj[k] + dUk
        end

    else
        for k in indices
            uej = uj[k] + ej[k]
            U[k] = uej
        end
    end

    #println("U=", U.data[1,1:end,1:end, 1:end])
    #println("norm(U)=", norm(U.data[1,1:end,1:end, 1:end]))

    Dmin .= Inf

    iter = true
    plusIt = true
    diffU = false
    j_iter = 0
    nf = 0

    @inbounds while (j_iter < maxiters && iter)
        iter = false
        j_iter += 1

        U_.data .= U.data
        nf += s
        f(F, U, p, tj + sdt * c)

        diffU = false

        #println("j_iter=",j_iter,",norm(F)=", norm(F.data[1,1:end,1:end, 1:end]))

        for k in indices
            Fk = F[k]
            Lk = sdt * (b * Fk)
            dUk = muladd(mu[1], Lk[1], ej[k])
            for is in 2:s
                dUk = muladd(mu[is], Lk[is], dUk)
            end
            Uk = uj[k] + dUk
            U[k] = Uk
            L[k] = Lk
            Uk_ = U_[k]
            DY = maximum(abs(Uk - Uk_))

            if DY > 0
                diffU = true
                if DY < Dmin[k]
                    Dmin[k] = DY
                    iter = true
                end
            end
        end

        #println("j_iter=",j_iter,",norm(L)=", norm(L.data[1,1:end,1:end, 1:end]))
        #println("j_iter=",j_iter,",norm(U)=", norm(U.data[1,1:end,1:end, 1:end]))
        #println("")

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
            nf += s
            f(F, U, p, tj + sdt * c)

            for k in indices
                Fk = F[k]
                Lk = sdt * (b * Fk)
                L[k] = Lk
            end
        end

        @inbounds for k in indices    #Equivalent to compensated summation
            Lk = L[k]
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

        dts[1] = min(abs(sdt), dtmax)
        dts[2] = dt

        stats.nfpiter += j_iter
        stats.nf += nf
    end

    #println("j=", step_number, ",nit=", j_iter)
    return step_retcode
end

function IRKNGLstep_SIMD_fixed_simpl!(
        ttj::Array{tType, 1},
        uj::uType,
        ej::uType,
        dts::Array{tType, 1},
        stats::SciMLBase.DEStats,
        coeffs::tcoeffs_SIMD{floatT},
        cache::IRKGL_SIMD_Cache{realuType, floatT, fType, pType, s_, dim_}
    ) where {
        uType, tType, realuType, floatT, fType, pType, s_, dim_,
    }
    @unpack mu, c, b, nu = coeffs
    @unpack p, U, U_, L, F, Dmin, tf = cache

    f = cache.odef
    initial_extrap = cache.initial_extrap
    step_number = cache.step_number[]
    maxiters = (step_number == 1 ? 10 + cache.maxiters : cache.maxiters)

    s = length(b)
    tj = ttj[1]
    te = ttj[2]

    len = length(uj)
    #lenq = cache.length_q
    lenq = div(len, 2)
    indices = eachindex(uj)
    indices1 = 1:lenq
    indices2 = (lenq + 1):len
    dtmax = abs((tf - ttj[1]) - ttj[2])

    dt = dts[1]
    dtprev = dts[2]
    signdt = dts[3]
    sdt = dt * signdt

    j_iter = 0  # counter of fixed_point iterations
    nf = 0
    nf2 = 0

    step_retcode = true

    if initial_extrap && step_number > 1
        for k in indices2
            Lk = L[k]
            dUk = muladd(nu[1], Lk[1], ej[k])
            for is in 2:s
                dUk = muladd(nu[is], Lk[is], dUk)
            end
            U[k] = uj[k] + dUk
        end

    else
        for k in indices2
            uej = uj[k] + ej[k]
            U[k] = uej
        end
    end

    nf += s
    for k in indices1
        Uk = U[k + lenq]
        F[k] = Uk
        Fk = F[k]
        Lk = sdt * (b * Fk)
        L[k] = Lk
        dUk = muladd(mu[1], Lk[1], ej[k])
        for is in 2:s
            dUk = muladd(mu[is], Lk[is], dUk)
        end
        U[k] = uj[k] + dUk
    end

    Dmin .= Inf

    iter = true # Initialize iter outside the for loop
    plusIt = true
    diffU = true

    @inbounds while (j_iter < maxiters && iter)
        iter = false
        j_iter += 1

        U_.data .= U.data

        nf2 += s
        f(F, U, p, tj + sdt * c)

        for k in indices2
            Fk = F[k]
            Lk = sdt * (b * Fk)
            L[k] = Lk
            dUk = muladd(mu[1], Lk[1], ej[k])
            for is in 2:s
                dUk = muladd(mu[is], Lk[is], dUk)
            end
            U[k] = uj[k] + dUk
        end

        nf += s
        for k in indices1
            Uk = U[k + lenq]
            F[k] = Uk
            Fk = F[k]
            Lk = sdt * (b * Fk)
            L[k] = Lk
            dUk = muladd(mu[1], Lk[1], ej[k])
            for is in 2:s
                dUk = muladd(mu[is], Lk[is], dUk)
            end
            U[k] = uj[k] + dUk
        end

        diffU = false

        for k in indices1
            Uk = U[k]
            Uk_ = U_[k]
            DY = maximum(abs(Uk - Uk_))

            if DY > 0
                diffU = true
                if DY < Dmin[k]
                    Dmin[k] = DY
                    iter = true
                end
            end
        end

        if (!iter && diffU && plusIt)  #
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
        @inbounds if (j_iter < maxiters && diffU)
            j_iter += 1

            nf2 += s
            f(F, U, p, tj + sdt * c)

            for k in indices2
                Fk = F[k]
                Lk = sdt * (b * Fk)
                dUk = muladd(mu[1], Lk[1], ej[k])
                for is in 2:s
                    dUk = muladd(mu[is], Lk[is], dUk)
                end
                U[k] = uj[k] + dUk
                L[k] = Lk
            end

            nf += s
            for k in indices1
                Uk = U[k + lenq]
                F[k] = Uk
                Fk = F[k]
                Lk = sdt * (b * Fk)
                L[k] = Lk
            end
        end

        @inbounds for k in indices    #Equivalent to compensated summation
            Lk = L[k]
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

        dts[1] = min(abs(sdt), dtmax)
        dts[2] = dt

        stats.nfpiter += j_iter
        stats.nf += nf
        stats.nf2 += nf2
    end

    return step_retcode
end
