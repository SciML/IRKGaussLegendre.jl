struct tcoeffs{tType}
    mu::Array{tType, 2}
    c::Array{tType, 1}
    b::Array{tType, 1}
    nu::Array{tType, 2}
    alpha::Array{tType, 1}
    X::Array{tType, 1}
    Y::Array{tType, 2}
    Z::Array{tType, 1}
    kappa::Array{tType, 1}
    X2::Array{tType, 1}
    Y2::Array{tType, 2}
    Z2::Array{tType, 1}
    eta::Array{tType, 2}
end

struct tcache{uType, realuType, tType, fT, pT}
    odef::fT # function defining the ODE system
    p::pT # parameters and so
    abstol::tType
    reltol::tType
    U::Array{uType, 1}
    U_::Array{uType, 1}
    L::Array{uType, 1}
    L_::Array{uType, 1}
    F::Array{uType, 1}
    Dmin::realuType
    maxiters::Int64
    maxtrials::Int64
    step_number::Array{Int64, 0}
    initial_extrap::Bool
    length_u::Int64
    length_q::Int64
    tf::tType
    lambdas::Array{tType, 1}
end

struct tcoeffs_SIMD{floatT, s_}
    b::Vec{s_, floatT}
    c::Vec{s_, floatT}
    mu::VecArray{s_, floatT, 2}
    nu::VecArray{s_, floatT, 2}
    alpha::Vec{s_, floatT}
    X::Array{floatT, 1}
    Y::Array{floatT, 2}
    Z::Array{floatT, 1}
    nu_::Array{floatT, 2}
    kappa::Vec{s_, floatT}
    kappa_::Array{floatT, 1}
    X2::Array{floatT, 1}
    Y2::Array{floatT, 2}
    Z2::Array{floatT, 1}
    eta::VecArray{s_, floatT, 2}
end

struct IRKGL_SIMD_Cache{realuType, floatT, fType, pType, s_, dim_}
    odef::fType # function defining the ODE system
    p::pType # parameters and so
    abstol::floatT
    reltol::floatT
    U::VecArray{s_, floatT, dim_}
    U_::VecArray{s_, floatT, dim_}
    L::VecArray{s_, floatT, dim_}
    L_::VecArray{s_, floatT, dim_}
    F::VecArray{s_, floatT, dim_}
    Dmin::realuType
    maxiters::Int64
    maxtrials::Int64
    step_number::Array{Int64, 0}
    initial_extrap::Bool
    length_u::Int64
    length_q::Int64
    tf::floatT
    lambdas::Array{floatT, 1}
end

abstract type IRKAlgorithm{
    s,
    second_order_ode,
    simd,
    fseq,
    maxtrials,
    initial_extrapolation,
} <: SciMLBase.AbstractODEAlgorithm end
struct IRKGL16{
        s,
        second_order_ode,
        simd,
        fseq,
        maxtrials,
        initial_extrapolation,
    } <: IRKAlgorithm{
        s,
        second_order_ode,
        simd,
        fseq,
        maxtrials,
        initial_extrapolation,
    } end
function IRKGL16(;
        s = 8,
        second_order_ode = false,
        #simd = eltype(prob.u0)<:Union{Float32,Float64} ? true : false,
        simd = true,
        fseq = true,
        maxtrials = 5,
        initial_extrapolation = true
    )
    return IRKGL16{
        s,
        second_order_ode,
        simd,
        fseq,
        maxtrials,
        initial_extrapolation,
    }()
end

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem{
            uType,
            tspanType,
            isinplace,
        },
        alg::IRKGL16{
            s,
            second_order_ode,
            simd,
            fseq,
            maxtrials,
            initial_extrapolation,
        },
        args...;
        dt = zero(eltype(tspanType)),
        maxiters = 100,
        save_everystep = true,
        saveat = nothing,
        adaptive = true,
        reltol = eltype(tspanType)(1.0e-6),
        abstol = eltype(tspanType)(1.0e-6),
        kwargs...
    ) where {
        uType,
        tspanType,
        isinplace,
        s,
        second_order_ode,
        simd,
        fseq,
        maxtrials,
        initial_extrapolation,
    }
    checks = true

    stats = SciMLBase.DEStats(0)
    #stats = DiffEqBase.DEStats(0)
    stats.nf = 0
    stats.nf2 = 0
    stats.nfpiter = 0
    stats.naccept = 0
    stats.nreject = 0

    if (prob.f isa ODEFunction)
        @unpack u0, tspan, p = prob
        f = SciMLBase.unwrapped_f(prob.f)
    else
        @warn("Error: incorrect ODEFunction")
        sol = SciMLBase.build_solution(
            prob, alg, [prob.tspan[1]], [u0], stats = stats, retcode = ReturnCode.Failure
        )
        return (sol)
    end

    t0 = prob.tspan[1]
    tf = prob.tspan[2]

    tType = eltype(tspanType)
    realuType = typeof(real(u0))

    step_fun::Function = empty

    if simd && eltype(u0) <: Union{Float32, Float64}
        floatType = eltype(u0)
        if !adaptive
            if !second_order_ode
                if !fseq
                    step_fun = IRKGLstep_SIMD_fixed!
                else
                    step_fun = IRKGLstep_HYBR_fixed!
                end
            else
                if !fseq
                    step_fun = IRKNGLstep_SIMD_fixed_2nd!
                else
                    step_fun = IRKNGLstep_HYBR_fixed_2nd!
                end
            end
        else
            if !second_order_ode
                if !fseq
                    step_fun = IRKstep_SIMD_adaptive!
                else
                    step_fun = IRKstep_HYBR_adaptive!
                end
            else
                if !fseq
                    step_fun = IRKNGLstep_SIMD_adaptive_2nd!
                else
                    step_fun = IRKNGLstep_HYBR_adaptive_2nd!
                end
            end
        end

    else
        if adaptive
            if !second_order_ode
                step_fun = IRKstep_adaptive!
            else
                step_fun = IRKNGLstep_adaptive_2nd!
            end
        else
            if !second_order_ode
                step_fun = IRKstep_fixed!
            else
                step_fun = IRKNGLstep_fixed_2nd!
            end
        end
    end

    abstol = convert(tType, abstol)
    reltol = convert(tType, reltol)

    if (dt == 0)
        d0 = MyNorm(u0, abstol, reltol)
        du0 = similar(u0)
        f(du0, u0, p, t0)
        d1 = MyNorm(du0, abstol, reltol)
        if (d0 < 1.0e-5 || d1 < 1.0e-5)
            dt = convert(tType, 1.0e-6)
        else
            dt = convert(tType, 0.01 * (d0 / d1))
        end
    end

    dt = min(abs(dt), abs(tf - t0))
    signdt = sign(tspan[2] - tspan[1])
    dts = Array{tType}(undef, 1)

    if !adaptive
        dtprev = dt
    else
        dtprev = zero(tType)
    end

    dts = [dt, dtprev, signdt]
    sdt = signdt * dt

    #   Memory preallocation (IRKL_Cache)

    Dmin = similar(real(u0))
    step_number = Array{Int64, 0}(undef)
    step_number[] = 0
    length_u = length(u0)
    length_q = div(length_u, 2)

    if simd && eltype(u0) <: Union{Float32, Float64}
        coeffs = tcoeffs{tType}(
            zeros(s, s), zeros(s), zeros(s), zeros(s, s),
            zeros(s), zeros(s + 1), zeros(s, s + 1), zeros(s),
            zeros(s), zeros(s + 1), zeros(s, s + 1), zeros(1),
            zeros(s, s)
        )
        mu_ = coeffs.mu
        c_ = coeffs.c
        b_ = coeffs.b
        nu_ = coeffs.nu
        alpha_ = coeffs.alpha
        X = coeffs.X
        Y = coeffs.Y
        Z = coeffs.Z
        kappa_ = coeffs.kappa
        X2 = coeffs.X2
        Y2 = coeffs.Y2
        Z2 = coeffs.Z2
        eta_ = coeffs.eta

        EstimateCoeffs!(s, alpha_, tType)

        GaussLegendreCoefficients!(s, mu_, c_, b_, eta_, tType)
        X .= vcat(-c_[end:-1:1], [0])
        Y .= hcat(mu_, zeros(tType, s))
        Z .= c_
        if dtprev == 0
            nu_ .= zeros(tType, s, s)
        else
            nu_ .= -(PolInterp(X, Y, Z))'
        end

        # Inteporlation
        X2 .= vcat([zero(tType)], c_[1:1:end])
        Y2 .= hcat(zeros(tType, s), mu_')

        dims = size(u0)

        c = vload(Vec{s, floatType}, c_, 1)
        b = vload(Vec{s, floatType}, b_, 1)
        nu = VecArray{s, floatType, 2}(nu_)
        mu = VecArray{s, floatType, 2}(mu_)
        alpha = vload(Vec{s, floatType}, alpha_, 1)
        kappa = vload(Vec{s, floatType}, kappa_, 1)
        eta = VecArray{s, floatType, 2}(eta_)

        zz = zeros(Float64, s, dims...)
        U = VecArray{s, Float64, length(dims) + 1}(zz)
        U_ = deepcopy(U)
        L = deepcopy(U)
        L_ = deepcopy(U)
        F = deepcopy(U)

        coeffs = tcoeffs_SIMD(
            b, c, mu, nu, alpha, X, Y, Z, nu_, kappa, kappa_, X2, Y2, Z2, eta
        )

        cache = IRKGL_SIMD_Cache(
            f, p, abstol, reltol,
            U, U_, L, L_, F,
            Dmin, maxiters, maxtrials, step_number,
            initial_extrapolation, length_u, length_q, tf,
            fill(zero(tType), 2)
        )

    else
        coeffs = tcoeffs{tType}(
            zeros(s, s), zeros(s), zeros(s),
            zeros(s, s), zeros(s), zeros(s + 1), zeros(s, s + 1), zeros(s),
            zeros(s), zeros(s + 1), zeros(s, s + 1), zeros(1),
            zeros(s, s)
        )

        @unpack mu, c, b, nu, alpha, X, Y, Z, kappa, X2, Y2, Z2, eta = coeffs
        EstimateCoeffs!(s, alpha, tType)

        GaussLegendreCoefficients!(s, mu, c, b, eta, tType)
        X .= vcat(-c[end:-1:1], [0])
        Y .= hcat(mu, zeros(tType, s))
        Z .= c
        if dtprev == 0
            nu .= zeros(tType, s, s)
        else
            nu .= -(PolInterp(X, Y, Z))'
        end

        # Inteporlation
        X2 .= vcat([zero(tType)], c[1:1:end])
        Y2 .= hcat(zeros(tType, s), mu')

        U = Array{uType}(undef, s)
        U_ = Array{uType}(undef, s)
        L = Array{uType}(undef, s)
        L_ = Array{uType}(undef, s)
        F = Array{uType}(undef, s)
        for i in 1:s
            U[i] = zero(u0)
            U_[i] = zero(u0)
            L[i] = zero(u0)
            L_[i] = zero(u0)
            F[i] = zero(u0)
        end

        cache = tcache(
            f,
            p,
            abstol,
            reltol,
            U,
            U_,
            L,
            L_,
            F,
            Dmin,
            maxiters,
            maxtrials,
            step_number,
            initial_extrapolation,
            length_u,
            length_q,
            tf,
            fill(zero(tType), 2)
        )
    end

    #if saveat !== nothing
    #    error("saveat is not currently supported in this algorithm.")
    #end
    # Define save times
    save_times = [0]
    next_save_index = 2  # start saving after t0
    if saveat != nothing
        if isa(saveat, Number)
            save_everystep = false
            if tf > t0
                save_times = sort(union(collect(t0:saveat:tf), tf), order = Base.Forward)
            else
                save_times = sort(union(collect(tf:saveat:t0), t0), order = Base.Reverse)
            end
        elseif isa(saveat, AbstractVector)
            save_everystep = false
            if tf > t0
                save_times = sort(union(t0, saveat, tf), order = Base.Forward)  # include t0 and tf to ensure full span
            else
                save_times = sort(union(tf, saveat, t0), order = Base.Reverse)
            end
        else
            error("Unsupported saveat format")
        end
    end

    #   initialization output variables
    uu = uType[]
    tt = tType[]

    push!(uu, copy(u0))
    push!(tt, t0)

    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej = zero(u0)
    ux = copy(u0)  # to compute saveat
    indices = eachindex(uj)
    uj_ = copy(uj)
    tj_ = tj[1]

    cont = true
    error_warn = 0

    step_retcode = true

    while cont
        tj_ = tj[1]

        if saveat != nothing
            uj_ .= uj
        end

        step_number[] += 1
        step_retcode = step_fun(
            tj,
            uj,
            ej,
            dts,
            stats,
            coeffs,
            cache
        )

        if !step_retcode
            error_warn = 1
            cont = false
            break
        end

        if tj[1] == tf
            cont = false
            break
        end

        if save_everystep
            push!(tt, tj[1])
            push!(uu, copy(uj))
        elseif saveat != nothing
            while signdt * (tj[1] - save_times[next_save_index]) > 0
                #interpolation
                theta = abs((save_times[next_save_index] - tj_) / dts[2])
                Z2 .= [theta]

                if !simd
                    kappa .= (PolInterp(X2, Y2, Z2))
                    for k in indices
                        dUik = kappa[1] * L[1][k]
                        for js in 2:s
                            dUik = muladd(kappa[js], L[js][k], dUik)
                        end
                        ux[k] = uj_[k] + dUik
                    end
                else
                    kappa_ .= (PolInterp(X2, Y2, Z2))
                    kappa = vload(Vec{s, floatType}, kappa_, 1)
                    for k in indices
                        Lk = L[k]
                        ux[k] = uj_[k] + sum(Lk * kappa)
                    end
                end

                push!(tt, save_times[next_save_index])
                push!(uu, copy(ux))
                next_save_index += 1
            end
            if tj[1] == save_times[next_save_index]
                push!(tt, tj[1])
                push!(uu, copy(uj))
                next_save_index += 1
            end
        end
    end

    stats.naccept = step_number[]

    if error_warn != 0
        @warn("Error during the integration warn=$error_warn")
        sol = SciMLBase.build_solution(
            prob, alg, tt, uu, stats = stats, retcode = ReturnCode.Failure
        )

    else
        if saveat != nothing
            while signdt * (tj[1] - save_times[next_save_index]) > 0
                #interpolation
                theta = abs((save_times[next_save_index] - tj_) / dts[2])
                Z2 .= [theta]

                if !simd
                    kappa .= (PolInterp(X2, Y2, Z2))
                    for k in indices
                        dUik = kappa[1] * L[1][k]
                        for js in 2:s
                            dUik = muladd(kappa[js], L[js][k], dUik)
                        end
                        ux[k] = uj_[k] + dUik
                    end
                else
                    kappa_ .= (PolInterp(X2, Y2, Z2))
                    kappa = vload(Vec{s, floatType}, kappa_, 1)
                    for k in indices
                        Lk = L[k]
                        ux[k] = uj_[k] + sum(Lk * kappa)
                    end
                end

                push!(tt, save_times[next_save_index])
                push!(uu, copy(ux))
                next_save_index += 1
            end
        end

        push!(uu, copy(uj))
        push!(tt, tj[1])

        sol = SciMLBase.build_solution(
            prob, alg, tt, uu, stats = stats, retcode = ReturnCode.Success
        )
    end

    return (sol)
end
