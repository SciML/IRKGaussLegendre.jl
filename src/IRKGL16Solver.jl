
struct tcoeffs{tType}
    mu::Array{tType, 2}
    c::Array{tType, 1}
    b::Array{tType, 1}
    nu::Array{tType, 2}
    alpha::Array{tType, 1}
    X::Array{tType, 1}
    Y::Array{tType, 2}
    Z::Array{tType, 1}
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
    second_order_ode,
    mstep,
    simd,
    maxtrials,
    initial_extrapolation,
    threading
} <: SciMLBase.AbstractODEAlgorithm end
struct IRKGL16{
    second_order_ode,
    mstep,
    simd,
    maxtrials,
    initial_extrapolation,
    threading
} <: IRKAlgorithm{
    second_order_ode,
    mstep,
    simd,
    maxtrials,
    initial_extrapolation,
    threading
} end
function IRKGL16(;
        second_order_ode = false,
        mstep = 1,
        simd = false,
        maxtrials = 5,
        initial_extrapolation = true,
        threading = false)
    IRKGL16{
        second_order_ode,
        mstep,
        simd,
        maxtrials,
        initial_extrapolation,
        threading
    }()
end

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem{
            uType,
            tspanType,
            isinplace
        },
        alg::IRKGL16{
            second_order_ode,
            mstep,
            simd,
            maxtrials,
            initial_extrapolation,
            threading
        },
        args...;
        dt = zero(eltype(tspanType)),
        maxiters = 100,
        save_everystep = true,
        adaptive = true,
        reltol = eltype(tspanType)(1e-6),
        abstol = eltype(tspanType)(1e-6),
        kwargs...) where {
        uType,
        tspanType,
        isinplace,
        second_order_ode,
        mstep,
        simd,
        maxtrials,
        initial_extrapolation,
        threading
}
    checks = true

    s = 8
    stats = SciMLBase.DEStats(0)
    #stats = DiffEqBase.DEStats(0)
    stats.nf = 0
    stats.nf2 = 0
    stats.nfpiter = 0
    stats.naccept = 0
    stats.nreject = 0

    if (prob.f isa DynamicalODEFunction)
        @unpack tspan, p = prob
        #        f = prob.f
        f = SciMLBase.unwrapped_f(prob.f)
        f1 = prob.f.f1
        f2 = prob.f.f2
        r0 = prob.u0.x[1]
        v0 = prob.u0.x[2]
        u0 = ArrayPartition(r0, v0)
    elseif (prob.f isa ODEFunction)
        @unpack u0, tspan, p = prob
        f = SciMLBase.unwrapped_f(prob.f)
    else
        @warn("Error: incorrect ODEFunction")
        sol = SciMLBase.build_solution(prob, alg, [], [], retcode = ReturnCode.Failure)
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
                step_fun = IRKGLstep_SIMD_fixed!
            else
                step_fun = IRKNGLstep_SIMD_fixed_simpl!
            end
        else
            if !second_order_ode
                step_fun = IRKstep_SIMD_adaptive!
            else
                step_fun = IRKNGLstep_SIMD_adaptive_simpl!
            end
        end

    else
        if (threading == true && Threads.nthreads() > 1)
            #    step_fun=IRKStep_par!
            if (prob.f isa ODEFunction)
                if adaptive
                    step_fun = IRKstep_par_adaptive!
                else
                    step_fun = IRKstep_par_fixed!
                end
            else # (prob.f<:DynamicalODEFunction)
                #
                # @warn("The DynamicalODEProblem problem type will be removed in a future version. Please, code it as ODEProblem with the second_order_ode=true keyword argument instead.")
                #
                if adaptive
                    step_fun = IRKstepDynODE_par_adaptive!
                else
                    step_fun = IRKstepDynODE_par_fixed!
                end
            end
            #
            #
        else
            #    step_fun=IRKStep_seq!
            if (prob.f isa ODEFunction)
                if adaptive
                    #
                    if !second_order_ode
                        step_fun = IRKstep_adaptive!
                    else
                        step_fun = IRKNGLstep_adaptive_simpl!
                    end
                else
                    if !second_order_ode
                        step_fun = IRKstep_fixed!
                    else
                        step_fun = IRKNGLstep_fixed_simpl!
                    end
                end
            else # (prob.f<:DynamicalODEFunction)
                #
                #@warn("The DynamicalODEProblem problem type will be removed in a future version. Please, code it as ODEProblem with the second_order_ode=true keyword argument instead.")
                #
                if adaptive
                    step_fun = IRKstepDynODE_adaptive!
                else
                    step_fun = IRKstepDynODE_fixed!
                end
            end
        end
    end

    abstol = convert(tType, abstol)
    reltol = convert(tType, reltol)

    if (dt == 0)
        d0 = MyNorm(u0, abstol, reltol)
        du0 = similar(u0)
        if (prob.f isa DynamicalODEFunction)
            f1(du0.x[1], u0.x[1], u0.x[2], p, t0)
            f2(du0.x[2], u0.x[1], u0.x[2], p, t0)
        else # ODEFunction
            f(du0, u0, p, t0)
        end
        d1 = MyNorm(du0, abstol, reltol)
        if (d0 < 1e-5 || d1 < 1e-5)
            dt = convert(tType, 1e-6)
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

    if simd
        coeffs = tcoeffs{tType}(zeros(s, s), zeros(s), zeros(s), zeros(s, s),
            zeros(s), zeros(s + 1), zeros(s, s + 1), zeros(s))
        mu_ = coeffs.mu
        c_ = coeffs.c
        b_ = coeffs.b
        nu_ = coeffs.nu
        alpha_ = coeffs.alpha
        X = coeffs.X
        Y = coeffs.Y
        Z = coeffs.Z

        EstimateCoeffs!(alpha_, tType)

        GaussLegendreCoefficients!(mu_, c_, b_, tType)
        #ExtrapolationCoefficients!(nu_, mu_, c_,  signdt*dt, signdt*dtprev, tType)
        X .= vcat(-c_[end:-1:1], [0])
        Y .= hcat(mu_, zeros(tType, s))
        Z .= c_
        if dtprev == 0
            nu_ .= zeros(tType, s, s)
        else
            nu_ .= -(PolInterp(sdt * X, Y, sdt * Z))'
        end

        dims = size(u0)

        c = vload(Vec{s, floatType}, c_, 1)
        b = vload(Vec{s, floatType}, b_, 1)
        nu = VecArray{s, floatType, 2}(nu_)
        mu = VecArray{s, floatType, 2}(mu_)
        alpha = vload(Vec{s, floatType}, alpha_, 1)

        zz = zeros(Float64, s, dims...)
        U = VecArray{s, Float64, length(dims) + 1}(zz)
        U_ = deepcopy(U)
        L = deepcopy(U)
        L_ = deepcopy(U)
        F = deepcopy(U)

        coeffs = tcoeffs_SIMD(b, c, mu, nu, alpha, X, Y, Z, nu_)

        cache = IRKGL_SIMD_Cache(f, p, abstol, reltol,
            U, U_, L, L_, F,
            Dmin, maxiters, maxtrials, step_number,
            initial_extrapolation, length_u, length_q, tf,
            fill(zero(tType), 2))

    else
        coeffs = tcoeffs{tType}(zeros(s, s), zeros(s), zeros(s), zeros(s, s),
            zeros(s), zeros(s + 1), zeros(s, s + 1), zeros(s))
        @unpack mu, c, b, nu, alpha, X, Y, Z = coeffs
        EstimateCoeffs!(alpha, tType)

        GaussLegendreCoefficients!(mu, c, b, tType)
        #ExtrapolationCoefficients!(nu, mu, c,  signdt*dt, signdt*dtprev, tType)
        X .= vcat(-c[end:-1:1], [0])
        Y .= hcat(mu, zeros(tType, s))
        Z .= c
        if dtprev == 0
            nu .= zeros(tType, s, s)
        else
            nu .= -(PolInterp(sdt * X, Y, sdt * Z))'
        end

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

        cache = tcache(f,
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
            fill(zero(tType), 2))
    end

    #   initialization output variables
    uu = uType[]
    tt = tType[]

    push!(uu, copy(u0))
    push!(tt, t0)

    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej = zero(u0)

    cont = true
    error_warn = 0

    step_retcode = true

    while cont
        @inbounds begin
            for i in 1:mstep
                step_number[] += 1
                step_retcode = step_fun(tj,
                    uj,
                    ej,
                    dts,
                    stats,
                    coeffs,
                    cache)

                if !step_retcode
                    error_warn = 1
                    cont = false
                    break
                end

                if tj[1] == tf
                    cont = false
                    break
                end
            end
        end

        if save_everystep
            push!(tt, tj[1])
            push!(uu, copy(uj))
        end
    end

    stats.naccept = step_number[]

    if error_warn != 0
        @warn("Error during the integration warn=$error_warn")
        sol = SciMLBase.build_solution(
            prob, alg, tt, uu, stats = stats, retcode = ReturnCode.Failure)

    else
        if !save_everystep
            push!(uu, copy(uj))
            push!(tt, tj[1])
        end

        sol = SciMLBase.build_solution(
            prob, alg, tt, uu, stats = stats, retcode = ReturnCode.Success)
    end

    return (sol)
end
