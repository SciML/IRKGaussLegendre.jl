
struct tcoeffs{T}
    mu::Array{T, 2}
    hc::Array{T, 1}
    hb::Array{T, 1}
    nu::Array{T, 2}
    alpha::Array{T, 1}
end

struct tcache{uType, realuType, realuiType, uLowType, low_prec_type}
    U::Array{uType, 1}
    Uz::Array{uType, 1}
    L::Array{uType, 1}
    Lz::Array{uType, 1}
    F::Array{uType, 1}
    #    Dmin::Array{uType, 1}
    Dmin::Array{realuType, 1}
    Eval::Array{Bool, 1}
    DY::Array{realuiType, 1}
    rejects::Array{Int64, 1}
    nfcn::Array{Int64, 1}
    lambdas::Array{realuiType, 1}
    nrmdigits::Array{Int64, 0}
end

struct tcacheMix{uType, realuType, realuiType, uLowType, low_prec_type}
    U::Array{uType, 1}
    Uz::Array{uType, 1}
    L::Array{uType, 1}
    Lz::Array{uType, 1}
    F::Array{uType, 1}
    #    Dmin::Array{uType, 1}
    Dmin::Array{realuType, 1}
    Eval::Array{Bool, 1}
    DY::Array{realuiType, 1}
    rejects::Array{Int64, 1}
    nfcn::Array{Int64, 1}
    lambdas::Array{realuiType, 1}
    Ulow::Array{uLowType, 1}
    DU::Array{uLowType, 1}
    DF::Array{uLowType, 1}
    DL::Array{uLowType, 1}
    delta::Array{uLowType, 1}
    Fa::Array{uLowType, 1}
    Fb::Array{uLowType, 1}
    Plow::Array{low_prec_type, 1}
    normU::Array{low_prec_type, 1}
    lhb::Array{low_prec_type, 1}
    lmu::Array{low_prec_type, 2}
    nrmdigits::Array{Int64, 0}
end

abstract type IRKAlgorithm{
    mstep,
    maxtrials,
    initial_interp,
    myoutputs,
    threading,
    mixed_precision,
    low_prec_type,
    nrmbits,
} <: OrdinaryDiffEqAlgorithm end
struct IRKGL16{
    mstep,
    maxtrials,
    initial_interp,
    myoutputs,
    threading,
    mixed_precision,
    low_prec_type,
    nrmbits,
} <: IRKAlgorithm{
    mstep,
    maxtrials,
    initial_interp,
    myoutputs,
    threading,
    mixed_precision,
    low_prec_type,
    nrmbits,
} end
function IRKGL16(;
        mstep = 1,
        maxtrials = 5,
        initial_interp = true,
        myoutputs = false,
        threading = false,
        mixed_precision = false,
        low_prec_type = Float64,
        nrmbits = 6)
    IRKGL16{
        mstep,
        maxtrials,
        initial_interp,
        myoutputs,
        threading,
        mixed_precision,
        low_prec_type,
        nrmbits,
    }()
end

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{
            uType,
            tspanType,
            isinplace,
        },
        alg::IRKGL16{
            mstep,
            maxtrials,
            initial_interp,
            myoutputs,
            threading,
            mixed_precision,
            low_prec_type,
            nrmbits,
        },
        args...;
        dt = 0.0,
        maxiters = 100,
        save_everystep = true,
        adaptive = true,
        reltol = 1e-6,
        abstol = 1e-6,
        kwargs...) where {
        uType,
        tspanType,
        isinplace,
        mstep,
        maxtrials,
        initial_interp,
        myoutputs,
        threading,
        mixed_precision,
        low_prec_type,
        nrmbits,
    }
    s = 8
    stats = DiffEqBase.Stats(0)

    if (prob.f isa DynamicalODEFunction)
        @unpack tspan, p = prob
        f1 = prob.f.f1
        f2 = prob.f.f2
        r0 = prob.u0.x[1]
        v0 = prob.u0.x[2]
        u0 = ArrayPartition(r0, v0)
    elseif (prob.f isa ODEFunction)
        @unpack f, u0, tspan, p = prob
    else
        println("Error: incorrect ODEFunction")
        sol = DiffEqBase.build_solution(prob, alg, [], [], retcode = ReturnCode.Failure)
        return (sol)
    end

    signdt = sign(tspan[2] - tspan[1])

    t0 = prob.tspan[1]
    tf = prob.tspan[2]

    tType = eltype(tspanType)
    uiType = eltype(u0)
    realuType = typeof(real(u0))
    realuiType = real(uiType)

    nrmdig = Array{Int64, 0}(undef)
    if (nrmbits > 0)
        nrmdig[] = 2^nrmbits
    else
        nrmdig[] = 0
    end

    if uiType <: Complex
        lu0 = convert.(Complex{low_prec_type}, u0)
        uLowType = typeof(lu0)
    else
        lu0 = convert.(low_prec_type, u0)
        uLowType = typeof(lu0)
    end

    coeffs = tcoeffs{realuiType}(zeros(s, s), zeros(s), zeros(s), zeros(s, s), zeros(s))

    @unpack mu, hc, hb, nu, alpha = coeffs

    Treltol = convert(realuiType, reltol)
    Tabstol = convert(realuiType, abstol)

    if (dt == 0)
        d0 = MyNorm(u0, Tabstol, Treltol)
        du0 = similar(u0)
        if (prob.f isa DynamicalODEFunction)
            f1(du0.x[1], u0.x[1], u0.x[2], p, t0)
            f2(du0.x[2], u0.x[1], u0.x[2], p, t0)
        else # ODEFunction
            f(du0, u0, p, t0)
        end
        d1 = MyNorm(du0, Tabstol, Treltol)
        if (d0 < 1e-5 || d1 < 1e-5)
            dt = convert(tType, 1e-6)
        else
            dt = convert(tType, 0.01 * (d0 / d1))
        end
    end

    #    dt = min(dt, tf - t0)
    dt = min(abs(dt), abs(tf - t0))

    EstimateCoeffs!(alpha, realuiType)
    MuCoefficients!(mu, realuiType)

    dts = Array{tType}(undef, 1)

    if (adaptive == false)
        dtprev = dt
    else
        dtprev = zero(tType)
    end

    dts = [dt, dtprev, signdt]
    #    signdt = sign(dt)
    HCoefficients!(mu, hc, hb, nu, signdt * dt, signdt * dtprev, realuiType)

    #   m: output saved at every m steps
    #   n: Number of macro-steps  (Output is saved for n+1 time values)

    if (save_everystep == false)
        m = 1
        n = Inf
    else
        m = Int64(mstep)
        n = Inf
    end

    U1 = Array{uType}(undef, s)
    U2 = Array{uType}(undef, s)
    U3 = Array{uType}(undef, s)
    U4 = Array{uType}(undef, s)
    U5 = Array{uType}(undef, s)
    #    U6 = Array{uType}(undef, s)
    U6 = Array{realuType}(undef, s)
    for i in 1:s
        U1[i] = zero(u0)
        U2[i] = zero(u0)
        U3[i] = zero(u0)
        U4[i] = zero(u0)
        U5[i] = zero(u0)
        #        U6[i] = zero(u0)
        U6[i] = zero(real(u0))
    end

    if (mixed_precision == true && prob.f isa ODEFunction)
        U11 = Array{uLowType}(undef, s)
        U12 = Array{uLowType}(undef, s)
        U13 = Array{uLowType}(undef, s)
        U14 = Array{uLowType}(undef, s)
        U15 = Array{uLowType}(undef, s)
        U16 = Array{uLowType}(undef, s)
        U17 = Array{uLowType}(undef, s)
        for i in 1:s
            U11[i] = zero(lu0)
            U12[i] = zero(lu0)
            U13[i] = zero(lu0)
            U14[i] = zero(lu0)
            U15[i] = zero(lu0)
            U16[i] = zero(lu0)
            U17[i] = zero(lu0)
        end

        lmu = convert.(low_prec_type, mu)
        lhb = convert.(low_prec_type, hb)

        if (typeof(p) != SciMLBase.NullParameters && p != nothing)
            Plow = convert.(low_prec_type, p)
        else
            Plow = []
        end

        cache = tcacheMix{uType, realuType, realuiType, uLowType, low_prec_type}(U1,
            U2,
            U3,
            U4,
            U5,
            U6,
            fill(true, s),
            fill(zero(realuiType), s),
            [0],
            [0, 0],
            fill(zero(realuiType), 2),
            U11,
            U12,
            U13,
            U14,
            U15,
            U16,
            U17,
            Plow,
            fill(zero(low_prec_type),
                s),
            lhb,
            lmu,
            nrmdig)

        @unpack U,
        Uz,
        L,
        Lz,
        F,
        Dmin,
        Eval,
        DY,
        rejects,
        nfcn,
        lambdas,
        Ulow,
        DU,
        DF,
        DL,
        delta,
        Fa,
        Fb,
        Plow,
        normU,
        lhb,
        lmu,
        nrmdigits = cache

    else
        cache = tcache{uType, realuType, realuiType, uLowType, low_prec_type}(U1,
            U2,
            U3,
            U4,
            U5,
            U6,
            fill(true, s),
            fill(zero(realuiType), s),
            [0],
            [0, 0],
            fill(zero(realuiType), 2),
            nrmdig)
        @unpack U, Uz, L, Lz, F, Dmin, Eval, DY, rejects, nfcn, lambdas, nrmdigits = cache
    end

    #   initialization output variables
    uu = uType[]
    tt = tType[]
    iters = Int[]
    steps = Int[]

    push!(uu, u0)
    push!(tt, t0)
    push!(iters, 0)
    push!(steps, 0)
    tj = [t0, zero(t0)]
    uj = copy(u0)
    ej = zero(u0)

    j = 0
    cont = true

    if (threading == true && Threads.nthreads() > 1)
        while cont
            tit = 0
            it = 0
            k = 0

            @inbounds begin
                for i in 1:m
                    j += 1
                    k += 1
                    (status, it) = IRKStep_par!(s,
                        j,
                        tj,
                        tf,
                        uj,
                        ej,
                        prob,
                        dts,
                        coeffs,
                        cache,
                        maxiters,
                        maxtrials,
                        initial_interp,
                        Tabstol,
                        Treltol,
                        adaptive,
                        threading,
                        mixed_precision,
                        low_prec_type)

                    if (status == "Failure")
                        #                    println("Fail")
                        sol = DiffEqBase.build_solution(prob, alg, tt, uu,
                            retcode = ReturnCode.Failure)
                        return (sol)
                    end
                    tit += it

                    if (dts[1] == 0)
                        break
                    end
                end
            end

            #            cont = (signdt * (tj[1] + tj[2]) < signdt * tf) && (j < n * m)

            if signdt == 1
                cont = ((tj[1] + tj[2]) < tf) && (j < n * m)
            else
                cont = ((tj[1] + tj[2]) > tf) && (j < n * m)
            end

            if (save_everystep == true) || (cont == false)
                push!(tt, tj[1] + tj[2])
                push!(uu, uj + ej)

                if (myoutputs == true)
                    push!(iters, convert(Int64, round(tit / k)))
                    #                    push!(steps, dts[2])
                    push!(steps, signdt * dts[2])
                end
            end
        end

    else # threading==false
        while cont
            tit = 0
            it = 0
            k = 0

            @inbounds begin
                for i in 1:m
                    j += 1
                    k += 1
                    (status, it) = IRKStep_seq!(s,
                        j,
                        tj,
                        tf,
                        uj,
                        ej,
                        prob,
                        dts,
                        coeffs,
                        cache,
                        maxiters,
                        maxtrials,
                        initial_interp,
                        Tabstol,
                        Treltol,
                        adaptive,
                        threading,
                        mixed_precision,
                        low_prec_type)

                    if (status == "Failure")
                        #                    println("Fail")
                        sol = DiffEqBase.build_solution(prob, alg, tt, uu,
                            retcode = ReturnCode.Failure)
                        return (sol)
                    end
                    tit += it

                    if (dts[1] == 0)
                        break
                    end
                end
            end

            #            cont = (signdt * (tj[1] + tj[2]) < signdt * tf) && (j < n * m)

            if signdt == 1
                cont = ((tj[1] + tj[2]) < tf) && (j < n * m)
            else
                cont = ((tj[1] + tj[2]) > tf) && (j < n * m)
            end

            if (save_everystep == true) || (cont == false)
                push!(tt, tj[1] + tj[2])
                push!(uu, uj + ej)

                if (myoutputs == true)
                    push!(iters, convert(Int64, round(tit / k)))
                    push!(steps, signdt * dts[2])
                end
            end
        end
    end

    sol = DiffEqBase.build_solution(prob, alg, tt, uu, stats = stats,
        retcode = ReturnCode.Success)

    sol.stats.nf = nfcn[1]
    sol.stats.nf2 = nfcn[2]
    sol.stats.nreject = rejects[1]
    sol.stats.naccept = j

    if (myoutputs == true)
        return (sol, iters, steps)
    else
        return (sol)
    end
end

function IRKStep_seq!(s,
        j,
        tj,
        tf,
        uj,
        ej,
        prob,
        dts,
        coeffs,
        cache,
        maxiters,
        maxtrials,
        initial_interp,
        abstol,
        reltol,
        adaptive,
        threading,
        mixed_precision,
        low_prec_type)
    if (prob.f isa ODEFunction)
        if (adaptive == true)
            if (mixed_precision == true)
                (status, it) = IRKstep_adaptive_Mix!(s,
                    j,
                    tj,
                    tf,
                    uj,
                    ej,
                    prob,
                    dts,
                    coeffs,
                    cache,
                    maxiters,
                    maxtrials,
                    initial_interp,
                    abstol,
                    reltol,
                    adaptive,
                    threading,
                    mixed_precision,
                    low_prec_type)
            else
                (status, it) = IRKstep_adaptive!(s,
                    j,
                    tj,
                    tf,
                    uj,
                    ej,
                    prob,
                    dts,
                    coeffs,
                    cache,
                    maxiters,
                    maxtrials,
                    initial_interp,
                    abstol,
                    reltol,
                    adaptive,
                    threading)
            end
        else
            if (mixed_precision == true)
                (status, it) = IRKstep_fixed_Mix!(s,
                    j,
                    tj,
                    tf,
                    uj,
                    ej,
                    prob,
                    dts,
                    coeffs,
                    cache,
                    maxiters,
                    initial_interp,
                    abstol,
                    reltol,
                    adaptive,
                    threading,
                    mixed_precision,
                    low_prec_type)
            else
                (status, it) = IRKstep_fixed!(s,
                    j,
                    tj,
                    tf,
                    uj,
                    ej,
                    prob,
                    dts,
                    coeffs,
                    cache,
                    maxiters,
                    initial_interp,
                    abstol,
                    reltol,
                    adaptive,
                    threading)
            end
        end

    else  # (typeof(prob.f<:DynamicalODEFunction))
        if (adaptive == true)
            (status, it) = IRKstepDynODE_adaptive!(s,
                j,
                tj,
                tf,
                uj,
                ej,
                prob,
                dts,
                coeffs,
                cache,
                maxiters,
                maxtrials,
                initial_interp,
                abstol,
                reltol,
                adaptive,
                threading)
        else
            (status, it) = IRKstepDynODE_fixed!(s,
                j,
                tj,
                tf,
                uj,
                ej,
                prob,
                dts,
                coeffs,
                cache,
                maxiters,
                initial_interp,
                abstol,
                reltol,
                adaptive,
                threading)
        end
    end

    return (status, it)
end

function IRKStep_par!(s,
        j,
        tj,
        tf,
        uj,
        ej,
        prob,
        dts,
        coeffs,
        cache,
        maxiters,
        maxtrials,
        initial_interp,
        abstol,
        reltol,
        adaptive,
        threading,
        mixed_precision,
        low_prec_type)
    if (prob.f isa ODEFunction)
        if (adaptive == true)
            if (mixed_precision == true)
                (status, it) = IRKstep_par_adaptive_Mix!(s,
                    j,
                    tj,
                    tf,
                    uj,
                    ej,
                    prob,
                    dts,
                    coeffs,
                    cache,
                    maxiters,
                    maxtrials,
                    initial_interp,
                    abstol,
                    reltol,
                    adaptive,
                    threading,
                    mixed_precision,
                    low_prec_type)
            else
                (status, it) = IRKstep_par_adaptive!(s,
                    j,
                    tj,
                    tf,
                    uj,
                    ej,
                    prob,
                    dts,
                    coeffs,
                    cache,
                    maxiters,
                    maxtrials,
                    initial_interp,
                    abstol,
                    reltol,
                    adaptive,
                    threading)
            end
        else
            if (mixed_precision == true)
                (status, it) = IRKstep_par_fixed_Mix!(s,
                    j,
                    tj,
                    tf,
                    uj,
                    ej,
                    prob,
                    dts,
                    coeffs,
                    cache,
                    maxiters,
                    initial_interp,
                    abstol,
                    reltol,
                    adaptive,
                    threading,
                    mixed_precision,
                    low_prec_type)
            else
                (status, it) = IRKstep_par_fixed!(s,
                    j,
                    tj,
                    tf,
                    uj,
                    ej,
                    prob,
                    dts,
                    coeffs,
                    cache,
                    maxiters,
                    initial_interp,
                    abstol,
                    reltol,
                    adaptive,
                    threading)
            end
        end

    else  # (typeof(prob.f<:DynamicalODEFunction))
        if (adaptive == true)
            (status, it) = IRKstepDynODE_par_adaptive!(s,
                j,
                tj,
                tf,
                uj,
                ej,
                prob,
                dts,
                coeffs,
                cache,
                maxiters,
                maxtrials,
                initial_interp,
                abstol,
                reltol,
                adaptive,
                threading)
        else
            (status, it) = IRKstepDynODE_par_fixed!(s,
                j,
                tj,
                tf,
                uj,
                ej,
                prob,
                dts,
                coeffs,
                cache,
                maxiters,
                initial_interp,
                abstol,
                reltol,
                adaptive,
                threading)
        end
    end

    return (status, it)
end
