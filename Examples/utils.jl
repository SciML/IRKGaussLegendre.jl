function run_many_NbodyIRKGL(
        alg, ddt0, tspan, u0, params, HAM; simdB = true, adaptiveB = true
    )
    nruns = 1

    u0_B = BigFloat.(u0)
    Gm_B = BigFloat.(params)

    cpus = similar(ddt0)
    iters = similar(ddt0)
    retcodes = [true for k in ddt0]
    MaxΔHglobal = [0.0 for i in ddt0]
    MaxΔHlocal = [0.0 for i in ddt0]

    sols = Array{IRKGaussLegendre.ODESolution}(undef, length(ddt0))
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(NbodyODE!, u0, tspan, Gm)
    H0 = HAM(u0_B, Gm_B)

    for i in 1:length(ddt0)
        print(",", ddt0[i])
        dt0 = ddt0[i]

        # save_everystep=true
        m0 = 1
        sols[i] = solve(prob, alg, dt = dt0, adaptive = false)
        if sols[i].retcode == ReturnCode.Success
            iters[i] = sols[i].stats.nfpiter
        else
            retcodes[i] = false
            iters[i] = Inf
        end

        m0 = max(1, div(Int64(ceil((tF - t0) / ddt0[i])), 1000))
        H = [HAM(BigFloat.(u), Gm_B) for u in sols[i].u]
        ΔH0 = @. Float64(abs(H / H0 - 1))
        H_lerr = @. Float64(abs((H[2:end] / H[1:(end - 1)]) - 1))
        MaxΔHglobal[i] = maximum(ΔH0)
        MaxΔHlocal[i] = maximum(H_lerr)

        # save_everystep=false
        solx = solve(prob, alg, dt = dt0, adaptive = false, save_everystep = false)
        if solx.retcode == ReturnCode.Success
            cpus[i] = 0.0
            for k in 1:nruns
                cpus[i] += @elapsed solve(
                    prob, alg, dt = dt0, adaptive = false, save_everystep = false
                )
            end
            cpus[i] = cpus[i] / nruns
        else
            cpus[i] = Inf
        end
    end

    return sols, retcodes, iters, cpus, MaxΔHglobal, MaxΔHlocal
end

function plots_NbodyIRKGL(title, ddt0, u0, params, HAM, sols, iters, cpus, MaxΔH)
    u0_B = BigFloat.(u0)
    Gm_B = BigFloat.(params)

    H0 = HAM(u0_B, Gm_B)
    MaxΔH = [0.0 for i in ddt0]

    pl1 = plot(
        ddt0, cpus, seriestype = :scatter, label = "",
        title = "CPU-time", xlabel = "dt", ylabel = "CPU"
    )
    pl2 = plot(
        ddt0, iters, seriestype = :scatter, label = "",
        title = "Iterations", xlabel = "dt", ylabel = "iter"
    )

    pl3 = plot(
        title = "Error in Ham", xlabel = "t ", ylabel = "log10(H/H0)",
        yscale = :log10, label = ""
    )

    for i in 1:length(ddt0)
        m0 = max(1, div(Int64(ceil((tF - t0) / ddt0[i])), 1000))
        ΔH0 = map(x -> HAM(BigFloat.(x), Gm_B), sols[i].u) ./ H0 .- 1
        pl3 = plot!(sols[i].t[2:m0:end], abs.(ΔH0[2:m0:end]), labels = "")
    end

    fig = plot(
        pl1, pl2, pl3, layout = (1, 3), size = (950, 300),
        plot_title = title, plot_titlevspan = 0.2
    )
    return fig
end
