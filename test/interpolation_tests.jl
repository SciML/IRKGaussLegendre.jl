using Test
using IRKGaussLegendre, SciMLBase, ForwardDiff

function lotka_volterra!(du, u, p, t)
    du[1] = 1.5 * u[1] - u[1] * u[2]
    du[2] = -3 * u[2] + u[1] * u[2]
    return nothing
end

function harmonic!(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[1]
    return nothing
end

exponential!(du, u, p, t) = (du .= u)

maxdiff(a, b) = maximum(maximum(abs.(x .- y)) for (x, y) in zip(a, b))

@testset "Dense output uses the collocation interpolant" begin
    prob = ODEProblem(lotka_volterra!, [1.0, 1.0], (0.0, 10.0))
    ts = 0.05:0.1:9.95
    ref = solve(prob, IRKGL16(), abstol = 1.0e-15, reltol = 1.0e-15, saveat = ts)

    @testset "simd = $simd" for simd in (false, true)
        sol = solve(prob, IRKGL16(simd = simd), abstol = 1.0e-12, reltol = 1.0e-12)
        sa = solve(prob, IRKGL16(simd = simd), abstol = 1.0e-12, reltol = 1.0e-12, saveat = ts)
        @test sol.dense
        @test maxdiff(sol.(ts), ref.u[2:(end - 1)]) < 1.0e-6
        @test maxdiff(sol.(ts), sa.u[2:(end - 1)]) < 1.0e-12
        @test maxdiff(sol(ts).u, sol.(ts)) == 0
        @test all(sol(sol.t[i]) == sol.u[i] for i in eachindex(sol.t))
        @test sol(3.3; idxs = 2) == sol(3.3)[2]
        @test sol(3.3; idxs = [2, 1]) == sol(3.3)[[2, 1]]
        @test_throws ErrorException sol(10.5)
    end

    @testset "backward integration" begin
        bprob = ODEProblem(lotka_volterra!, ref.u[end], (10.0, 0.0))
        sol = solve(bprob, IRKGL16(), abstol = 1.0e-12, reltol = 1.0e-12)
        @test maxdiff(sol.(ts), ref.u[2:(end - 1)]) < 1.0e-6
    end

    @testset "fixed step" begin
        sol = solve(prob, IRKGL16(), dt = 0.25, adaptive = false)
        sa = solve(prob, IRKGL16(), dt = 0.25, adaptive = false, saveat = ts)
        @test maxdiff(sol.(ts), sa.u[2:(end - 1)]) < 1.0e-12
        @test maxdiff(sol.(ts), ref.u[2:(end - 1)]) < 1.0e-6
    end

    @testset "BigFloat" begin
        bprob = ODEProblem(lotka_volterra!, BigFloat[1, 1], (big(0.0), big(10.0)))
        sol = solve(bprob, IRKGL16(), abstol = big(1.0e-20), reltol = big(1.0e-20))
        @test eltype(sol(big(3.3))) == BigFloat
        @test maxdiff(sol.(big.(ts)), ref.u[2:(end - 1)]) < 1.0e-6
    end

    @testset "no dense output when not saving every step" begin
        sol = solve(prob, IRKGL16(), save_everystep = false)
        @test !sol.dense
        sa = solve(prob, IRKGL16(), saveat = ts)
        @test !sa.dense
    end
end

@testset "Dense output of second_order_ode solutions" begin
    prob = ODEProblem(harmonic!, [1.0, 0.0], (0.0, 2.0))
    @testset "simd=$simd fseq=$fseq adaptive=$adaptive" for simd in (false, true),
            fseq in (false, true), adaptive in (false, true)

        sol = solve(
            prob, IRKGL16(; simd, fseq, second_order_ode = true);
            dt = 0.2, adaptive, abstol = 1.0e-12, reltol = 1.0e-12
        )
        ts = 0.05:0.1:1.95
        err = maximum(maximum(abs.(sol(t) .- [cos(t), -sin(t)])) for t in ts)
        @test err < 1.0e-7
        # The collocation polynomial of each step must meet the saved endpoint
        # values on both sides.
        for i in 2:length(sol.t)
            @test sol(prevfloat(sol.t[i])) ≈ sol.u[i] atol = 1.0e-10
        end
        for i in 1:(length(sol.t) - 1)
            @test sol(nextfloat(sol.t[i])) ≈ sol.u[i] atol = 1.0e-10
        end
        # Position components of Val{1} carry the derivative of the position
        # interpolant, which approximates the true velocity.
        for t in (0.3, 0.7, 1.3)
            @test sol(t, Val{1}) ≈ [-sin(t), -cos(t)] rtol = 1.0e-9
        end
    end
end

@testset "Val{1} differentiates the second-order position interpolant" begin
    poly7!(du, u, p, t) = (du[1] = u[2]; du[2] = t^7; nothing)
    @testset "polynomial forcing, simd=$simd" for simd in (false, true)
        sol = solve(
            ODEProblem(poly7!, [0.0, 0.0], (0.0, 1.0)),
            IRKGL16(; simd, second_order_ode = true); dt = 1.0, adaptive = false
        )
        for t in (0.13, 0.37, 0.63, 0.87)
            @test sol(t, Val{1}) ≈ ForwardDiff.derivative(t -> sol(t), t) rtol =
                1.0e-8 atol = 1.0e-12
        end
    end

    # A single oscillator step of size 2 makes the position interpolant's
    # derivative differ measurably from the velocity interpolant.
    @testset "oscillator h=2, simd=$simd backwards=$backwards" for simd in (false, true),
            backwards in (false, true)

        a, b = backwards ? (2.0, 0.0) : (0.0, 2.0)
        sol = solve(
            ODEProblem(harmonic!, [cos(a + 0.7), -sin(a + 0.7)], (a, b)),
            IRKGL16(; simd, second_order_ode = true);
            dt = 2.0, adaptive = false, maxiters = 1000
        )
        for t in a .+ (b - a) .* (0.13, 0.37, 0.63, 0.87)
            @test sol(t, Val{1}) ≈ ForwardDiff.derivative(t -> sol(t), t) rtol =
                1.0e-8 atol = 1.0e-12
        end
    end
end

@testset "First derivative of dense output" begin
    prob = ODEProblem(exponential!, [1.0], (0.0, 1.0))
    @testset "simd=$simd" for simd in (false, true)
        sol = solve(prob, IRKGL16(; simd), dt = 0.25, adaptive = false)
        for t in (0.1, 0.3, 0.6, 0.9)
            @test sol(t, Val{1})[1] ≈ exp(t) rtol = 1.0e-9
        end
        # Derivative at a saved endpoint uses the step polynomial, not the
        # saved value shortcut.
        @test sol(0.5, Val{1})[1] ≈ exp(0.5) rtol = 1.0e-9
        @test sol(0.5, Val{1}; continuity = :right)[1] ≈ exp(0.5) rtol = 1.0e-9
        # Vector-time, component selection, and in-place calls.
        @test sol([0.3, 0.7], Val{1}).u ≈ [[exp(0.3)], [exp(0.7)]] rtol = 1.0e-9
        @test sol(0.3, Val{1}; idxs = 1) ≈ exp(0.3) rtol = 1.0e-9
        out = zeros(1)
        sol(out, 0.3, Val{1})
        @test out[1] ≈ exp(0.3) rtol = 1.0e-9
        @test_throws ArgumentError sol(0.3, Val{2})
    end
end

@testset "AD through interpolation time" begin
    prob = ODEProblem(exponential!, [1.0], (0.0, 1.0))
    sol = solve(prob, IRKGL16(), dt = 0.25, adaptive = false)
    @test ForwardDiff.derivative(t -> sol(t)[1], 0.3) ≈ exp(0.3) rtol = 1.0e-9
    # A dual-valued query at a saved time still flows through the polynomial.
    @test ForwardDiff.derivative(t -> sol(t)[1], 0.5) ≈ exp(0.5) rtol = 1.0e-9
end

@testset "AD honors `continuity` at saved knots" begin
    # u' = t^9 makes the one-sided interpolant derivatives differ measurably at
    # the interior saved points, so a dual-numbered query must pick the same
    # step as a Float64 `Val{1}` query.
    poly9!(du, u, p, t) = (du[1] = t^9; nothing)
    @testset "simd=$simd backwards=$backwards" for simd in (false, true),
            backwards in (false, true)

        a, b = backwards ? (2.0, 0.0) : (0.0, 2.0)
        sol = solve(
            ODEProblem(poly9!, [0.0], (a, b)), IRKGL16(; simd);
            dt = 1.0, adaptive = false
        )
        for tk in sol.t, continuity in (:left, :right)
            @test ForwardDiff.derivative(t -> sol(t; continuity)[1], tk) ≈
                sol(tk, Val{1}; continuity)[1] rtol = 1.0e-11 atol = 1.0e-12
            @test ForwardDiff.derivative(t -> sol(t, Val{1}; continuity)[1], tk) ≈
                ForwardDiff.derivative(
                t -> ForwardDiff.derivative(w -> sol(w; continuity)[1], t), tk
            ) rtol = 1.0e-11 atol = 1.0e-12
        end
        tdir = sign(sol.t[end] - sol.t[1])
        @test_throws ErrorException ForwardDiff.derivative(
            t -> sol(t)[1], sol.t[end] + tdir
        )
    end
end
