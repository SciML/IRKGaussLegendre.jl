using IRKGaussLegendre
using IRKGaussLegendre: IRKstep_fixed!, IRKstep_adaptive!, tcoeffs, tcache,
    GaussLegendreCoefficients!, EstimateCoeffs!, PolInterp!, MyNorm
using Test
using AllocCheck

# Test ODE function that avoids global variable boxing
struct TestParams
    g::Float64
    L::Float64
end

function test_pendulum!(du, u, p::TestParams, t)
    du[1] = u[2]
    du[2] = -(p.g / p.L) * sin(u[1])
    return nothing
end

@testset "AllocCheck - Core Functions" begin
    # Setup
    s = 8
    u0 = [0.0, pi / 2]
    params = TestParams(9.81, 1.0)

    # Create coefficients
    coeffs = tcoeffs{Float64}(
        zeros(s, s), zeros(s), zeros(s), zeros(s, s),
        zeros(s), zeros(s + 1), zeros(s, s + 1), zeros(s),
        zeros(s), zeros(s + 1), zeros(s, s + 1), zeros(1)
    )

    # Initialize coefficients
    GaussLegendreCoefficients!(coeffs.mu, coeffs.c, coeffs.b, Float64)
    EstimateCoeffs!(coeffs.alpha, Float64)

    # Create cache
    U = [zero(u0) for _ in 1:s]
    U_ = [zero(u0) for _ in 1:s]
    L_ = [zero(u0) for _ in 1:s]
    L__ = [zero(u0) for _ in 1:s]
    F = [zero(u0) for _ in 1:s]
    Dmin = fill(Inf, 2)
    step_number = Array{Int64, 0}(undef)
    step_number[] = 2  # Not first step to avoid extra iterations

    cache = tcache(
        test_pendulum!, params,
        1.0e-8, 1.0e-8,
        U, U_, L_, L__, F,
        Dmin, 100, 5, step_number,
        true, length(u0), div(length(u0), 2), 10.0,
        fill(0.0, 2)
    )

    # Test vectors
    ttj = [0.0, 0.0]
    uj = copy(u0)
    ej = zero(u0)
    dts = [0.1, 0.1, 1.0]  # dt == dtprev to avoid PolInterp allocation
    stats = SciMLBase.DEStats(0)

    @testset "PolInterp! zero allocations" begin
        X = zeros(s + 1)
        Y = zeros(s, s + 1)
        Z = zeros(s)
        pz = zeros(s, s)

        X .= [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.0]
        Y .= rand(s, s + 1)
        Z .= [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

        # Warm up
        PolInterp!(pz, X, Y, Z)

        # Test zero allocations
        alloc = @allocated PolInterp!(pz, X, Y, Z)
        @test alloc == 0
    end

    @testset "MyNorm zero allocations" begin
        # Warm up
        MyNorm(u0, 1.0e-8, 1.0e-8)

        # Test zero allocations
        alloc = @allocated MyNorm(u0, 1.0e-8, 1.0e-8)
        @test alloc == 0
    end

    @testset "IRKstep_fixed! zero allocations" begin
        # Warm up the step function
        for _ in 1:10
            ttj .= [0.0, 0.0]
            uj .= copy(u0)
            ej .= zero(u0)
            dts .= [0.1, 0.1, 1.0]
            step_number[] = 2
            IRKstep_fixed!(ttj, uj, ej, dts, stats, coeffs, cache)
        end

        # Test zero allocations
        ttj .= [0.0, 0.0]
        uj .= copy(u0)
        ej .= zero(u0)
        dts .= [0.1, 0.1, 1.0]
        step_number[] = 2

        alloc = @allocated IRKstep_fixed!(ttj, uj, ej, dts, stats, coeffs, cache)
        @test alloc == 0
    end

    @testset "IRKstep_adaptive! zero allocations" begin
        # Warm up the step function
        for _ in 1:10
            ttj .= [0.0, 0.0]
            uj .= copy(u0)
            ej .= zero(u0)
            dts .= [0.1, 0.1, 1.0]  # dt == dtprev to avoid PolInterp allocation
            step_number[] = 2
            IRKstep_adaptive!(ttj, uj, ej, dts, stats, coeffs, cache)
        end

        # Test zero allocations
        ttj .= [0.0, 0.0]
        uj .= copy(u0)
        ej .= zero(u0)
        dts .= [0.1, 0.1, 1.0]
        step_number[] = 2

        alloc = @allocated IRKstep_adaptive!(ttj, uj, ej, dts, stats, coeffs, cache)
        @test alloc == 0
    end
end

@testset "AllocCheck - Full Solver Allocation Bounds" begin
    # Test that the solver's allocation overhead is reasonable
    # Most allocations should be for solution storage

    u0 = [0.0, pi / 2]
    params = TestParams(9.81, 1.0)
    tspan = (0.0, 10.0)
    prob = ODEProblem(test_pendulum!, u0, tspan, params)

    # Warm up thoroughly
    for _ in 1:20
        solve(prob, IRKGL16(), reltol = 1.0e-8, abstol = 1.0e-8)
    end

    # Measure allocations
    alloc = @allocated sol = solve(prob, IRKGL16(), reltol = 1.0e-8, abstol = 1.0e-8)
    nsteps = length(sol.t)

    # Check that per-step overhead is reasonable
    # We expect:
    # - Solution storage: ~24 bytes per step (1 Float64 time + 2 Float64 state)
    # - Some overhead for vectors
    # We allow up to 5KB per step as a generous bound
    per_step = alloc ÷ nsteps
    @test per_step < 5000  # Less than 5KB per step

    # Total allocation should be bounded
    @test alloc < 150_000  # Less than 150KB total for this problem
end
