using IRKGaussLegendre, BenchmarkTools

const SUITE = BenchmarkGroup()

# Stiff-ish test problem
function rober!(du, u, p, t)
    du[1] = -0.04u[1] + 1.0e4 * u[2] * u[3]
    du[2] = 0.04u[1] - 1.0e4 * u[2] * u[3] - 3.0e7 * u[2]^2
    du[3] = 3.0e7 * u[2]^2
    return nothing
end
prob_stiff = ODEProblem(rober!, [1.0, 0.0, 0.0], (0.0, 1.0e4))

# Smooth non-stiff problem
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
    return nothing
end
prob = ODEProblem(lorenz!, [1.0, 0.0, 0.0], (0.0, 20.0))

# =============================================================================
# IRKGL16 solves
# =============================================================================

SUITE["solve"] = BenchmarkGroup()

SUITE["solve"]["adaptive_nonstiff"] = @benchmarkable solve(
    $prob, IRKGL16(); reltol = 1.0e-8, abstol = 1.0e-8
)
SUITE["solve"]["adaptive_stiff"] = @benchmarkable solve(
    $prob_stiff, IRKGL16(); reltol = 1.0e-6, abstol = 1.0e-6, saveat = 1.0e2
)
SUITE["solve"]["fixed"] = @benchmarkable solve(
    $prob, IRKGL16(); adaptive = false, dt = 0.05
)
