using SciMLTesting, IRKGaussLegendre, JET, Test

# The SciML common interface IRKGaussLegendre deliberately reexports so that
# `using IRKGaussLegendre` is enough to build and solve an ODE problem, as the README's
# Quick Start does. Owned and documented upstream; kept in sync with the reexport
# `export` block in src/IRKGaussLegendre.jl.
const REEXPORTS = (
    :DEStats, :EnsembleAnalysis, :EnsembleDistributed, :EnsembleProblem, :EnsembleSerial,
    :EnsembleSolution, :EnsembleSplitThreads, :EnsembleSummary, :EnsembleThreads,
    :NullParameters, :ODEFunction, :ODEProblem, :ODESolution, :ReturnCode, :remake,
    :solve, :successful_retcode,
)

run_qa(IRKGaussLegendre; reexports_allow = REEXPORTS)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using IRKGaussLegendre`,
    # so the allow-list cannot drift into approving names the package no longer provides.
    @testset "$name" for name in REEXPORTS
        @test name in names(IRKGaussLegendre)
        @test isdefined(@__MODULE__, name)
    end
end
