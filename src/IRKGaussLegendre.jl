__precompile__()

module IRKGaussLegendre

    import SciMLBase
    using DiffEqBase: DEVerbosity

    # The SciML common interface that IRKGaussLegendre reexports (see the second
    # `export` below), so that `using IRKGaussLegendre` on its own is enough to build an
    # ODE problem, solve it with `IRKGL16`, and inspect the result -- which is exactly
    # what the README's Quick Start and the docs examples do. Callbacks and the iterator
    # interface are deliberately absent: this solver implements only `__solve` and
    # honors no `callback` keyword. Every name stays owned and documented upstream.
    using SciMLBase: DEStats, EnsembleAnalysis, EnsembleDistributed, EnsembleProblem,
        EnsembleSerial, EnsembleSolution, EnsembleSplitThreads, EnsembleSummary,
        EnsembleThreads, NullParameters, ODEFunction, ODEProblem, ODESolution,
        ReturnCode, remake, solve, successful_retcode
    using SciMLLogging: AbstractVerbosityPreset, Standard, @SciMLMessage
    import LinearAlgebra
    using Parameters: @unpack
    using PrecompileTools: @compile_workload, @setup_workload
    using SIMD: Vec, vload

    """
            CompiledFloats

        Floating-point element types supported by the compiled SIMD IRK implementation.

        ```jldoctest
        julia> using IRKGaussLegendre

        julia> Float64 <: CompiledFloats
        true
        ```
    """
    const CompiledFloats = Union{Float32, Float64}

    include("IRKCoefficients.jl")
    include("./simd/VecArray_def.jl")
    include("compensated_sum.jl")
    include("IRKGL16Solver.jl")
    include("IRKGL16AuxFunctions.jl")
    include("IRKGL16step_fixed_seq.jl")
    include("IRKGL16step_adaptive_seq.jl")
    include("./simd/IRKGL16step_fixed_SIMD.jl")
    include("./simd/IRKGL16step_adaptive_SIMD.jl")
    include("./hybrid/IRKGL16step_fixed_HYBR.jl")
    include("./hybrid/IRKGL16step_adaptive_HYBR.jl")

    SciMLBase.allows_arbitrary_number_types(::IRKAlgorithm) = true
    SciMLBase.isautodifferentiable(::IRKAlgorithm) = true
    SciMLBase.forwarddiffs_model(::IRKAlgorithm) = false
    SciMLBase.forwarddiffs_model_time(::IRKAlgorithm) = false
    SciMLBase.allowscomplex(::IRKAlgorithm) = true
    SciMLBase.isadaptive(::IRKAlgorithm) = true
    SciMLBase.isdiscrete(::IRKAlgorithm) = false

    export IRKGL16, IRKAlgorithm
    export tcoeffs, CompiledFloats

    # Reexported SciML common interface; approved via `reexports_allow` in test/qa/qa.jl.
    export DEStats, EnsembleAnalysis, EnsembleDistributed, EnsembleProblem,
        EnsembleSerial, EnsembleSolution, EnsembleSplitThreads, EnsembleSummary,
        EnsembleThreads, NullParameters, ODEFunction, ODEProblem, ODESolution,
        ReturnCode, remake, solve, successful_retcode

    @setup_workload begin
        @compile_workload begin
            function precompile_rhs!(du, u, p, t)
                du[1] = u[1]
                return nothing
            end

            prob = SciMLBase.ODEProblem(precompile_rhs!, [1.0], (0.0, 0.1))
            SciMLBase.solve(prob, IRKGL16(); adaptive = false, dt = 0.05)
        end
    end

end # module
