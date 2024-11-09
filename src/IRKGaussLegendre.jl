__precompile__()

module IRKGaussLegendre

using Reexport
@reexport using SciMLBase
using DiffEqBase
import FastBroadcast

using LinearAlgebra
using Parameters
using RecursiveArrayTools
using SIMD

const CompiledFloats = Union{Float32, Float64}

include("IRKCoefficients.jl")
include("./simd/VecArray_def.jl")
include("IRKGL16Solver.jl")
include("IRKGL16AuxFunctions.jl")
include("IRKGL16step_fixed_seq.jl")
include("IRKGL16step_adaptive_seq.jl")
include("./par/IRKGL16step_fixed_par.jl")
include("./par/IRKGL16step_adaptive_par.jl")
include("./simd/IRKGL16step_fixed_SIMD.jl")
include("./simd/IRKGL16step_adaptive_SIMD.jl")

SciMLBase.allows_arbitrary_number_types(::IRKAlgorithm) = true
SciMLBase.isautodifferentiable(::IRKAlgorithm) = true
SciMLBase.forwarddiffs_model(::IRKAlgorithm) = false
SciMLBase.forwarddiffs_model_time(::IRKAlgorithm) = false
SciMLBase.allowscomplex(::IRKAlgorithm) = true
SciMLBase.isadaptive(::IRKAlgorithm) = true
SciMLBase.isdiscrete(::IRKAlgorithm) = false

export IRKGL16, IRKAlgorithm
export tcoeffs, CompiledFloats

end # module
