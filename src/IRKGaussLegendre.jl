__precompile__()

module IRKGaussLegendre

using Reexport
@reexport using DiffEqBase

using LinearAlgebra
using Parameters
using OrdinaryDiffEq
using RecursiveArrayTools

const CompiledFloats = Union{Float32, Float64}

include("IRKCoefficients.jl")
include("IRKGL16AuxFunctions.jl")
include("IRKGL16Solver.jl")

export IRKGL16, IRKAlgorithm
export tcoeffs, CompiledFloats

end # module
