__precompile__()

module IRKGaussLegendre

using Reexport
@reexport using DiffEqBase

using LinearAlgebra,StaticArrays
using Parameters
using OrdinaryDiffEq
using Printf
using DoubleFloats

const CompiledFloats = Union{Float32,Float64}

include("IRKGL16Solver.jl")
include("IRKGL16Solver2.jl")

export IRKGL16,IRKAlgorithm
export IRKGL162,IRKAlgorithm2
export tcoeffs, CompiledFloats

end # module
