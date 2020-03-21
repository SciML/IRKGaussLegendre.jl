__precompile__()

module IRKGaussLegendre

using Reexport
@reexport using DiffEqBase

using LinearAlgebra,StaticArrays
using Parameters
using OrdinaryDiffEq
using Printf

include("IRKGL16Solver.jl")

export IRKGL16,IRKAlgorithm
export tcoeffs

end # module
