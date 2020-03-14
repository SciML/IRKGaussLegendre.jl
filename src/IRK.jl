__precompile__()

module IRK

using Reexport
@reexport using DiffEqBase

using LinearAlgebra,StaticArrays
using Parameters
using OrdinaryDiffEq

include("IRK8.jl")

export IRK8,IRKAlgorithm
export tcoeffs

end # module
