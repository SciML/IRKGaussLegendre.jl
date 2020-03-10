__precompile__()

module IRK

using Reexport
@reexport using DiffEqBase

using LinearAlgebra
using Parameters
using OrdinaryDiffEq

include("IRK8.jl")

export IRK8,IRKAlgorithm
export IRK8Coefficients

end # module
