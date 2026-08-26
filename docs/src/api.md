# Public API

IRKGaussLegendre exports the following solver configuration and low-level
coefficient types. Problem construction and `solve` are provided by
[SciMLBase](https://docs.sciml.ai/SciMLBase/stable/).

```@docs
IRKGL16
IRKAlgorithm
CompiledFloats
tcoeffs
```

## Reexported SciML common interface

`using IRKGaussLegendre` also brings in the parts of the SciML common interface needed to
build an ODE problem, solve it with `IRKGL16`, and inspect the result, so they do not
have to be imported separately — this is what the Quick Start example in the README
relies on. These names are owned and documented by
[SciMLBase](https://docs.sciml.ai/SciMLBase/stable/):

  - Problems: `ODEProblem`, `EnsembleProblem`
  - Functions: `ODEFunction`
  - Solutions: `ODESolution`, `EnsembleSolution`, `EnsembleSummary`, `DEStats`
  - Ensemble algorithms: `EnsembleSerial`, `EnsembleThreads`, `EnsembleDistributed`,
    `EnsembleSplitThreads`, and the `EnsembleAnalysis` module
  - Solving: `solve`, `remake`
  - Return status: `ReturnCode`, `successful_retcode`
  - `NullParameters`

Callbacks and the iterator interface are deliberately *not* reexported: `IRKGL16`
implements only `solve` and honors no `callback` keyword, so those names would resolve to
methods that cannot be used with this solver. Second-order systems are handled by passing
`second_order_ode = true` to `IRKGL16` alongside an ordinary `ODEProblem`, so
`SecondOrderODEProblem` is not part of this surface either. Anything else from SciMLBase
must be imported from SciMLBase directly.
