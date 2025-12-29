# IRKGaussLegendre.jl

IRKGaussLegendre.jl is an efficient Julia implementation of a 16th order implicit Runge-Kutta Gauss-Legendre method for high-precision numerical integration of non-stiff ODE systems.

## Features

- **High-order accuracy**: 16th order convergence (8-stage IRK scheme based on Gauss-Legendre nodes)
- **Symplectic**: Preserves geometric structure for Hamiltonian systems
- **Adaptive time-stepping**: Automatic step size control for efficiency and accuracy
- **SIMD-vectorization**: Optimized implementation for Float32/Float64 computations
- **Multi-threading support**: Parallel execution for stage-wise computations
- **DifferentialEquations.jl integration**: Works seamlessly with the SciML ecosystem

## Installation

```julia
using Pkg
Pkg.add("IRKGaussLegendre")
```

## Quick Start

```julia
using IRKGaussLegendre, OrdinaryDiffEq

# Define your ODE
function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - β * u[3]
end

# Set up and solve the problem
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
p = (10.0, 28.0, 8/3)
prob = ODEProblem(lorenz!, u0, tspan, p)

sol = solve(prob, IRKGL16(), reltol=1e-12, abstol=1e-12)
```

## API Reference

```@docs
IRKGL16
```

## Solver Options

### Common Keyword Arguments

- `dt`: Step size
- `saveat`: Specific times to save the solution. If given a number, expands to `tspan[1]:saveat:tspan[2]`
- `save_everystep`: Save result at every step (default: `true`)
- `adaptive`: Enable adaptive time-stepping (default: `true`)
- `maxiters`: Maximum number of fixed-point iterations
- `abstol`: Absolute tolerance for adaptive time-stepping
- `reltol`: Relative tolerance for adaptive time-stepping

### IRKGL16-specific Options

- `second_order_ode`: Set to `true` for second-order differential equations (default: `false`)
- `simd`: Enable SIMD-vectorized implementation for Float32/Float64 (default: `false`)
- `initial_extrapolation`: Use extrapolation from previous step for stage initialization (default: `true`)
- `maxtrials`: Maximum attempts to accept adaptive step size (default: `5`)
- `threading`: Enable multi-threading for stage-wise parallelization (default: `false`)

## Return Codes

- `ReturnCode.Success`: Integration completed successfully
- `ReturnCode.Failure`: Integration failed

## See Also

- [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) - The SciML differential equations ecosystem
- [README](https://github.com/SciML/IRKGaussLegendre.jl) - Full examples including the Pythagorean three-body problem
