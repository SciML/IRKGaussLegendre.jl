# IRKGaussLegendre.jl

[![Build Status](https://github.com/SciML/IRKGaussLegendre.jl/workflows/CI/badge.svg)](https://github.com/SciML/IRKGaussLegendre.jl/actions?query=workflow%3ACI)


IRKGaussLegendre.jl provides an efficient Julia implementation of a 16th-order implicit Runge–Kutta method based on 8-stage Gauss–Legendre collocation.



Requires Julia 1.5 version or higher

## Installation

This package can be installed using

```julia
julia>using Pkg
julia>Pkg.add("IRKGaussLegendre")
julia>using IRKGaussLegendre
```

## Quick Start example: Hamiltonian System

This example solves the Hamiltonian system with Hamiltonian function,

$$H(q,p) = \tfrac12(p_1^2 + p_2^2 + q_1^2 + q_2^2).$$

This is a harmonic oscillator, whose equations of motion are
$$\dot{q} = p, \qquad \dot{p} = -q.$$


```julia
using IRKGaussLegendre
function harmonic!(du, u, p, t)
    q1, q2, p1, p2 = u
    du[1] = p1
    du[2] = p2
    du[3] = -q1
    du[4] = -q2
end

u0 = [1.0, 0.0, 0.0, 1.0]
tspan = (0.0, 100.0)

prob = ODEProblem(harmonic!, u0, tspan)
sol = solve(prob, IRKGL16(second_order_ode=true); dt = 1.0, adaptive = false)
```

## Overview of IRKGaussLegendre.jl

IRKGaussLegendre.jl provides high-order implicit Runge–Kutta solvers based on Gauss–Legendre collocation, implemented within the `DifferentialEquations.jl` framework. The package is intended for the accurate time integration of non-stiff ordinary differential equations, with particular emphasis on
Hamiltonian and other conservative systems.

The implemented Gauss–Legendre methods are symplectic, A-stable, and superconvergent, achieving order \(2s\) for an \(s\)-stage scheme. These properties make them well suited for long-time simulations where preservation of geometric structure and invariants is important.

IRKGaussLegendre.jl supports systems written in ODEs  formulations. Several implementation variants are provided, including fully SIMD-vectorized, Hybrid, and sequential execution paths, following the SIMD-based approach described in [2].Also a variant for numerically integrating system of second-order ODE structure more efficiently is available.

The solvers integrate seamlessly into the SciML ecosystem and therefore inherit standard features such as problem definitions via ODEProblem, flexible output control, and optional adaptive time stepping.

### Available IRKGL16 variants:

These presets determine how much of the computation is SIMD-vectorized.

- **Fully vectorized implementation**: (`simd=true, fseq=false`)
- **Hybrid implementation** : It is intended for cases where the user-supplied function is not compatible with `SIMD.jl`. All computations except the user-supplied function evaluations are SIMD-vectorized (`simd=true, fseq=true`)
- **Fully sequential implementation**: It is used for state variables of type other than Float32 or Float64 (e.g., BigFloat) or whenever the **simd** option is explicitly disabled (`simd=false`) by the user

### Mathematical problem setting

- **General first-order problem**

We consider initial value problems for systems of ordinary differential equations (ODEs) of the form:

$$\frac{d}{dt}u=f(t,u),\quad  u(t_0)=u_0 \tag{1}$$

where $f: \mathbb{R}^{D+1} \to \mathbb{R}^D$ is a sufficiently smooth function and $u_0 \in \mathbb{R}^D$.

- **Second-order systems**

Second-order differential equations of the form

$$\frac{d^2q}{dt^2}=g(t,q)$$

are particularly common in Hamiltonian and mechanical systems where $q \in \mathbb{R}^d$. Such systems can be numerically integrated 
by rewriting them as a first-order system of the form (1) with dimension $D=2d$, state vector $u=(q,v)$, where $v=\frac{d}{dt} q$, and

$$f(t,u)=(v,g(t,q)).$$

The user may indicate that the ODE problem has this special structure by setting the option `second_order_ode = true`, as shown in the harmonic oscillator Quick Start example above. If this option is not explicitly enabled (it is disabled by default), the solver treats the problem as a general first-order ODE system, which typically results in higher CPU time.



### Solver options

#### Available common keyword arguments

- **`dt`**: stepsize
- **`saveat`**: denotes specific times to save the solution at, during the solving phase. If saveat is given a number, then it will automatically expand to tspan[1]:saveat:tspan[2]. The default value is []
- **`save_everystep`**: saves the result at every step. Default is true
- **`adaptive`**: =true (adaptive timestepping); =false (constant timestepping)
- **`maxiters`**: maximum number of fixed-point iterations before stopping
- **`abstol`**: absolute tolerance in adaptive timestepping (default abstol=1e-6)
- **`reltol`**: relative tolerance in adaptive timestepping (default reltol=1e-6)


#### Solver-specific keyword arguments

- **`s`**: the number of stages of the method, which determines its order 2s. Valid choices are s=8 (default value), s=6, s=4 and s=2.

* **`second_order_ode`** (boolean):

  * `true`: Indicates that the `ODEProblem` has a second-order ODE structure.
  * `false` (default): Indicates that the problem should be treated as a general first-order system.

- **`simd`** (boolean):
  - `true`: SIMD-vectorized implementation (only  for Float32 or Float64)
  - `false` (default):  use the sequential implementation (works with arbitrary number types) 

- **`fseq`** (boolean):
  - `true` (default): user-supplied function evaluations are not SIMD-vectorized
  - `false`: user-supplied function evaluations are SIMD-vectorized

- **`initial_extrapolation`**: initialization method for stages.
  - `false`: basic initialization method
  - `true` (default): extrapolation from the stage values of the previous step


- **`maxtrials`**: maximum number of attempts to accept adaptive step size





### Return Codes

The solver uses the standard ReturnCode symbols from SciMLBase for compatibility.


## Extended Example: Pythagorean Three-Body Problem

A full demonstration of IRKGL16 on the Pythagorean three-body problem  
(including performance comparisons and energy error plots) is available here:

**[Pythagorean three-body example](docs/src/examples/pythagorean_three_body.md)**


## More Examples

[Benchmark examples](https://github.com/mikelehu/Implicit_symplectic_can_outperform_explicit_symplectic/tree/main/Splitting%20Solvers%20Implementation/Other%20Benchmarks)

## Implementation details

[1]  [Antoñana, M., Makazaga, J., and Murua, A. "Reducing and monitoring round-off error propagation
for symplectic implicit Runge-Kutta schemes."  Numerical Algorithms. 2017.](https://doi.org/10.1007/s11075-017-0287-z)

[2] [Antoñana, M., Makazaga,J., and Murua, A. "SIMD-vectorized implicit symplectic integrators can outperform explicit symplectic ones."  2025.](https://doi.org/10.48550/arXiv.2511.03655)



## Contact

If you have any questions or suggestions, feel free to open an issue or contact us at mikel.antonana@ehu.eus.

Updated December 17th, 2025
