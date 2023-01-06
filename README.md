# IRKGaussLegendre.jl

[![Build Status](https://github.com/SciML/IRKGaussLegendre.jl/workflows/CI/badge.svg)](https://github.com/SciML/IRKGaussLegendre.jl/actions?query=workflow%3ACI)

IRKGaussLegendre.jl is an efficient Julia implementation of an implicit Runge-Kutta Gauss-Legendre 16th order method.
The method is fully integrated into the **DifferentialEquations.jl ecosystem** for high-performance high-precision
integration.

Required Julia 1.5 version or higher

## Description

We present a Julia implementation of a 16th order Implicit Runge-Kutta integrator IRKGL16 (a 8-stage
IRK scheme based on Gauss-Legendre nodes) for **high accuracy** numerical integration of non-stiff
ODE systems. Our algorithm supports **adaptive timesteping, mixed precision and multithreading** to
solve problems fast and accuracy

The family of implicit Runge-Kutta schemes based on collocation with Gauss-Legendre nodes are
known to be symplectic and super-convergent (order 2s for the method with s internal nodes), and
thus very convenient for the high precision numerical integration of Hamiltonian systems with
constant time-step size. For **non-stiff problems**, implementations based on fixed-point iterations are
recommended

We believe that, for general (non-necessarily Hamiltonian) non-stiff ODE systems, such implicit
Runge-Kutta methods (implemented with fixed point iteration) can be very competitive for high
precision computations (for accuracy requirements that exceeds double precision arithmetic)


## Installation

This package can be installed using

```julia
julia>using Pkg
julia>Pkg.add("IRKGaussLegendre.jl")
julia>using IRKGaussLegendre
```

## Solver options

### Available common arguments

- dt: stepsize
- save_everystep: default is true
- adaptive: =true (adaptive timestepping); =false (fixed timestepping)
- maxiters: maximum number of iterations before stopping


### No-common arguments


- initial_interp: initialization method for stages.
        - =false  simplest initialization
        - =true (default) interpolating from the stage values of previous step
- mstep: output saved at every 'mstep' steps. Default 1.
- myoutputs: default false
- maxtrials: maximum number of attempts to accept adaptive step size
- threading
      - =false (default): sequential execution of the numerical integration
      - =true: parallel execution (stage-wise parallelization)
- mixed_precision
      - =false (default)
      - =true: combine 'base precision arithmetic' with precision specified
      in low_prec_type variable
- low_prec_type: (Float64, Float32,...)
- nrmbits: number of bits to remove when applying the stop criterion

## Return Codes

The solution types have a retcode field which returns a symbol signifying the error state of the solution. The retcodes are as follows:

- ReturnCode.Success: The integration completed without erroring.
- ReturnCode.Failure: General uncategorized failures or errors.

## Example: Burrau's problem of three bodies

Three point masses attract each other according to the Newtonian law of gravitation. The masses of the particles are
m1=3, m2=4, and m3=5; they are initially located at the apexes of a right triangle with sides 3, 4, and 5, so that the
corresponding masses and sides are opposite. The particles are free to move in the plane of the triangle and are at rest initially.

Szebehely, V. 1967, "Burrau's Problem of Three Bodies", Proceedings of the National Academy of Sciences of the United States of America, vol. 58, Issue 1, pp. 60-65 [postscript file](http://www.ucolick.org/~laugh/oxide/projects/szebehely1.ps)


### Step 1: Defining  the problem

To solve this numerically, we define a problem type by giving it the equation, the initial
condition, and the timespan to solve over:

```julia
using IRKGaussLegendre
using Plots, LinearAlgebra, LaTeXStrings
```

```julia
function NbodyODE!(du,u,Gm,t)
     N = length(Gm)
     du[1,:,:] .= 0
     for i in 1:N
        qi = u[2,:,i]
        Gmi = Gm[i]
        du[2,:,i] = u[1,:,i]
        for j in (i+1):N
           qj = u[2,:,j]
           Gmj = Gm[j]
           qij = qi - qj
           auxij = (qij[1]*qij[1]+qij[2]*qij[2]+qij[3]*qij[3])^(-3/2)
           du[1,:,i] -= Gmj*auxij*qij
           du[1,:,j] += Gmi*auxij*qij
        end
     end

    return
end
```

```julia
Gm = [5, 4, 3]
N=length(Gm)
q=[1,-1,0,-2,-1,0,1,3,0]
v=zeros(size(q))
q0 = reshape(q,3,:)
v0 = reshape(v,3,:)
u0 = Array{Float64}(undef,2,3,N)
u0[1,:,:] = v0
u0[2,:,:] = q0
tspan = (0.0,63.0)
prob=ODEProblem(NbodyODE!,u0,tspan,Gm);
```

### Step 2: Solving the problem


After defining a problem, you solve it using solve

```julia
sol1=solve(prob,IRKGL16(),adaptive=true, reltol=1e-12, abstol=1e-12);
```

### Step 3: Analyzing the solution


#### Orbits


```julia
bodylist = ["Body-1", "Body-2", "Body-3"]
pl = plot(title="Burrau problem (Adaptive)",aspect_ratio=1)

ulist1 = sol1.u[1:end]
tlist1 = sol1.t[1:end]

for j = 1:3
 xlist  = map(u->u[2,1,j], ulist1)
 ylist  = map(u->u[2,2,j], ulist1)
 pl = plot!(xlist,ylist, label = bodylist[j])   
end  
plot(pl)
```
![Burrau problem](/Examples/BurrauOrbits.png)


#### Step Size

```julia
plot(xlabel="t", ylabel="step size",title="Adaptive step size")
steps1 =sol1.t[2:end]-sol1.t[1:end-1]
plot!(sol1.t[2:end],steps1)
```

![Burrau problem](/Examples/BurrauStepSize.png)


#### Energy-Error

```julia
function NbodyEnergy(u,Gm)
     N = length(Gm)
     zerouel = zero(eltype(u))
     T = zerouel
     U = zerouel
     for i in 1:N
        qi = u[2,:,i]
        vi = u[1,:,i]
        Gmi = Gm[i]
        T += Gmi*(vi[1]*vi[1]+vi[2]*vi[2]+vi[3]*vi[3])
        for j in (i+1):N
           qj = u[2,:,j]  
           Gmj = Gm[j]
           qij = qi - qj
           U -= Gmi*Gmj/norm(qij)
        end
     end
    1/2*T + U
end
```

```julia
setprecision(BigFloat, 256)
u0Big=BigFloat.(u0)
GmBig=BigFloat.(Gm)

E0=NbodyEnergy(u0Big,GmBig)
ΔE1 = map(x->NbodyEnergy(BigFloat.(x),GmBig), sol1.u)./E0.-1
plot(title="Energy error", xlabel="t", ylabel=L"\Delta E")
plot!(sol1.t,log10.(abs.(ΔE1)), label="")
```

![Burrau problem](/Examples/BurrauEnergyError.png)


## More Examples

[Benchmark examples](https://github.com/SciML/IRKGaussLegendre.jl/tree/master/Benchmarks)

## Implementation details

[Antoñana, M., Makazaga, J., Murua, Ander. "Reducing and monitoring round-off error propagation
for symplectic implicit Runge-Kutta schemes."  Numerical Algorithms. 2017.](https://doi.org/10.1007/s11075-017-0287-z)
