# IRKGaussLegendre.jl

[![Build Status](https://github.com/SciML/IRKGaussLegendre.jl/workflows/CI/badge.svg)](https://github.com/SciML/IRKGaussLegendre.jl/actions?query=workflow%3ACI)

IRKGaussLegendre.jl is an efficient Julia implementation of an implicit Runge-Kutta Gauss-Legendre 16th order method.
The method is fully integrated into the **DifferentialEquations.jl ecosystem** for  high-precision
integration.

Required Julia 1.5 version or higher

## Description

We present a Julia implementation of a 16th order Implicit Runge-Kutta integrator IRKGL16 (a 8-stage
IRK scheme based on Gauss-Legendre nodes) for **high accuracy** numerical integration of non-stiff
ODE systems. Our algorithm supports **adaptive time-steping and SIMD-vectorization** to
solve problems fast and accurately

The family of implicit Runge-Kutta schemes based on collocation with Gauss-Legendre nodes is
known to be symplectic and super-convergent (order 2s for the method with s internal nodes), making them very convenient for  high-precision numerical integration of Hamiltonian systems with
constant time-step size. For **non-stiff problems**, implementations based on fixed-point iterations are
recommended

We believe that, for general (non-necessarily Hamiltonian) non-stiff ODE systems, such implicit
Runge-Kutta methods (implemented with fixed point iteration) can be very competitive for high
precision computations. We show that a vectorized implementation of IRKGL16 that exploits the SIMD-based parallelism offered by modern processor can be more efficient than high order explicit Runge-Kutta methods even for double precision computations (**see the following experiment as an example**).


## Installation

This package can be installed using

```julia
julia>using Pkg
julia>Pkg.add("IRKGaussLegendre.jl")
julia>using IRKGaussLegendre
```

## Solver options

### Available common keyword arguments

- dt: stepsize
- saveat: denotes specific times to save the solution at, during the solving phase. If saveat is given a number, then it will automatically expand to tspan[1]:saveat:tspan[2]. For methods where interpolation is not possible, saveat may be equivalent to tstops. The default value is [].
- save_everystep: saves the result at every step. Default is true (see keyword argument mstep)
- adaptive: =true (adaptive timestepping); =false (fixed timestepping)
- maxiters: maximum number of fixed-point iterations before stopping
- abstol: absolute tolerance in adaptive timestepping
- reltol: relative tolerance in adaptive timestepping


### No-common keyword arguments


- second_order_ode (boolean):

      - =false (default): for a ODEProblem type 
      - =true: for a second-order differential equation 

- simd (boolean):

      - =true: SIMD-vectorized implementation only available for Float32 or Float64 computations
      - =false (default):  generic implementation that can use with arbitrary Julia-defined number systems 




- initial_extrapolation: initialization method for stages.

        - =false: simplest initialization
        - =true (default): extrapolation from the stage values of previous step


- maxtrials: maximum number of attempts to accept adaptive step size

- threading

      - =false (default): sequential execution of the numerical integration
      - =true: computations using threads (shared memory multi-threading) for stage-wise parallelization
      


## Return Codes

The solution types have a retcode field which returns a symbol signifying the error state of the solution. The retcodes are as follows:

- ReturnCode.Success: The integration completed without erroring.
- ReturnCode.Failure: General uncategorized failures or errors.

## Example: Pythagorean three body problem

Three point masses attract each other according to the Newtonian law of gravitation. The masses of the particles are
m1=3, m2=4, and m3=5; they are initially located at the apexes of a right triangle with sides 3, 4, and 5, so that the
corresponding masses and sides are opposite. The particles are free to move in the plane of the triangle and are at rest initially.

Szebehely, V. 1967, "Burrau's Problem of Three Bodies", Proceedings of the National Academy of Sciences of the United States of America, vol. 58, Issue 1, pp. 60-65 [postscript file](http://www.ucolick.org/~laugh/oxide/projects/szebehely1.ps)


### Step 1: Defining  the problem

To solve this numerically, we define a problem type by giving it the equation, the initial
condition, and the timespan to solve over:

```julia
using IRKGaussLegendre
using OrdinaryDiffEq
using Plots, LinearAlgebra, LaTeXStrings
using BenchmarkTools
```

```julia
function NbodyODE!(F,u,Gm,t)
     N = length(Gm)
     for i in 1:N
        for k in 1:3
            F[k, i, 2] = 0
        end
     end
     for i in 1:N
        xi = u[1,i,1]
        yi = u[2,i,1]
        zi = u[3,i,1]
        Gmi = Gm[i]
        for j in i+1:N
            xij = xi - u[1,j,1]
            yij = yi - u[2,j,1]
            zij = zi - u[3,j,1]
            Gmj = Gm[j]
            dotij = (xij*xij+yij*yij+zij*zij)
            auxij = 1/(sqrt(dotij)*dotij)
            Gmjauxij = Gmj*auxij
            F[1,i,2] -= Gmjauxij*xij
            F[2,i,2] -= Gmjauxij*yij
            F[3,i,2] -= Gmjauxij*zij
            Gmiauxij = Gmi*auxij
            F[1,j,2] += Gmiauxij*xij
            F[2,j,2] += Gmiauxij*yij
            F[3,j,2] += Gmiauxij*zij
        end
     end
     for i in 1:3, j in 1:N
        F[i,j,1] = u[i,j,2]
     end
    return nothing
end
```

```julia
Gm = [5, 4, 3]
N=length(Gm)
q=[1,-1,0,-2,-1,0,1,3,0]
v=zeros(size(q))
q0 = reshape(q,3,:)
v0 = reshape(v,3,:)
u0 = Array{Float64}(undef,3,N,2)
u0[:,:,1] = q0
u0[:,:,2] = v0
tspan = (0.0,63.0)
prob=ODEProblem(NbodyODE!,u0,tspan,Gm);
```

### Step 2: Solving the problem


After defining a problem, you solve it using solve:

#### Generic implementation

```julia
sol1=solve(prob,IRKGL16(second_order_ode=true, simd=false),adaptive=true, reltol=1e-14, abstol=1e-14);
```

```julia
@btime solve(prob,IRKGL16(second_order_ode=true, simd=false),adaptive=true, reltol=1e-14, abstol=1e-14);
```
24.989 ms (6018 allocations: 1.98 MiB)


#### SIMD-vectorized implementation 
- **6 x faster than generic implementation !!**
- **Faster than DPRKN12 !!**

```julia
sol2=solve(prob,IRKGL16(second_order_ode=true, simd=true),adaptive=true, reltol=1e-14, abstol=1e-14)
```

```julia
@btime solve(prob,IRKGL16(second_order_ode=true, simd=true),adaptive=true, reltol=1e-14, abstol=1e-14);
```
4.386 ms (6000 allocations: 2.31 MiB)


#### DPRKN12

```julia
function NbodyODE2nd!(ddu,du,u,Gm,t)

     N = length(Gm)

     for i in 1:N
         for k in 1:3
             ddu[k,i]= 0
         end
     end

     for i in 1:N
        xi = u[1,i]
        yi = u[2,i]
        zi = u[3,i]
        Gmi = Gm[i]
        for j in (i+1):N
           xij = xi - u[1,j]
           yij = yi - u[2,j]
           zij = zi - u[3,j]
           Gmj = Gm[j]
           dotij = (xij*xij+yij*yij+zij*zij)
           auxij = 1/(sqrt(dotij)*dotij)
           Gmjauxij = Gmj*auxij
           ddu[1,i] -= Gmjauxij*xij
           ddu[2,i] -= Gmjauxij*yij
           ddu[3,i] -= Gmjauxij*zij
           Gmiauxij = Gmi*auxij
           ddu[1,j] += Gmiauxij*xij
           ddu[2,j] += Gmiauxij*yij
           ddu[3,j] += Gmiauxij*zij
        end
     end

    return nothing

end
```

```julia
q0=u0[:,:,1]
v0=u0[:,:,2]
prob2nd = SecondOrderODEProblem(NbodyODE2nd!,v0,q0,tspan,Gm)
sol3 =solve(prob2nd,DPRKN12(),abstol=1e-14,reltol=1e-14);
```

```julia
@btime solve(prob2nd,DPRKN12(),abstol=1e-14,reltol=1e-14, save_everystep=false);
```
 6.082 ms (58656 allocations: 1.05 MiB)


**Convert solution data** to u0 format
```julia
nk=length(sol3.t)
etype=eltype(u0)
sol3u=Vector{Array{etype, 3}}(undef,nk)
sol3t=Vector{etype}(undef,nk)
uk=copy(u0)

for k in 1:nk
    uk[:,:,1]=sol3.u[k].x[2]
    uk[:,:,2]=sol3.u[k].x[1]
    sol3u[k]=copy(uk)
    sol3t[k]=sol3.t[k]
end
```




### Step 3: Analyzing the solution




#### Orbits


```julia
bodylist = ["Body-1", "Body-2", "Body-3"]
pl = plot(title="Pythagorean problem",xlabel="x", ylabel="y",aspect_ratio=1)

ulist1 = sol1.u[1:end]
tlist1 = sol1.t[1:end]

for j = 1:3
 xlist  = map(u->u[1,j,1], ulist1)
 ylist  = map(u->u[2,j,1], ulist1)
 pl = plot!(xlist,ylist, label = bodylist[j])   
end  
plot(pl)
```
![Burrau problem](/Examples/PythagoreanOrbits.png)


#### Step Size

```plot(xlabel="t", ylabel="step size",title="Adaptive step size")
steps1 =sol1.t[2:end]-sol1.t[1:end-1]
plot!(sol1.t[2:end],steps1, label="generic")
steps2 =sol2.t[2:end]-sol2.t[1:end-1]
plot!(sol2.t[2:end],steps2, label="simd")
```

![Burrau problem](/Examples/PythagoreanStepsize.png)


#### Energy-Error

```julia
function NbodyEnergy(u,Gm)
     N = length(Gm)
     zerouel = zero(eltype(u))
     T = zerouel
     U = zerouel
     for i in 1:N
        qi = u[:,i,1]
        vi = u[:,i,2]
        Gmi = Gm[i]
        T += Gmi*(vi[1]*vi[1]+vi[2]*vi[2]+vi[3]*vi[3])
        for j in (i+1):N
           qj = u[:,j,1]
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
ΔE2 = map(x->NbodyEnergy(BigFloat.(x),GmBig), sol2.u)./E0.-1
plot(title="Error in energy", legend=:bottomright,
     xlabel="t", ylabel=L"log10(\Delta E)")
plot!(sol1.t[2:end], abs.(ΔE1[2:end]), yscale=:log10, label="IRKGL16-generic")
plot!(sol2.t[2:end], abs.(ΔE2[2:end]), yscale=:log10, label="IRKGL16-simd")
plot!(sol3t[2:end], abs.(ΔE3[2:end]), yscale=:log10, label="DPRKN12")
```

![Burrau problem](/Examples/PythagoreanEnergyError.png)


## More Examples

[Benchmark examples](https://github.com/mikelehu/Implicit_symplectic_can_outperform_explicit_symplectic/tree/main/Other%20Benchmarks)

## Implementation details

- [Antoñana, M., Makazaga, J., Murua, Ander. "Reducing and monitoring round-off error propagation
for symplectic implicit Runge-Kutta schemes."  Numerical Algorithms. 2017.](https://doi.org/10.1007/s11075-017-0287-z)

- [Antoñana, M., Murua, Ander. "SIMD-vectorized implicit symplectic integrators
can outperform explicit symplectic ones."  2025.](https://github.com/mikelehu/Implicit_symplectic_can_outperform_explicit_symplectic)



## Contact

If you have any questions or suggestions, feel free to open an issue or contact us at mikel.antonana@ehu.eus.

Updated July 29th, 2025
