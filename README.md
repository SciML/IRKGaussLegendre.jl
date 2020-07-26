# IRKGL16
Implicit Runge-Kutta Gauss-Legendre 16th order
(implementation in Julia)

## Installation

This package can be installed using

```julia
using Pkg
Pkg.add("https://github.com/mikelehu/IRKGaussLegendre.jl")
using IRKGaussLegendre
```

## Example: Burrau's problem of three bodies

Three point masses attract each other according to the Newtonian law of gravitation. The masses of the particles are
m1=3, m2=4, and m3=5; they are initially located at the apexes of a right triangle with sides 3, 4, and 5, so that the
corresponding masses and sides are opposite. The particles are free to move in the plane of the triangle and are at rest initially.


## Step 1: Defining  the problem

To solve this numerically, we define a problem type by giving it the equation, the initial
condition, and the timespan to solve over:


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

## Step 2: Solving the problem


After defining a problem, you solve it using solve

```julia
sol1=solve(prob,IRKGL16(),reltol=1e-12, abstol=1e-12);
```

## Step 3: Analyzing the solution

```julia
bodylist = ["Body-1", "Body-2", "Body-3"]
pl = plot(title="Burrau problem (Adaptive)",aspect_ratio=1)

ulist1 = sol1.u[1:end]
tlist1 = sol1.t[1:end]

for j = 1:3
 xlist  = map(u->u[2,1,j], ulist1)
 ylist  = map(u->u[2,2,j], ulist1)
 pl2 = plot!(xlist,ylist, label = bodylist[j])   
end  
plot(pl2)
```
![Burrau problem](/ODEProblems/Burrau.png)

## Implementation details

[Antoñana, M., Makazaga, J., Murua, Ander. "Reducing and monitoring round-off error propagation
for symplectic implicit Runge-Kutta schemes."  Numerical Algorithms. 2017.](https://doi.org/10.1007/s11075-017-0287-z)
