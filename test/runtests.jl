using IRKGaussLegendre, ODEProblemLibrary, DiffEqDevTools, Test
import ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear

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
sol1 = solve(prob, IRKGL16(), reltol = 1e-12, abstol = 1e-12)

# Analytical tests

sol = solve(prob_ode_2Dlinear, IRKGL16())
@test sol.errors[:l2] < 1e-16

# dts = (1 // 2) .^ (4:-1:2)
dts = BigFloat.( (1 // 2) .^ (3:-1:1))
sim = test_convergence(dts, prob_ode_bigfloat2Dlinear, IRKGL16())
@test abs(sim.𝒪est[:final] - 16) < 0.5

# Backward integrations tests

#Define the problem
const g = 9.81
L = 1.0
u0 = [0, pi / 2]
function simplependulum(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
end

# adaptive=true

tspan = (6.3, 2.0)
prob = ODEProblem(simplependulum, u0, tspan)
sol = solve(prob, IRKGL16(), reltol = 1e-14, abstol = 1e-14)

@test sol.t[end] == tspan[2]

# adaptive=false

tspan = (6.3, 2.0)
dt0 = -0.1
prob = ODEProblem(simplependulum, u0, tspan)
sol = solve(prob, IRKGL16(), dt = dt0, adaptive = false)

@test sol.t[end] == tspan[2]
