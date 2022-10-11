using IRKGaussLegendre, ODEProblemLibrary, DiffEqDevTools, Test
import ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear

function NbodyODE!(du, u, Gm, t)
    N = length(Gm)
    du[1, :, :] .= 0
    for i in 1:N
        qi = u[2, :, i]
        Gmi = Gm[i]
        du[2, :, i] = u[1, :, i]
        for j in (i + 1):N
            qj = u[2, :, j]
            Gmj = Gm[j]
            qij = qi - qj
            auxij = (qij[1] * qij[1] + qij[2] * qij[2] + qij[3] * qij[3])^(-3 / 2)
            du[1, :, i] -= Gmj * auxij * qij
            du[1, :, j] += Gmi * auxij * qij
        end
    end

    return
end

Gm = [5, 4, 3]
N = length(Gm)
q = [1, -1, 0, -2, -1, 0, 1, 3, 0]
v = zeros(size(q))
q0 = reshape(q, 3, :)
v0 = reshape(v, 3, :)
u0 = Array{Float64}(undef, 2, 3, N)
u0[1, :, :] = v0
u0[2, :, :] = q0
tspan = (0.0, 63.0)
prob = ODEProblem(NbodyODE!, u0, tspan, Gm)
sol1 = solve(prob, IRKGL16(), reltol = 1e-12, abstol = 1e-12)

# Analytical tests

sol = solve(prob_ode_2Dlinear, IRKGL16())
@test sol.errors[:l2] < 1e-16

dts = (1 // 2) .^ (4:-1:2)
sim = test_convergence(dts, prob_ode_bigfloat2Dlinear, IRKGL16())
@test abs(sim.𝒪est[:final] - 16) < 0.5
