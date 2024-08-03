#
# Equations for DynamicalODEProblem
#
function NbodyODEv!(dv, q, v, Gm, t)
    #
    #    dotv
    #
    N = length(Gm)

    for i in 1:3, j in 1:N
        dv[i, j] = 0
    end

    for i in 1:N
        xi = q[1, i]
        yi = q[2, i]
        zi = q[3, i]
        Gmi = Gm[i]
        for j in (i + 1):N
            xij = xi - q[1, j]
            yij = yi - q[2, j]
            zij = zi - q[3, j]
            Gmj = Gm[j]
            dotij = (xij * xij + yij * yij + zij * zij)
            auxij = 1 / (sqrt(dotij) * dotij)
            Gmjauxij = Gmj * auxij
            dv[1, i] -= Gmjauxij * xij
            dv[2, i] -= Gmjauxij * yij
            dv[3, i] -= Gmjauxij * zij
            Gmiauxij = Gmi * auxij
            dv[1, j] += Gmiauxij * xij
            dv[2, j] += Gmiauxij * yij
            dv[3, j] += Gmiauxij * zij
        end
    end

    return nothing
end

function NbodyODEq!(dq, q, v, Gm, t)
    #
    #    dotq
    #
    N = length(Gm)
    for i in 1:3, j in 1:N
        dq[i, j] = v[i, j]
    end

    return nothing
end

function solcopyNbody(u0, sol, prob, alg)
    #    
    #   return a copy of sol, converting sol.u to original u0's format 
    #
    #    
    n = length(sol.u)
    ttype = eltype(sol.t)
    uu = Vector{typeof(u0)}(undef, n)
    tt = Vector{ttype}(undef, n)
    uk = copy(u0)

    for k in 1:n
        uk[:, :, 1] = sol.u[k].x[1]
        uk[:, :, 2] = sol.u[k].x[2]
        uu[k] = copy(uk)
        tt[k] = sol.t[k]
    end

    sol_new = DiffEqBase.build_solution(
        prob, alg, tt, uu, stats = sol.stats, retcode = sol.retcode)

    return sol_new
end
