function NbodyEnergy(u, Gm)
    N = length(Gm)
    zerouel = zero(eltype(u))
    T = zerouel
    U = zerouel
    for i = 1:N
        qi = u[2, :, i]
        vi = u[1, :, i]
        Gmi = Gm[i]
        T += Gmi * (vi[1] * vi[1] + vi[2] * vi[2] + vi[3] * vi[3])
        for j = (i+1):N
            qj = u[2, :, j]  # qj = u[2,:,j]
            Gmj = Gm[j]
            qij = qi - qj
            U -= Gmi * Gmj / norm(qij)
        end
    end
    1 / 2 * T + U
end

function NbodyEnergy2(u, Gm)
    #
    #   special implementantion: considering:
    #              u[1,:,:]: position
    #              u[2,:,:]: velocity
    #
    N = length(Gm)
    zerouel = zero(eltype(u))
    T = zerouel
    U = zerouel
    for i = 1:N
        qi = u[1, :, i]    # qi = u[2,:,i]
        vi = u[2, :, i]    # vi = u[1,:,i]
        Gmi = Gm[i]
        T += Gmi * (vi[1] * vi[1] + vi[2] * vi[2] + vi[3] * vi[3])
        for j = (i+1):N
            qj = u[1, :, j]  # qj = u[2,:,j]
            Gmj = Gm[j]
            qij = qi - qj
            U -= Gmi * Gmj / norm(qij)
        end
    end
    1 / 2 * T + U
end

# OdeProblem

function NbodyODE!(du, u, Gm, t)
    N = length(Gm)
    du[1, :, :] .= zero(eltype(u))
    for i = 1:N
        qi = u[2, :, i]
        Gmi = Gm[i]
        du[2, :, i] = u[1, :, i]
        for j = (i+1):N
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

# DynamicalODEProblem

function NbodyODEv!(dv, q, v, Gm, t)
    #
    #    dotv
    #
    N = length(Gm)
    dv[:, :] .= zero(eltype(q))
    for i = 1:N
        qi = q[:, i]
        Gmi = Gm[i]
        for j = (i+1):N
            qj = q[:, j]
            Gmj = Gm[j]
            qij = qi - qj
            auxij = (qij[1] * qij[1] + qij[2] * qij[2] + qij[3] * qij[3])^(-3 / 2)
            dv[:, i] -= Gmj * auxij * qij
            dv[:, j] += Gmi * auxij * qij
        end
    end

    return
end

function NbodyODEq!(dq, q, v, Gm, t)
    #
    #    dotq
    #
    N = length(Gm)
    for i = 1:N
        dq[:, i] = v[:, i]
    end

    return
end

## 2ndOrderProblem

function NbodyODE2nd!(ddu, du, u, Gm, t)
    N = length(Gm)
    ddu[:, :] .= zero(eltype(u))
    for i = 1:N
        qi = u[:, i]
        Gmi = Gm[i]
        for j = (i+1):N
            qj = u[:, j]
            Gmj = Gm[j]
            qij = qi - qj
            auxij = (qij[1] * qij[1] + qij[2] * qij[2] + qij[3] * qij[3])^(-3 / 2)
            ddu[:, i] -= Gmj * auxij * qij
            ddu[:, j] += Gmi * auxij * qij
        end
    end

    return
end
