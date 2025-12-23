function NbodyEnergy(u, Gm)
    N = length(Gm)
    zerouel = zero(eltype(u))
    T = zerouel
    U = zerouel
    for i in 1:N
        qi = u[:, i, 1]
        vi = u[:, i, 2]
        Gmi = Gm[i]
        T += Gmi * (vi[1] * vi[1] + vi[2] * vi[2] + vi[3] * vi[3])
        for j in (i + 1):N
            qj = u[:, j, 1]
            Gmj = Gm[j]
            qij = qi - qj
            U -= Gmi * Gmj / norm(qij)
        end
    end
    1 / 2 * T + U
end

function NbodyEnergy3(u, Gm)
    #
    #   berria 2024-07-19
    #          Recursive arrays !!!    
    #              
    #
    N = length(Gm)
    zerouel = zero(eltype(u))
    T = zerouel
    U = zerouel
    q = u.x[1]
    v = u.x[2]
    for i in 1:N
        qi = q[:, i]
        vi = v[:, i]
        Gmi = Gm[i]
        T += Gmi * (vi[1] * vi[1] + vi[2] * vi[2] + vi[3] * vi[3])
        for j in (i + 1):N
            qj = q[:, j]  # qj = u[2,:,j]
            Gmj = Gm[j]
            qij = qi - qj
            U -= Gmi * Gmj / norm(qij)
        end
    end
    1 / 2 * T + U
end

function NbodyEnergy4(u, Gm)
    #
    #   berria 2024-07-19
    #          Array partition !!!    
    #              
    #
    N = length(Gm)
    zerouel = zero(eltype(u))
    T = zerouel
    U = zerouel
    q = u.x[2]
    v = u.x[1]
    for i in 1:N
        qi = q[:, i]
        vi = v[:, i]
        Gmi = Gm[i]
        T += Gmi * (vi[1] * vi[1] + vi[2] * vi[2] + vi[3] * vi[3])
        for j in (i + 1):N
            qj = q[:, j]  # qj = u[2,:,j]
            Gmj = Gm[j]
            qij = qi - qj
            U -= Gmi * Gmj / norm(qij)
        end
    end
    1 / 2 * T + U
end

###

function NbodyBarycenter(u, Gm)
    N = length(Gm)
    dim = size(u, 1)
    A = zeros(dim)
    B = zeros(dim)
    for i in 1:N
        qi = u[:, i, 1]
        vi = u[:, i, 2]
        Gmi = Gm[i]
        A += Gmi * qi
        B += Gmi * vi
    end
    return A, B
end

function NbodyODE!(F, u, Gm, t)
    N = length(Gm)
    for i in 1:N
        for k in 1:3
            F[k, i, 2] = 0
        end
    end
    for i in 1:(N - 1)
        xi = u[1, i, 1]
        yi = u[2, i, 1]
        zi = u[3, i, 1]
        Gmi = Gm[i]
        for j in (i + 1):N
            xij = xi - u[1, j, 1]
            yij = yi - u[2, j, 1]
            zij = zi - u[3, j, 1]
            Gmj = Gm[j]
            dotij = (xij * xij + yij * yij + zij * zij)
            auxij = 1 / (sqrt(dotij) * dotij)
            Gmjauxij = Gmj * auxij
            F[1, i, 2] -= Gmjauxij * xij
            F[2, i, 2] -= Gmjauxij * yij
            F[3, i, 2] -= Gmjauxij * zij
            Gmiauxij = Gmi * auxij
            F[1, j, 2] += Gmiauxij * xij
            F[2, j, 2] += Gmiauxij * yij
            F[3, j, 2] += Gmiauxij * zij
        end
    end
    for i in 1:3, j in 1:N

        F[i, j, 1] = u[i, j, 2]
    end
    return nothing
end
