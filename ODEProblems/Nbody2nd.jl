
function NbodyODE2nd!(ddu, du, u, Gm, t)
    N = length(Gm)

    for i in 1:N
        for k in 1:3
            ddu[k, i] = 0
        end
    end

    for i in 1:N
        xi = u[1, i]
        yi = u[2, i]
        zi = u[3, i]
        Gmi = Gm[i]
        for j in (i + 1):N
            xij = xi - u[1, j]
            yij = yi - u[2, j]
            zij = zi - u[3, j]
            Gmj = Gm[j]
            dotij = (xij * xij + yij * yij + zij * zij)
            auxij = 1 / (sqrt(dotij) * dotij)
            Gmjauxij = Gmj * auxij
            ddu[1, i] -= Gmjauxij * xij
            ddu[2, i] -= Gmjauxij * yij
            ddu[3, i] -= Gmjauxij * zij
            Gmiauxij = Gmi * auxij
            ddu[1, j] += Gmiauxij * xij
            ddu[2, j] += Gmiauxij * yij
            ddu[3, j] += Gmiauxij * zij
        end
    end

    return nothing
end
