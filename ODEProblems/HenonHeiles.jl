function Potential_U(q1, q2)
    return 1 / 2 * (q1 * q1 + q2 * q2) + q1 * q1 * q2 - 1 / 3 * q2 * q2 * q2
end

function HenonHeliesHam(u, params)
    q1 = u[1]
    q2 = u[2]
    p1 = u[3]
    p2 = u[4]

    U = Potential_U(q1, q2)

    return 1 / 2 * (p1 * p1 + p2 * p2) + U
end

function HenonHeliesHam(u)
    return HenonHeliesHam(u, [])
end

function H1(u)
    p1 = u[3]
    p2 = u[4]

    return 1 / 2 * (p1 * p1 + p2 * p2)
end

function H2(u)
    q1 = u[1]
    q2 = u[2]

    U = Potential_U(q1, q2)

    return U
end

function HenonHeilesODE!(du, u, params, t)
    q1 = u[1]
    q2 = u[2]
    p1 = u[3]
    p2 = u[4]

    du[1] = p1
    du[2] = p2

    du[3] = -q1 - 2 * q1 * q2
    du[4] = -q2 - q1 * q1 + q2 * q2

    return nothing
end

function flowH1HenonHeiles!(uj, ej, h, params)
    q1 = uj[1]
    q2 = uj[2]
    p1 = uj[3]
    p2 = uj[4]

    uj[1] = q1 + h * p1
    uj[2] = q2 + h * p2

    return nothing
end

function flowH2HenonHeiles!(uj, ej, h, params)
    q1 = uj[1]
    q2 = uj[2]
    p1 = uj[3]
    p2 = uj[4]

    uj[3] = p1 + h * (-q1 - 2 * q1 * q2)
    uj[4] = p2 + h * (-q2 - q1 * q1 + q2 * q2)

    return nothing
end

#
# functions for SecondOrderODEProblem
#

function HenonHeilesODE2nd!(ddu, du, u, par, t)
    q1 = u[1]
    q2 = u[2]
    ddu[1] = -q1 - 2 * q1 * q2
    ddu[2] = -q2 - q1 * q1 + q2 * q2

    return nothing
end

#
# functions for DynamicalODEProblem
#

function HenonHeilesODEv!(dv, q, v, par, t)
    q1 = q[1]
    q2 = q[2]

    dv[1] = -q1 - 2 * q1 * q2
    dv[2] = -q2 - q1 * q1 + q2 * q2

    return nothing
end

function HenonHeilesODEq!(dq, q, v, par, t)
    p1 = v[1]
    p2 = v[2]

    dq[1] = p1
    dq[2] = p2

    return nothing
end
