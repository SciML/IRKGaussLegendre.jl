#
#  Initial values from:
#         "Achieving Brouwer law with implicit Runge-Kutta methods"
#                 E. Hairer, R.I. McLachlan and A. Razakarivony
#
#          (solution is chaotic)

function InitialHenon(T = Float64)
    q = [parse(BigFloat, "0.0"), parse(BigFloat, "0.3")]

    p = [parse(BigFloat, "0.0"), parse(BigFloat, "0.2")]

    H0 = 1 / 8
    V0 = V = 1 / 2 * (q[1]^2 + q[2]^2) + q[1]^2 * q[2] - 1 / 3 * q[2]^3
    p[1] = sqrt(2 * H0 - 2 * V0 - p[2]^2)

    u0 = convert.(T, [q; p])

    return u0
end
