#
#  NLS: Non Linear Schrödinger
#

function NLSHam(u, par)
    N = 5

    @inbounds begin
        q = view(u, 1:N)
        p = view(u, (N + 1):(2 * N))

        H1 = zero(eltype(u))
        H2 = zero(eltype(u))

        for i in 1:N
            H1 += 1 / 4 * (q[i]^2 + p[i]^2)^2
        end

        for i in 2:N
            H2 += (
                p[i - 1]^2 * p[i]^2 + q[i - 1]^2 * q[i]^2 - q[i - 1]^2 * p[i]^2 -
                    p[i - 1]^2 * q[i]^2 + 4 * p[i - 1] * p[i] * q[i - 1] * q[i]
            )
        end

        return (H1 - H2)
    end
end

function NLSODE_!(du, u, par, t)
    #
    #   wrong for simd !!!
    #
    N = 5

    q = view(u, 1:N)
    p = view(u, (N + 1):(2 * N))

    du[1] = p[1] * (q[1]^2 + p[1]^2) -
        (2 * p[1] * p[2]^2 - 2 * p[1] * q[2]^2 + 4 * p[2] * q[1] * q[2])
    du[N + 1] = -q[1] * (q[1]^2 + p[1]^2) +
        (2 * q[1] * q[2]^2 - 2 * q[1] * p[2]^2 + 4 * p[1] * p[2] * q[2])

    for i in 2:(N - 1)
        du[i] = p[i] * (q[i]^2 + p[i]^2) -
            (
            2 * p[i - 1]^2 * p[i] - 2 * q[i - 1]^2 * p[i] +
                4 * p[i - 1] * q[i - 1] * q[i]
        ) -
            (
            2 * p[i] * p[i + 1]^2 - 2 * p[i] * q[i + 1]^2 +
                4 * p[i + 1] * q[i] * q[i + 1]
        )

        du[i + N] = -q[i] * (q[i]^2 + p[i]^2) +
            (
            2 * q[i - 1]^2 * q[i] - 2 * p[i - 1]^2 * q[i] +
                4 * p[i - 1] * p[i] * q[i - 1]
        ) +
            (
            2 * q[i] * q[i + 1]^2 - 2 * q[i] * p[i + 1]^2 +
                4 * p[i] * p[i + 1] * q[i + 1]
        )
    end

    du[N] = p[N] * (q[N]^2 + p[N]^2) -
        (2 * p[N - 1]^2 * p[N] - 2 * q[N - 1]^2 * p[N] + 4 * p[N - 1] * q[N - 1] * q[N])
    du[2 * N] = -q[N] * (q[N]^2 + p[N]^2) +
        (
        2 * q[N - 1]^2 * q[N] - 2 * p[N - 1]^2 * q[N] +
            4 * p[N - 1] * p[N] * q[N - 1]
    )

    return nothing
end

function NLSODE!(du, u, par, t)

    #  2024-07-30 New
    #   correct for simd
    #
    N = 5

    #q = view(u, 1:N)
    #p = view(u, (N + 1):(2 * N))

    iq = 0
    ip = N

    du[1] = u[ip + 1] * (u[iq + 1]^2 + u[ip + 1]^2) -
        (
        2 * u[ip + 1] * u[ip + 2]^2 - 2 * u[ip + 1] * u[iq + 2]^2 +
            4 * u[ip + 2] * u[iq + 1] * u[iq + 2]
    )
    du[N + 1] = -u[iq + 1] * (u[iq + 1]^2 + u[ip + 1]^2) +
        (
        2 * u[iq + 1] * u[iq + 2]^2 - 2 * u[iq + 1] * u[ip + 2]^2 +
            4 * u[ip + 1] * u[ip + 2] * u[iq + 2]
    )

    for i in 2:(N - 1)
        du[i] = u[ip + i] * (u[iq + i]^2 + u[ip + i]^2) -
            (
            2 * u[ip + (i - 1)]^2 * u[ip + i] - 2 * u[iq + (i - 1)]^2 * u[ip + i] +
                4 * u[ip + (i - 1)] * u[iq + (i - 1)] * u[iq + i]
        ) -
            (
            2 * u[ip + i] * u[ip + (i + 1)]^2 - 2 * u[ip + i] * u[iq + (i + 1)]^2 +
                4 * u[ip + (i + 1)] * u[iq + i] * u[iq + (i + 1)]
        )

        du[i + N] = -u[iq + i] * (u[iq + i]^2 + u[ip + i]^2) +
            (
            2 * u[iq + (i - 1)]^2 * u[iq + i] - 2 * u[ip + (i - 1)]^2 * u[iq + i] +
                4 * u[ip + (i - 1)] * u[ip + i] * u[iq + (i - 1)]
        ) +
            (
            2 * u[iq + i] * u[iq + (i + 1)]^2 - 2 * u[iq + i] * u[ip + (i + 1)]^2 +
                4 * u[ip + i] * u[ip + (i + 1)] * u[iq + (i + 1)]
        )
    end

    du[N] = u[ip + N] * (u[iq + N]^2 + u[ip + N]^2) -
        (
        2 * u[ip + (N - 1)]^2 * u[ip + N] - 2 * u[iq + (N - 1)]^2 * u[ip + N] +
            4 * u[ip + (N - 1)] * u[iq + (N - 1)] * u[iq + N]
    )
    du[2 * N] = -u[iq + N] * (u[iq + N]^2 + u[ip + N]^2) +
        (
        2 * u[iq + (N - 1)]^2 * u[iq + N] - 2 * u[ip + (N - 1)]^2 * u[iq + N] +
            4 * u[ip + (N - 1)] * u[ip + N] * u[iq + (N - 1)]
    )

    return nothing
end
