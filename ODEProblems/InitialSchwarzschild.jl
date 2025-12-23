
#
#  Initial values for Schwarzschild problem
#  (polar coordinates)
#  

function InitialSchwarzschild(T = Float64)
    params = BigFloat[0.995, 4.6, 8.9e-5]

    r = 11
    θ = pi/2
    pr = 0
    u0 = BigFloat[r, θ, pr, 0]

    pθ = r * sqrt(-1 - 2*Ham_Schwarzschild(u0, params))

    params = convert.(T, params)

    u0 = convert.(T, [r, θ, pr, pθ])

    return u0, params
end
