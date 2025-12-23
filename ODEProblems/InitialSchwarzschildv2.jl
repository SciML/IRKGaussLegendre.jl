
#
#  Initial values for Schwarzschild problem
#  (2025-07-09 cartesian coordinates)
#  

function InitialSchwarzschildv2(T = Float64)
    params = BigFloat[0.995, 4.6, 8.9e-5]

    r = 11
    θ = pi/2
    pr = 0

    u0 = BigFloat[r, θ, pr, 0]

    pθ = r * sqrt(-1 - 2*Ham_Schwarzschild(u0, params))
    params = convert.(T, params)

    c=cos(θ)
    s=sin(θ)
    x=r*c
    y=r*s
    px=c*pr-s*pθ/r
    py=s*pr+c*pθ/r
    u0 = convert.(T, [x, y, px, py])

    return u0, params
end
