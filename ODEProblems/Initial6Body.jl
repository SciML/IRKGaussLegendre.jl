include("Initial10Body.jl")

function Initial6Body(T = Float64)

    # DE430
    N = 6

    u_, Gm_, bodylist_ = Initial10Body()

    # Recalcular datos del Sistema, tomando el Sistema Solar interior como una masa puntual en su centro de masas

    q0 = zeros(3)
    v0 = zeros(3)
    Gm0 = 0.0
    for i in 1:5
        Gmi = Gm_[i]
        qi = u_[:, i, 1]
        vi = u_[:, i, 2]
        q0 += Gmi * qi
        v0 += Gmi * vi
        Gm0 += Gmi
    end
    q0 /= Gm0
    v0 /= Gm0

    u0 = Array{T}(undef, 3, N, 2)
    u0[:, 1, 1] = q0
    u0[:, 1, 2] = v0
    u0[:, 2:end, 1] = u_[:, 6:10, 1]
    u0[:, 2:end, 2] = u_[:, 6:10, 2]

    Gm = [Gm0; Gm_[6:10]]
    bodylist = ["Sun" "Jupiter" "Saturn" "Uranus" "Neptune" "Pluto"]

    return u0, Gm, bodylist
end
