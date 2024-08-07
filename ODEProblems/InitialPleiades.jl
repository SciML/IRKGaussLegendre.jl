
function InitialPleiades(T)
    N = 7
    Gm = T.([1, 2, 3, 4, 5, 6, 7])

    q = [
        parse(BigFloat, "3.0"),
        parse(BigFloat, "3.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "3.0"),
        parse(BigFloat, "-3.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "-1.0"),
        parse(BigFloat, "2.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "-3.0"),
        parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "2.0"),
        parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "-2.0"),
        parse(BigFloat, "-4.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "2.0"),
        parse(BigFloat, "4.0"),
        parse(BigFloat, "0.0")]

    v = [
        parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "0.0"),
        parse(BigFloat, "-1.25"),
        parse(BigFloat, "0.0"), parse(BigFloat, "0.0"),
        parse(BigFloat, "1.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "1.75"),
        parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "-1.5"),
        parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0")
    ]

    q0 = reshape(q, 3, :)
    v0 = reshape(v, 3, :)

    u0 = Array{T}(undef, 3, N, 2)
    u0[:, :, 2] .= convert.(T, v0)
    u0[:, :, 1] .= convert.(T, q0)

    return u0, Gm
end
