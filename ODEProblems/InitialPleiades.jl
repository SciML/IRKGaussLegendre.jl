
function InitialPleiades(T)
    q = [parse(BigFloat, "3.0"), parse(BigFloat, "3.0"), parse(BigFloat, "-1.0"),
        parse(BigFloat, "-3.0"), parse(BigFloat, "2.0"), parse(BigFloat, "-2.0"),
        parse(BigFloat, "2.0"), parse(BigFloat, "3.0"), parse(BigFloat, "-3.0"),
        parse(BigFloat, "2.0"), parse(BigFloat, "0.0"), parse(BigFloat, "0.0"),
        parse(BigFloat, "-4.0"), parse(BigFloat, "4.0")]

    v = [parse(BigFloat, "0.0"), parse(BigFloat, "0.0"), parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "0.0"), parse(BigFloat, "1.75"),
        parse(BigFloat, "-1.5"), parse(BigFloat, "0.0"), parse(BigFloat, "0.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "-1.25"), parse(BigFloat, "1.0"),
        parse(BigFloat, "0.0"), parse(BigFloat, "0.0")]

    u0 = Array{T}(undef, 1)

    u0 = convert.(T, [q; v])

    return u0
end
