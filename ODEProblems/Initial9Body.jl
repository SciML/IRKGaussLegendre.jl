function Initial9Body(T = Float64)
    GmSun = parse(BigFloat, "0.295912208285591100e-3")
    GmMercury = parse(BigFloat, "0.491248045036476000e-10")
    GmVenus = parse(BigFloat, "0.724345233264412000e-9")
    GmEarth = parse(BigFloat, "0.888769244512563400e-9")
    GmMoon = parse(BigFloat, "0.109318945074237400e-10")
    GmMars = parse(BigFloat, "0.954954869555077000e-10")

    Gm = [GmSun, GmMercury, GmVenus,
        GmEarth + GmMoon, GmMars,
        parse(BigFloat, "0.282534584083387000e-6"),
        parse(BigFloat, "0.845970607324503000e-7"),
        parse(BigFloat, "0.129202482578296000e-7"),
        parse(BigFloat, "0.152435734788511000e-7")]

    q = [
        parse(BigFloat, "0.00450250878464055477"),
        parse(BigFloat, "0.00076707642709100705"),
        parse(BigFloat, "0.00026605791776697764"), parse(
            BigFloat, "0.36176271656028195477"),
        parse(BigFloat, "-0.09078197215676599295"),
        parse(BigFloat, "-0.08571497256275117236"), parse(
            BigFloat, "0.61275194083507215477"),
        parse(BigFloat, "-0.34836536903362219295"),
        parse(BigFloat, "-0.19527828667594382236"), parse(
            BigFloat, "0.12051741410138465477"),
        parse(BigFloat, "-0.92583847476914859295"),
        parse(BigFloat, "-0.40154022645315222236"), parse(
            BigFloat, "-0.11018607714879824523"),
        parse(BigFloat, "-1.32759945030298299295"),
        parse(BigFloat, "-0.60588914048429142236"), parse(
            BigFloat, "-5.37970676855393644523"),
        parse(BigFloat, "-0.83048132656339789295"),
        parse(BigFloat, "-0.22482887442656542236"), parse(
            BigFloat, "7.89439068290953155477"),
        parse(BigFloat, "4.59647805517127300705"),
        parse(BigFloat, "1.55869584283189997764"), parse(
            BigFloat, "-18.26540225387235944523"),
        parse(BigFloat, "-1.16195541867586999295"),
        parse(BigFloat, "-0.25010605772133802236"), parse(
            BigFloat, "-16.05503578023336944523"),
        parse(BigFloat, "-23.94219155985470899295"),
        parse(BigFloat, "-9.40015796880239402236")
    ]

    v = [
        parse(BigFloat, "-0.00000035174953607552"),
        parse(BigFloat, "0.00000517762640983341"),
        parse(BigFloat, "0.00000222910217891203"), parse(
            BigFloat, "0.00336749397200575848"),
        parse(BigFloat, "0.02489452055768343341"),
        parse(BigFloat, "0.01294630040970409203"), parse(
            BigFloat, "0.01095206842352823448"),
        parse(BigFloat, "0.01561768426786768341"),
        parse(BigFloat, "0.00633110570297786403"), parse(
            BigFloat, "0.01681126830978379448"),
        parse(BigFloat, "0.00174830923073434441"),
        parse(BigFloat, "0.00075820289738312913"), parse(
            BigFloat, "0.01448165305704756448"),
        parse(BigFloat, "0.00024246307683646861"),
        parse(BigFloat, "-0.00028152072792433877"), parse(
            BigFloat, "0.00109201259423733748"),
        parse(BigFloat, "-0.00651811661280738459"),
        parse(BigFloat, "-0.00282078276229867897"), parse(
            BigFloat, "-0.00321755651650091552"),
        parse(BigFloat, "0.00433581034174662541"),
        parse(BigFloat, "0.00192864631686015503"), parse(
            BigFloat, "0.00022119039101561468"),
        parse(BigFloat, "-0.00376247500810884459"),
        parse(BigFloat, "-0.00165101502742994997"), parse(
            BigFloat, "0.00264276984798005548"),
        parse(BigFloat, "-0.00149831255054097759"),
        parse(BigFloat, "-0.00067904196080291327")]

    qq = copy(q)
    vv = copy(v)

    q0 = reshape(qq, 3, :)
    v0 = reshape(vv, 3, :)

    q0bar = [sum(Gm .* q0[j, :]) / sum(Gm) for j in 1:3]
    v0bar = [sum(Gm .* v0[j, :]) / sum(Gm) for j in 1:3]

    q0 = q0 .- q0bar
    v0 = v0 .- v0bar

    Gm = Gm[1:9]
    N = length(Gm)

    u0 = Array{T}(undef, 3, N, 2)
    u0[:, :, 2] .= v0[:, 1:N]
    u0[:, :, 1] .= q0[:, 1:N]

    Gm0 = Array{T}(undef, N)
    Gm0 .= Gm

    bodylist = ["Sun" "Mercury" "Venus" "EMB" "Mars" "Jupiter" "Saturn" "Uranus" "Neptune"]

    return u0, Gm0, bodylist
end
