function Initial16Body(T = Float64)
    nbody = 16
    iE = nbody - 1  # Earth indice
    iM = nbody    # Moon  indice

    # ********************    Masses  ******************************* #

    GmMoon = parse(BigFloat, "0.109318945074237400e-10")
    GmEarth = parse(BigFloat, "0.888769244512563400e-9")

    Gm = [parse(BigFloat, "0.295912208285591100e-3"),
        parse(BigFloat, "0.491248045036476000e-10"),
        parse(BigFloat, "0.724345233264412000e-9"),
        parse(BigFloat, "0.954954869555077000e-10"),
        parse(BigFloat, "0.282534584083387000e-6"),
        parse(BigFloat, "0.845970607324503000e-7"),
        parse(BigFloat, "0.129202482578296000e-7"),
        parse(BigFloat, "0.152435734788511000e-7"),
        parse(BigFloat, "0.217844105197418000e-11"),
        #Ceres, Pallas, Vesta, Iris, Bamberga
        parse(BigFloat, "0.140047655617234400e-12"),
        parse(BigFloat, "0.310444819893871300e-13"),
        parse(BigFloat, "0.385475018780881000e-13"),
        parse(BigFloat, "0.213643444257140700e-14"),
        parse(BigFloat, "0.138862658985619900e-14"),
        # Earth and Moon
        GmMoon,
        GmEarth]

    # ********************    Position  ******************************* #

    #    Table 13. Initial positions (au) and velocities (au/day) of the asteroids with respect to the Sun at
    #Julian day (TDB) 2440400.5 (June 28, 1969) in the ICRF2 frame (1 of 15).

    Sunq = [parse(BigFloat, "0.00450250878464055477"),
        parse(BigFloat, "0.00076707642709100705"),
        parse(BigFloat, "0.00026605791776697764")]

    Ceresq = Sunq + [parse(BigFloat, "1.438681809676469747"),
        parse(BigFloat, "-2.204373633189407045"),
        parse(BigFloat, "-1.326397853361325874")]

    Pallasq = Sunq + [parse(BigFloat, "0.203832272462290465"),
        parse(BigFloat, "-3.209619436062307152"),
        parse(BigFloat, "0.623843179079393351")]

    Vestaq = Sunq + [parse(BigFloat, "0.182371836377417107"),
        parse(BigFloat, "2.386628211277654010"),
        parse(BigFloat, "0.924596062836265498")]

    Irisq = Sunq + [parse(BigFloat, "1.892475267790300286"),
        parse(BigFloat, "-0.848414748075139946"),
        parse(BigFloat, "-0.157159319044464590")]

    Bambergaq = Sunq + [parse(BigFloat, "1.398759064223541682"),
        parse(BigFloat, "-1.287476729008325105"),
        parse(BigFloat, "-0.669098428660833799")]

    asteroidsq = vcat(Ceresq, Pallasq, Vestaq, Irisq, Bambergaq)

    #Table 6. Initial position (au) and velocity (au/day) of the Moon at Julian day (TDB) 2440400.5
    #(June 28, 1969), given with respect to Earth in the ICRF2 frame.        

    Moonq_Earth = [
        parse(BigFloat, "-0.00080817735147818490"),
        parse(BigFloat, "-0.0019946299854970130"),
        parse(BigFloat, "-0.00108726268307068900")]

    Moonv_Earth = [
        parse(BigFloat, "0.00060108481561422370"),
        parse(BigFloat, "-0.00016744546915764980"),
        parse(BigFloat, "-0.00008556214140094871")]

    #Table 5. Initial positions (au) and velocities (au/day) of the Sun and planets at Julian day (TDB) 2440400.5
    #(June 28, 1969), given with respect to the integration origin in the ICRF2 frame.    

    q = vcat(
        [parse(BigFloat, "0.00450250878464055477"),
            parse(BigFloat, "0.00076707642709100705"),
            parse(BigFloat, "0.00026605791776697764"), parse(
                BigFloat, "0.36176271656028195477"),
            parse(BigFloat, "-0.09078197215676599295"),
            parse(BigFloat, "-0.08571497256275117236"), parse(
                BigFloat, "0.61275194083507215477"),
            parse(BigFloat, "-0.34836536903362219295"),
            parse(BigFloat, "-0.19527828667594382236"), parse(
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
            parse(BigFloat, "-9.40015796880239402236"), parse(
                BigFloat, "-30.48331376718383944523"),
            parse(BigFloat, "-0.87240555684104999295"),
            parse(BigFloat, "8.91157617249954997764")], asteroidsq,

        # EM Bary and Moon
        [parse(BigFloat, "0.12051741410138465477"),
            parse(BigFloat, "-0.92583847476914859295"),
            parse(BigFloat, "-0.40154022645315222236")],

        #Table 6. Initial position (au) and velocity (au/day) of the Moon at Julian day (TDB) 2440400.5
        #(June 28, 1969), given with respect to Earth in the ICRF2 frame.        

        Moonq_Earth)

    # ********************    Velocity  ******************************* #

    Sunv = [parse(BigFloat, "-0.00000035174953607552"),
        parse(BigFloat, "0.00000517762640983341"),
        parse(BigFloat, "0.00000222910217891203")]

    Ceresv = Sunv + [parse(BigFloat, "0.008465406136316316"),
        parse(BigFloat, "0.004684247977335608"),
        parse(BigFloat, "0.000466157738595739")]

    Pallasv = Sunv + [parse(BigFloat, "0.008534313855651248"),
        parse(BigFloat, "-0.000860659210123161"),
        parse(BigFloat, "-0.000392901992572746")]

    Vestav = Sunv + [parse(BigFloat, "-0.010174496747119257"),
        parse(BigFloat, "0.000041478190529952"),
        parse(BigFloat, "0.001344157634155624")]

    Irisv = Sunv + [parse(BigFloat, "0.002786950314570632"),
        parse(BigFloat, "0.011314057384917047"),
        parse(BigFloat, "0.004975132577079665")]

    Bambergav = Sunv + [parse(BigFloat, "0.007164363244556328"),
        parse(BigFloat, "0.009219958777618218"),
        parse(BigFloat, "0.006857861727407507")]

    asteroidsv = vcat(Ceresv, Pallasv, Vestav, Irisv, Bambergav)

    v = vcat(
        [parse(BigFloat, "-0.00000035174953607552"),
            parse(BigFloat, "0.00000517762640983341"),
            parse(BigFloat, "0.00000222910217891203"), parse(
                BigFloat, "0.00336749397200575848"),
            parse(BigFloat, "0.02489452055768343341"),
            parse(BigFloat, "0.01294630040970409203"), parse(
                BigFloat, "0.01095206842352823448"),
            parse(BigFloat, "0.01561768426786768341"),
            parse(BigFloat, "0.00633110570297786403"), parse(
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
            parse(BigFloat, "-0.00067904196080291327"), parse(
                BigFloat, "0.00032220737349778078"),
            parse(BigFloat, "-0.00314357639364532859"),
            parse(BigFloat, "-0.00107794975959731297")], asteroidsv,

        #EM Bary and Moon*)
        [parse(BigFloat, "0.01681126830978379448"),
            parse(BigFloat, "0.00174830923073434441"),
            parse(BigFloat, "0.00075820289738312913")], Moonv_Earth)

    N = length(Gm)

    qq = copy(q)
    vv = copy(v)

    qEM = q[43:45]
    qEarth = qEM - GmMoon / (GmEarth + GmMoon) * Moonq_Earth
    qMoon = qEarth + Moonq_Earth
    qq[43:45] = qEarth
    qq[46:48] = qMoon

    vEM = v[43:45]
    vEarth = vEM - GmMoon / (GmEarth + GmMoon) * Moonv_Earth
    vMoon = vEarth + Moonv_Earth
    vv[43:45] = vEarth
    vv[46:48] = vMoon

    q0 = reshape(qq, 3, :)
    v0 = reshape(vv, 3, :)

    q0bar = [sum(Gm .* q0[j, :]) / sum(Gm) for j in 1:3]
    v0bar = [sum(Gm .* v0[j, :]) / sum(Gm) for j in 1:3]

    q0 = q0 .- q0bar
    v0 = v0 .- v0bar

    u0 = Array{T}(undef, 3, N, 2)
    u0[:, :, 2] = v0
    u0[:, :, 1] = q0

    Gm0 = Array{T}(undef, N)
    Gm0 .= Gm

    bodylist = ["Sun" "Mercury" "Venus" "Mars" "Jupiter" "Saturn" "Uranus" "Neptune" "Pluto" "Ceres" "Pallas" "Vesta" "Iris" "Bamberga" "Earth" "Moon"]

    return u0, Gm0, bodylist
end
