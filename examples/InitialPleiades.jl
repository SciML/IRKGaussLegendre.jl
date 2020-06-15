function InitialPleiades(T=Float64)

    u0=convert.(T,[3.0,3.0,-1.0,-3.0,2.0,-2.0,2.0,
        3.0,-3.0,2.0,0,0,-4.0,4.0,
        0,0,0,0,0,1.75,-1.5,
        0,0,0,-1.25,1,0,0])

    Gm=convert.(T,[1,2,3,4,5,6,7])
    return u0,Gm

end
