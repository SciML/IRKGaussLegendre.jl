function InitialBurrau(T=Float64)

    Gm = [parse(BigFloat,"5"), parse(BigFloat,"4"),parse(BigFloat,"3")]

    N=length(Gm)
    x1 = parse(BigFloat,"1")
    y1 = parse(BigFloat,"-1")
    z1 = parse(BigFloat,"0")
    x2 = parse(BigFloat,"-2")
    y2 = parse(BigFloat,"-1")
    z2 = parse(BigFloat,"0")
    x3 = parse(BigFloat,"1")
    y3 = parse(BigFloat,"3")
    z3 = parse(BigFloat,"0")

    q =  [x1,y1,z1,x2,y2,z2,x3,y3,z3]
    v = zeros(size(q))
    q0 = reshape(q,3,:)
    v0 = reshape(v,3,:)
    u0 = Array{T}(undef,2,3,N)
    u0[1,:,:] = convert.(T,v0)
    u0[2,:,:] = convert.(T,q0)

    return u0, convert.(T,Gm)

end

function InitialBurrau2(T=Float64)

    Gm = [parse(BigFloat,"5"), parse(BigFloat,"4"),parse(BigFloat,"3")]

    N=length(Gm)
    x1 = parse(BigFloat,"1")
    y1 = parse(BigFloat,"-1")
    z1 = parse(BigFloat,"0")
    x2 = parse(BigFloat,"-2")
    y2 = parse(BigFloat,"-1")
    z2 = parse(BigFloat,"0")
    x3 = parse(BigFloat,"1")
    y3 = parse(BigFloat,"3")
    z3 = parse(BigFloat,"0")

    q =  [x1,y1,z1,x2,y2,z2,x3,y3,z3]
    v = zeros(size(q))
    q0 = reshape(q,3,:)
    v0 = reshape(v,3,:)
    u0 = Array{T}(undef,2,3,N)
    u0[2,:,:] = convert.(T,v0)
    u0[1,:,:] = convert.(T,q0)

    return u0, convert.(T,Gm)

end
