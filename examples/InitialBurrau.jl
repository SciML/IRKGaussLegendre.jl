function InitialBurrau(T=Float64)

    Gm = [5, 4, 3]

    N=length(Gm)
    x1 = 1
    y1 = -1
    z1 = 0
    x2 = -2
    y2 = -1
    z2 = 0
    x3 = 1
    y3 = 3
    z3 = 0

    q =  [x1,y1,z1,x2,y2,z2,x3,y3,z3]
    v = zeros(size(q))
    q0 = reshape(q,3,:)
    v0 = reshape(v,3,:)
    u0 = Array{T}(undef,2,3,N)
    u0[1,:,:] = v0
    u0[2,:,:] = q0

    return u0, Gm

end
