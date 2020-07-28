#
#  Hénon-Heiles
#

function HenonHam(u,params)

    @inbounds begin
    q = view(u,1:2)
    p = view(u,3:4)

    T=zero(eltype(u))
    V=zero(eltype(u))

    T=1/2*(p[1]^2+p[2]^2)
    V=1/2*(q[1]^2+q[2]^2)+q[1]^2*q[2]-1/3*q[2]^3

    return(T+V)
    end

end

#
#  Dynamical Pro
#


function dotq(dq,q,p,params,t)
    @inbounds begin
    dq[1] = p[1]
    dq[2] = p[2]
    end
end

function dotp(dp,q,p,params,t)
    @inbounds begin
    dp[1] = -q[1]-2*q[1]*q[2]
    dp[2] = -q[2]-q[1]^2+q[2]^2
    end
end


#
#  Second Order Problem
#

function f2nd!(ddu,du,u,p,t)

  @inbounds begin
  q = view(u,1:2)   # x

  ddu[1] = -q[1]-2*q[1]*q[2]
  ddu[2] = -q[2]-q[1]^2+q[2]^2

  end
end
