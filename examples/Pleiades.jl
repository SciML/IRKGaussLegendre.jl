
# Energy

function NbodyEnergy(u, Gm)
"""
     Nbody problem Hamiltonian (Cartesian Coordinates)
"""

    # Declarations

    dim=2
    nbody=length(Gm)

    # Implementation

 @inbounds begin
    x = view(u,1:7)   # x
    y = view(u,8:14)  # y
    v = view(u,15:21) # x′
    w = view(u,22:28) # y′

    H=0.
    P=0.

    for i in 1:nbody
        H+=Gm[i]*(v[i]*v[i]+w[i]*w[i])
        for j in i+1:nbody
            r = ((x[i]-x[j])^2+(y[i]-y[j])^2)^(1/2)
            P+=(Gm[i]/r)*Gm[j]
        end
    end

    return(H/2-P)
    end

end


# OdeProblem


function f(du,u,p,t)
  @inbounds begin
  x = view(u,1:7)   # x
  y = view(u,8:14)  # y
  v = view(u,15:21) # x′
  w = view(u,22:28) # y′
  du[1:7] .= v
  du[8:14].= w
  for i in 15:28
    du[i] = zero(u[1])
  end
  for i=1:7,j=1:7
    if i != j
      r = ((x[i]-x[j])^2 + (y[i] - y[j])^2)^(3/2)
      du[14+i] += j*(x[j] - x[i])/r
      du[21+i] += j*(y[j] - y[i])/r
    end
  end
  end
end


# DynamicalOdeProblem


function dotv(dv,q,v,par,t)
@inbounds begin
  x = view(q,1:7)   # x
  y = view(q,8:14)  # y
  vx = view(v,1:7)   # x′
  vy = view(v,8:14)  # y′
  for i in 1:14
    dv[i] = zero(x[1])
  end
  for i=1:7,j=1:7
    if i != j
      r = ((x[i]-x[j])^2 + (y[i] - y[j])^2)^(3/2)
      dv[i] += j*(x[j] - x[i])/r
      dv[7+i] += j*(y[j] - y[i])/r
    end
  end
  end
end


function dotq(dq,q,v,par,t)
@inbounds begin
  x = view(q,1:7)   # x
  y = view(q,8:14)  # y
  vx = view(v,1:7)   # x′
  vy = view(v,8:14)  # y′
  dq[1:7] .= vx
  dq[8:14].= vy
  end
end
