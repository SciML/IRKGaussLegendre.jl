function NbodyEnergy(u,Gm)
     N = length(Gm)
     zerouel = zero(eltype(u))
     T = zerouel
     U = zerouel
     for i in 1:N
        qi = u[2,:,i]   #  x = view(u,1:7)   # x
        vi = u[1,:,i]
        Gmi = Gm[i]
        T += Gmi*(vi[1]*vi[1]+vi[2]*vi[2]+vi[3]*vi[3])
        for j in (i+1):N
           qj = u[2,:,j]
           Gmj = Gm[j]
           qij = qi - qj
           U -= Gmi*Gmj/norm(qij)
        end
     end
    1/2*T + U
end


function NbodyODE!(du,u,Gm,t)
     N = length(Gm)
     du[1,:,:] .= 0
     for i in 1:N
        qi = u[2,:,i]
        Gmi = Gm[i]
        du[2,:,i] = u[1,:,i]
        for j in (i+1):N
           qj = u[2,:,j]
           Gmj = Gm[j]
           qij = qi - qj
           auxij = (qij[1]*qij[1]+qij[2]*qij[2]+qij[3]*qij[3])^(-3/2)
           du[1,:,i] -= Gmj*auxij*qij
           du[1,:,j] += Gmi*auxij*qij
        end
     end

    return
    
end
