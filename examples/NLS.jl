#
#  NLS: Non Linear Schodinger
#

function NLSHam(u,p)

    N = 5

    @inbounds begin
    q = view(u,1:N)
    p = view(u,N+1:2*N)

    H1=zero(eltype(u))
    H2=zero(eltype(u))

    for i in 1:N
        H1+=1/4*(q[i]^2+p[i]^2)^2
    end

    for i in 2:N
        H2+=(p[i-1]^2*p[i]^2+q[i-1]^2*q[i]^2-q[i-1]^2*p[i]^2-p[i-1]^2*q[i]^2+
             4*p[i-1]*p[i]*q[i-1]*q[i])
    end

    return(H1-H2)
    end

end


function NLSODE!(du,u,p,t)

    N = 5

    q = view(u,1:N)
    p = view(u,N+1:2*N)


    du[1]=p[1]*(q[1]^2+p[1]^2)-(2*p[1]*p[2]^2-2*p[1]*q[2]^2+4*p[2]*q[1]*q[2])

    du[N+1]=-q[1]*(q[1]^2+p[1]^2)+(2*q[1]*q[2]^2-2*q[1]*p[2]^2+4*p[1]*p[2]*q[2])


    for i in 2:N-1

        du[i]=p[i]*(q[i]^2+p[i]^2)-(2*p[i-1]^2*p[i]-2*q[i-1]^2*p[i]+4*p[i-1]*q[i-1]*q[i])-(2*p[i]*p[i+1]^2-2*p[i]*q[i+1]^2+4*p[i+1]*q[i]*q[i+1])

        du[i+N]=-q[i]*(q[i]^2+p[i]^2)+(2*q[i-1]^2*q[i]-2*p[i-1]^2*q[i]+4*p[i-1]*p[i]*q[i-1])+(2*q[i]*q[i+1]^2-2*q[i]*p[i+1]^2+4*p[i]*p[i+1]*q[i+1])

   end

   du[N]=p[N]*(q[N]^2+p[N]^2)-(2*p[N-1]^2*p[N]-2*q[N-1]^2*p[N]+4*p[N-1]*q[N-1]*q[N])

   du[2*N]=-q[N]*(q[N]^2+p[N]^2)+(2*q[N-1]^2*q[N]-2*p[N-1]^2*q[N]+4*p[N-1]*p[N]*q[N-1])


end
