#
#
#   Auxiliar functions:
#       ErrorEst:
#       MyNorm
#       Rdigits

function ErrorEst(U,F,dt,beta,abstol,reltol)

    (s,)=size(F)
	D=length(U[1])

	est=zero(typeof(dt))

    @inbounds begin
	for k in eachindex(U[1])

		sum=zero(typeof(dt))
		maxU=zero(typeof(dt))
		for is in 1:s
        	sum+=beta[is]*F[is][k]
			maxU=max(maxU,abs(U[is][k]))
        end

		est+=(abs(dt*sum))^2/(abstol+maxU*reltol)
	end
    end

    return(est/D)

end


function MyNorm(u,abstol,reltol)

	norm=zero(eltype(u))

    @inbounds begin
	for k in eachindex(u)
	    aux=u[k]/(abstol+abs(u[k])*reltol)
		norm+=aux*aux
	end
    end

	norm=sqrt(norm/(length(u)))

    return(norm)

end


function Rdigits(x::Real,r::Integer)

	aux=copy(x)
	mx=r*aux
	mxx=mx+aux
	return (mxx-mx)

end
