include("IRKCoefficients.jl")

#abstract type OrdinaryDiffEqAdaptiveAlgorithm <: OrdinaryDiffEqAlgorithm end
#struct IRK8 <: OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm end
# solve ( alg::IRK8;)
#

mutable struct tcoeffs{T}
       mu::Array{T,2}
       hc::Array{T,1}
	   hb::Array{T,1}
	   nu::Array{T,2}
	   beta::Array{T,2}
end

mutable struct tcache{utype}
	U::Array{utype,1}
	Uz::Array{utype,1}
    L::Array{utype,1}
	F::Array{utype,1}
	Dmin::Array{utype,1}
    uz::utype
	ez::utype
end


abstract type IRKAlgorithm  <: OrdinaryDiffEqAlgorithm end
struct IRK8 <: IRKAlgorithm end

function DiffEqBase.solve(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
     alg::IRKAlgorithm;dt=(prob.tspan[2]-prob.tspan[1])/100,
     saveat=dt,
     maxiter=100,
     save_everystep=true,
     initial_interp=true,
	 adaptive=true,
	 reltol=1e-3,
	 abstol=1e-6,
	 myoutputs=false,
     kwargs...) where{uType,tType,isinplace}

#    println("IRK8....")

	s = 8

    reltol2s=sqrt(reltol)
	abstol2s=sqrt(abstol)

    coeffs=tcoeffs{typeof(dt)}(zeros(s,s),zeros(s),zeros(s),
	                           zeros(s,s),zeros(s,s))
	@unpack mu,hc,hb,nu,beta = coeffs

	EstimateCoeffs!(beta,typeof(dt))
	MuCoefficients!(mu,typeof(dt))

	if (adaptive==false)
		HCoefficients!(mu,hc,hb,nu,dt,dt)
	else
		HCoefficients!(mu,hc,hb,nu,dt,0.)
	end

    dtprev=0.
	rejects=0
	nfcn=0
	dts=[dt,dtprev,rejects,nfcn]
    sdt = sign(dt)

    @unpack f,u0,tspan,p=prob
    t0=tspan[1]
	tf=tspan[2]
	utype = typeof(u0)
    ttype = typeof(t0)

#   m: output saved at every m steps
#   n: Number of macro-steps  (Output is saved for n+1 time values)
	if (adaptive==true)
		m=1
		n=Inf
    elseif (save_everystep==false)
			m=convert(Int64,ceil((tf-t0)/(dt)))
			n=1
    	else
			m=convert(Int64,(round(saveat/dt)))
    		n=convert(Int64,ceil((tf-t0)/(m*dt)))
	end

    U1 = Array{typeof(u0)}(undef, s)
	U2 = Array{typeof(u0)}(undef, s)
    U3 = Array{typeof(u0)}(undef, s)
    U4 = Array{typeof(u0)}(undef, s)
    U5 = Array{typeof(u0)}(undef, s)

	for i in 1:s
		U1[i] = zeros(eltype(u0), size(u0))
		U2[i] = zeros(eltype(u0), size(u0))
		U3[i] = zeros(eltype(u0), size(u0))
		U4[i] = zeros(eltype(u0), size(u0))
		U5[i] = zeros(eltype(u0), size(u0))
	end

    cache=tcache{typeof(u0)}(U1,U2,U3,U4,U5,zeros(size(u0)),zeros(size(u0)))

	@unpack U,Uz,L,F,Dmin,uz,ez=cache

    ej=zeros(eltype(u0), size(u0))

	uu = Array{typeof(u0)}[]
	tt = Array{ttype}[]
    iters = Array{Int}[]
	steps = Array{Int}[]

	uu=[]
	tt=[]
	iters=[]
	steps=[]

	push!(uu,u0)
	push!(tt,t0)
	push!(iters,0)
	push!(steps,0)
    tj = [t0, zero(t0)]
    uj = copy(u0)

    j=0
    cont=true

    while cont
        j+=1
		tit=0
		it=0

        for i in 1:m
#         println("step:", j, " time=",tj[1]+tj[2]," dt=", dts[1], " dtprev=", dts[2])

		 (it) = IRKstep!(s,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
				 	     initial_interp,abstol2s,reltol2s,adaptive)
          tit+=it
        end

        cont = (sdt*(tj[1]+tj[2]) < sdt*tf) && (j<n)

		if (save_everystep==true) || (cont==false)
			push!(iters,convert(Int64,round(tit/m)))
			push!(uu,uj+ej)
			push!(tt,tj[1]+tj[2])
			push!(steps,dts[2])
		end

    end

#    println("End IRK8")
	sol=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= :Success)
	if (myoutputs==true)
		rejects=dts[3]
		nfcn=dts[4]
    	return(sol,iters,steps,rejects,nfcn)
	else
		return(sol)
	end

  end


function IRKstep!(s,ttj,uj,ej,prob,dts,coeffs,cache,maxiter,
		           initial_interp,abstol,reltol,adaptive)


		@unpack mu,hc,hb,nu,beta = coeffs
        @unpack f,u0,p,tspan=prob
		@unpack U,Uz,L,F,Dmin,uz,ez=cache

		dt=dts[1]
		dtprev=dts[2]
		rejects=dts[3]
		nfcn=dts[4]
		tf=tspan[2]

        elems = s*length(uj)
		pow=1/(s-1)

        tj = ttj[1]
        te = ttj[2]

        accept=false
		estimate=zero

		j=0   # outside while-accept

		while (!accept)

			if (adaptive == true)
				HCoefficients!(mu,hc,hb,nu,dt,dtprev)
				@unpack mu,hc,hb,nu,beta = coeffs
			end

        	if initial_interp
            	for is in 1:s
                	for k in eachindex(uj)
                    	aux=0.
                    	for js in 1:s
                        	aux+=nu[is,js]*L[js][k]
                    	end
                    	U[is][k]=(uj[k]+ej[k])+aux
                	end
            	end
        	else
            	for is in 1:s
                	@. U[is] = uj + ej
            	end
        	end

        	for is in 1:s
				nfcn+=1
            	f(F[is], U[is], p, tj + hc[is])
            	@. L[is] = hb[is]*F[is]
        	end

        	iter = true # Initialize iter outside the for loop
        	plusIt=true

        	j=1
			for is in 1:s Dmin[is] .= Inf end

        	while (j<maxiter && iter)

            	j+=1
            	iter=false
            	D0=0

        		for is in 1:s
            		Uz[is] .= U[is]
            		@. U[is] = uj + (ej+mu[is,1]*L[1] + mu[is,2]*L[2]+
				                    mu[is,3]*L[3] + mu[is,4]*L[4]+
                                 	mu[is,5]*L[5] + mu[is,6]*L[6]+
									mu[is,7]*L[7] + mu[is,8]*L[8])
        		end

            	for is in 1:s
                	eval=false
                	for k in eachindex(uj)
                        DY=abs(U[is][k]-Uz[is][k])
                        if DY>0.
                           eval=true
                           if DY< Dmin[is][k]
                              Dmin[is][k]=DY
                              iter=true
                           end
                       else
                           D0+=1
                       end
                	end

               		if eval==true
						nfcn+=1
                  		f(F[is], U[is], p,  tj + hc[is])
                  		@. L[is] = hb[is]*F[is]
               		end
           		end

            	if (iter==false && D0<elems && plusIt)
                	iter=true
                	plusIt=false
            	else
                	plusIt=true
            	end

        	end # while iter

        	iter = j
        	indices = eachindex(uj)

            uz.=uj
			ez.=ej

        	for k in indices    #Compensated summation
            	e0 = ej[k]
            	for is in 1:s
                	e0 += muladd(F[is][k], hb[is], -L[is][k])
            	end
            	res = Base.TwicePrecision(uj[k], e0)
            	for is in 1:s
	            	res += L[is][k]
            	end
            	uj[k] = res.hi
            	ej[k] = res.lo
         	end

            if (adaptive==false)
				accept=true
			else
				estimate=ErrorEstMax(uj+ej,uz+ez,beta,F,abstol,reltol)
				if (estimate < 2)
					accept=true
				else
#				println("Rejected step. dt=",dt, " estmate=", estimate)
				    rejects+=1
					uj.=uz
					ej.=ez
					dt=dt*(1/(estimate)^pow)
				end
			end

		end # while accept

		res = Base.TwicePrecision(tj, te) + dt
		ttj[1] = res.hi
		ttj[2] = res.lo

        if (adaptive==true)
        	dts[2]=dt
	    	dts[1]=min(dt*max(0.5, min(2,(1/estimate)^pow)),tf-(ttj[1]+ttj[2]))
			dts[3]=rejects
		end

       dts[4]=nfcn
#        println("New step size:  dt=",dts[1])
#        println("")

        return  (j)


end



function ErrorEstMax(uj,uz,beta,F,abstol,reltol)

    (s,)=size(F)
	est=0.

	for k in eachindex(uj)

		sum=0.
		for is in 1:s
        	sum+=beta[is]*F[is][k]
        end

		est=max(est,abs(sum)/(abstol+max(abs(uj[k]),abs(uz[k]))*reltol))

	end

    return(est)

end



function ErrorEst(uj,uz,abstol,reltol)

	est=0.

	for k in eachindex(uj)
	    aux=(uj[k]-uz[k])/(abstol+max(abs(uj[k]),abs(uz[k]))*reltol)
		est+=aux*aux
	end

	est=sqrt(est/(length(uj)))

    return(est)

end
