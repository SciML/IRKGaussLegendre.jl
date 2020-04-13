include("IRKCoefficients.jl")


mutable struct tcoeffs{T}
       mu::Array{T,2}
       hc::Array{T,1}
	   hb::Array{T,1}
	   nu::Array{T,2}
	   beta::Array{T,2}
	   beta2::Array{T,2}
end

mutable struct tcache{utype,eltypeu}
	U::Array{utype,1}
	Uz::Array{utype,1}
    L::Array{utype,1}
	Lz::Array{utype,1}
	F::Array{utype,1}
	Dmin::Array{utype,1}
	rejects::Array{Int64,1}
	nfcn::Array{Int64, 1}
	lambdas::Array{eltypeu,1}
end


abstract type IRKAlgorithm  <: OrdinaryDiffEqAlgorithm end
struct IRKGL16 <: IRKAlgorithm end

function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
     alg::IRKAlgorithm,args...;dt=0, #(prob.tspan[2]-prob.tspan[1])/10000,
     saveat=dt,
     maxiter=9,              # Float64 !!!
	 maxtrials=3,            # maximum of unsuccessful trials
     save_everystep=true,
     initial_interp=true,
	 adaptive=true,
	 reltol=1e-6,
	 abstol=1e-6,
	 myoutputs=false,
     kwargs...) where{uType,tType,isinplace}

#    println("IRKGL16....")


	s = 8
    destats = DiffEqBase.DEStats(0)

    reltol2s=sqrt(reltol)
	abstol2s=sqrt(abstol)

	@unpack f,u0,tspan,p=prob
    t0=tspan[1]
	tf=tspan[2]
	utype = typeof(u0)
	uitype = eltype(u0)
    ttype = typeof(t0)

	reltol2s=uitype(sqrt(reltol))
	abstol2s=uitype(sqrt(abstol))

    coeffs=tcoeffs{uitype}(zeros(s,s),zeros(s),zeros(s),
	                       zeros(s,s),zeros(s,s),zeros(s,s))
	@unpack mu,hc,hb,nu,beta,beta2 = coeffs


    if (dt==0)
		d0=MyNorm(u0,abstol,reltol)
		du0=similar(u0)
		f(du0, u0, p, t0)
    	d1=MyNorm(du0,abstol,reltol)
		if (d0<1e-5 || d1<1e-5)
			dt=convert(uitype,1e-6)
		else
			dt=convert(uitype,0.01*(d0/d1))
		end
		saveat=dt
#		println("dt=0 !!! dt=",dt, " typeof(dt)=",typeof(dt))
	end

	EstimateCoeffs!(beta,typeof(dt))
	EstimateCoeffs2!(beta2,typeof(dt))
	MuCoefficients!(mu,typeof(dt))

    dts = Array{typeof(dt)}(undef, 1)
    dtprev=zero(typeof(dt))

	if (adaptive==false)
		dtprev=dt
	else
		dtprev=zero(typeof(dt))
	end

	dts=[dt,dtprev]
	sdt = sign(dt)
	HCoefficients!(mu,hc,hb,nu,dt,dtprev,uitype)

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
	U6 = Array{typeof(u0)}(undef, s)
	for i in 1:s
		U1[i] = zeros(eltype(u0), size(u0))
		U2[i] = zeros(eltype(u0), size(u0))
		U3[i] = zeros(eltype(u0), size(u0))
		U4[i] = zeros(eltype(u0), size(u0))
		U5[i] = zeros(eltype(u0), size(u0))
		U6[i] = zeros(eltype(u0), size(u0))
	end

    cache=tcache{typeof(u0),eltype(u0)}(U1,U2,U3,U4,U5,U6,[0],[0],[0.,0.])
	@unpack U,Uz,L,Lz,F,Dmin,rejects,nfcn,lambdas=cache

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
		tit=0
		it=0

        for i in 1:m
		  j+=1
#         println("step:", j, " time=",tj[1]+tj[2]," dt=", dts[1], " dtprev=", dts[2])
   	      (status,it) = IRKStep!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
	                             maxtrials,initial_interp,abstol2s,reltol2s,adaptive)

         if (status=="Failure")
			 println("Fail")
			 sol=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= :Failure)
		     return(sol)
		 end
          tit+=it

        end

        cont = (sdt*(tj[1]+tj[2]) < sdt*tf) && (j<n*m)

		if (save_everystep==true) || (cont==false)
			push!(iters,convert(Int64,round(tit/m)))
			push!(uu,uj+ej)
			push!(tt,tj[1]+tj[2])
			push!(steps,dts[2])
		end

    end

#    println("End IRKGL16")

	sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= :Success)
	sol.destats.nf=nfcn[1]
	sol.destats.nreject=rejects[1]
	sol.destats.naccept=j

	if (myoutputs==true)
    	return(sol,iters,steps)
	else
		return(sol)
	end

  end


function IRKStep!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,maxtrials,
  		               initial_interp,abstol,reltol,adaptive)


   if (adaptive==true)

      if (maxiter>10)
#		  println("IRKstep_adaptive. maxiter=", maxiter)
	      (status,it)= IRKstep_adaptive!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
 		     	             maxtrials,initial_interp,abstol,reltol,adaptive)
      else
#		  println("IRKstep_adaptiveNEW maxiter=", maxiter)
		  (status,it)= IRKstep_adaptiveNEW!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
 		     	             maxtrials,initial_interp,abstol,reltol,adaptive)
	  end
   else
#	  println("IRKstep_fixed")
	  (status,it)= IRKstep_fixed!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
 			                      initial_interp,abstol,reltol,adaptive)
   end

   return(status,it)

end


function IRKstep_fixed!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiter,
		               initial_interp,abstol,reltol,adaptive)

        @unpack mu,hc,hb,nu,beta,beta2 = coeffs
        @unpack f,u0,p,tspan=prob
		@unpack U,Uz,L,Lz,F,Dmin,rejects,nfcn,lambdas=cache

		uitype = eltype(uj)

		dt=dts[1]
		dtprev=dts[2]
		tf=tspan[2]

        elems = s*length(uj)

        tj = ttj[1]
        te = ttj[2]

        nit=0

    	if (dt != dtprev)
			HCoefficients!(mu,hc,hb,nu,dt,dtprev,uitype)
			@unpack mu,hc,hb,nu,beta,beta2 = coeffs
		end

    	if initial_interp
        	for is in 1:s
            	for k in eachindex(uj)
                	aux=zero(eltype(uj))
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
			nfcn[1]+=1
        	f(F[is], U[is], p, tj + hc[is])
        	@. L[is] = hb[is]*F[is]
    	end

    	iter = true # Initialize iter outside the for loop
    	plusIt=true

    	nit=1
		for is in 1:s Dmin[is] .= Inf end

    	while (nit<maxiter && iter)

        	nit+=1
        	iter=false
        	D0=0

    		for is in 1:s
        		Uz[is] .= U[is]
        		DiffEqBase.@.. U[is] = uj + (ej+mu[is,1]*L[1] + mu[is,2]*L[2]+
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
					nfcn[1]+=1
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

		indices = eachindex(uj)

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

		res = Base.TwicePrecision(tj, te) + dt
		ttj[1] = res.hi
		ttj[2] = res.lo

    	dts[1]=min(dt,tf-(ttj[1]+ttj[2]))
		dts[2]=dt

#        println("New step size:  dt=",dts[1])
#        println("")

        return("Success",nit)


end


function IRKstep_adaptive!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiter,maxtrials,
		                   initial_interp,abstol,reltol,adaptive)


		@unpack mu,hc,hb,nu,beta,beta2 = coeffs
        @unpack f,u0,p,tspan=prob
		@unpack U,Uz,L,Lz,F,Dmin,rejects,nfcn,lambdas=cache

		uitype = eltype(uj)

		lambda=lambdas[1]
		lambdaprev=lambdas[2]

		dt=dts[1]
		dtprev=dts[2]
		tf=tspan[2]

        elems = s*length(uj)
		pow=eltype(uj)(1/(2*s))

        tj = ttj[1]
        te = ttj[2]

        accept=false
		estimate=zero(eltype(uj))

        nit=0
		ntrials=0

        if (j==1)
			mtrials=2*maxtrials
		else
			mtrials=maxtrials
		end

		for is in 1:s Lz[is].=L[is] end

		while (!accept && ntrials<mtrials)

			if (dt != dtprev)
				HCoefficients!(mu,hc,hb,nu,dt,dtprev,uitype)
				@unpack mu,hc,hb,nu,beta,beta2 = coeffs
			end

        	if initial_interp
            	for is in 1:s
                	for k in eachindex(uj)
                    	aux=zero(eltype(uj))
                    	for js in 1:s
                        	aux+=nu[is,js]*Lz[js][k]
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
				nfcn[1]+=1
            	f(F[is], U[is], p, tj + hc[is])
            	@. L[is] = hb[is]*F[is]
        	end

        	iter = true # Initialize iter outside the for loop
        	plusIt=true

        	nit=1
			for is in 1:s Dmin[is] .= Inf end

        	while (nit<maxiter && iter)

            	nit+=1
            	iter=false
            	D0=0

        		for is in 1:s
            		Uz[is] .= U[is]
            		DiffEqBase.@.. U[is] = uj + (ej+mu[is,1]*L[1] + mu[is,2]*L[2]+
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
						nfcn[1]+=1
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

            ntrials+=1

     		estimate=ErrorEst(U,F,dt,beta2,abstol,reltol)
			lambda=(estimate)^pow
			if (estimate < 2)
				accept=true
			else
			    rejects[1]+=1
				dt=dt/lambda
			end

		end # while accept


        if (!accept && ntrials==maxtrials)
			println("Fail !!!  Step=",j, " dt=", dts[1])
			return("Failure",0)
		end

		indices = eachindex(uj)

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

		res = Base.TwicePrecision(tj, te) + dt
		ttj[1] = res.hi
		ttj[2] = res.lo

		if (j==1)
            dts[1]=min(max(dt/2,min(2*dt,dt/lambda)),tf-(ttj[1]+ttj[2]))
		else
            barh=dt/lambda*(dt*lambdaprev/(dtprev*lambda))^((lambda+1)/(lambda+lambdaprev))
			dts[1]= min(max(dt/2,min(2*dt,barh)),tf-(ttj[1]+ttj[2]))
		end
		dts[2]=dt
        lambdas[2]=lambda

#        println("New step size:  dt=",dts[1])
#        println("")

        return("Success",nit)


end


function IRKstep_adaptiveNEW!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiter,maxtrials,
		                      initial_interp,abstol,reltol,adaptive)


		@unpack mu,hc,hb,nu,beta,beta2 = coeffs
        @unpack f,u0,p,tspan=prob
		@unpack U,Uz,L,Lz,F,Dmin,rejects,nfcn,lambdas=cache

		uitype = eltype(uj)

		lambda=lambdas[1]
		lambdaprev=lambdas[2]

		dt=dts[1]
		dtprev=dts[2]
		tf=tspan[2]

        elems = s*length(uj)
		pow=eltype(uj)(1/(2*s))

        tj = ttj[1]
        te = ttj[2]

        accept=false
		estimate=zero(eltype(uj))

        nit=0
		ntrials=0

        if (j==1)
			mtrials=2*maxtrials
		else
			mtrials=maxtrials
		end

		for is in 1:s Lz[is].=L[is] end

		while (!accept && ntrials<mtrials)

			if (dt != dtprev)
				HCoefficients!(mu,hc,hb,nu,dt,dtprev,uitype)
				@unpack mu,hc,hb,nu,beta,beta2 = coeffs
			end

        	for is in 1:s
            	for k in eachindex(uj)
                	aux=zero(eltype(uj))
                	for js in 1:s
                		aux+=nu[is,js]*Lz[js][k]
                	end
                	U[is][k]=(uj[k]+ej[k])+aux
            	end
        	end

        	for is in 1:s
				nfcn[1]+=1
        		f(F[is], U[is], p, tj + hc[is])
        		@. L[is] = hb[is]*F[is]
        	end

        	nit=1
			for is in 1:s Dmin[is] .= Inf end

			while (nit<maxiter)
            	nit+=1
        		for is in 1:s
            		Uz[is] .= U[is]
            		DiffEqBase.@.. U[is] = uj+(ej+mu[is,1]*L[1] + mu[is,2]*L[2]+
				                           mu[is,3]*L[3] + mu[is,4]*L[4]+
                                           mu[is,5]*L[5] + mu[is,6]*L[6]+
									       mu[is,7]*L[7] + mu[is,8]*L[8])
        		end

				for is in 1:s
					eval=false
					for k in eachindex(uj)
						if (abs(U[is][k]-Uz[is][k])>0.)
							eval=true
						end
					end
					if (eval==true)
						nfcn[1]+=1
	            		f(F[is], U[is], p, tj + hc[is])
	            		@. L[is] = hb[is]*F[is]
					end
	       	   end

        	end # while iter

            ntrials+=1

     		estimate=ErrorEst(U,F,dt,beta2,abstol,reltol)
			lambda=(estimate)^pow
			if (estimate < 2)
				accept=true
			else
			    rejects[1]+=1
				dt=dt/lambda
			end

		end # while accept


        if (!accept && ntrials==maxtrials)
			println("Fail !!!  Step=",j, " dt=", dts[1])
			return("Failure",0)
		end

		indices = eachindex(uj)

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

		res = Base.TwicePrecision(tj, te) + dt
		ttj[1] = res.hi
		ttj[2] = res.lo

		if (j==1)
            dts[1]=min(max(dt/2,min(2*dt,dt/lambda)),tf-(ttj[1]+ttj[2]))
		else
            barh=dt/lambda*(dt*lambdaprev/(dtprev*lambda))^((lambda+1)/(lambda+lambdaprev))
			dts[1]= min(max(dt/2,min(2*dt,barh)),tf-(ttj[1]+ttj[2]))
		end
		dts[2]=dt
        lambdas[2]=lambda

#        println("New step size:  dt=",dts[1])
#        println("")

        return("Success",nit)


end


function ErrorEst(U,F,dt,beta,abstol,reltol)

    (s,)=size(F)
	D=length(U[1])

	est=zero(typeof(dt))

	for k in eachindex(U[1])

		sum=zero(typeof(dt))
		maxU=zero(typeof(dt))
		for is in 1:s
        	sum+=beta[is]*F[is][k]
			maxU=max(maxU,abs(U[is][k]))
        end

		est+=(abs(dt*sum))^2/(abstol+maxU*reltol)

	end

    return(est/D)

end


function MyNorm(u,abstol,reltol)

	norm=zero(eltype(u))

	for k in eachindex(u)
	    aux=u[k]/(abstol+abs(u[k])*reltol)
		norm+=aux*aux
	end

	norm=sqrt(norm/(length(u)))

    return(norm)

end
