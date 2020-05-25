include("IRKCoefficients.jl")


mutable struct tcoeffs{T}
       mu::Array{T,2}
       hc::Array{T,1}
	   hb::Array{T,1}
	   nu::Array{T,2}
	   beta::Array{T,2}
	   beta2::Array{T,2}
end

mutable struct tcache{uType,elTypeu}
	U::Array{uType,1}
	Uz::Array{uType,1}
    L::Array{uType,1}
	Lz::Array{uType,1}
	F::Array{uType,1}
	Dmin::Array{uType,1}
	Eval::Array{Bool,1}
	rejects::Array{Int64,1}
	nfcn::Array{Int64, 1}
	lambdas::Array{elTypeu,1}
end


abstract type IRKAlgorithm  <: OrdinaryDiffEqAlgorithm end
struct IRKGL16 <: IRKAlgorithm end


function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
     alg::IRKAlgorithm,args...;
#     dt::tType=0.,
     dt=0.,
     saveat=1,
     maxiter=12,
	 maxtrials=3,            # maximum of unsuccessful trials
     save_everystep=true,
     initial_interp=true,
	 adaptive=true,
	 reltol=1e-6,
	 abstol=1e-6,
	 myoutputs=false,
     mixed_precision=false,
	 low_prec_type = Float64,
     kwargs...) where{uType,tType,isinplace}

#     println("IRKGL16....ODEproblem")


	s = 8
    destats = DiffEqBase.DEStats(0)

	@unpack f,u0,tspan,p,kwargs=prob
    t0=tspan[1]
	tf=tspan[2]
	tType2=eltype(tspan)
#	low_prec_type=typeof(kwargs[:lpp])
	uiType = eltype(u0)

    reltol2s=uiType(sqrt(reltol))
	abstol2s=uiType(sqrt(abstol))

    coeffs=tcoeffs{uiType}(zeros(s,s),zeros(s),zeros(s),
	                       zeros(s,s),zeros(s,s),zeros(s,s))
	@unpack mu,hc,hb,nu,beta,beta2 = coeffs


    if (dt==0)
		d0=MyNorm(u0,abstol,reltol)
		du0=similar(u0)
		f(du0, u0, p, t0)
    	d1=MyNorm(du0,abstol,reltol)
		if (d0<1e-5 || d1<1e-5)
			dt=convert(tType2,1e-6)
		else
			dt=convert(tType2,0.01*(d0/d1))
		end
	end

	EstimateCoeffs!(beta,uiType)
	EstimateCoeffs2!(beta2,uiType)
	MuCoefficients!(mu,uiType)

    dts = Array{tType2}(undef, 1)

	if (adaptive==false)
		dtprev=dt
	else
		dtprev=zero(tType2)
	end

	dts=[dt,dtprev]
	sdt = sign(dt)
	HCoefficients!(mu,hc,hb,nu,dt,dtprev,uiType)

#   m: output saved at every m steps
#   n: Number of macro-steps  (Output is saved for n+1 time values)


	if (adaptive==true)
	   maxiter=12
	   if (save_everystep==false)
    		m=1
		    n=Inf
       else
			m=Int64(saveat)
		    n=Inf
      end
    end

    if (adaptive==false)
		maxiter=100
	    if (save_everystep==false)
##			m=convert(Int64,ceil((tf-t0)/(dt)))
            m=1
    		n=convert(Int64,ceil((tf-t0)/(m*dt)))
    	else
##			m=convert(Int64,(round(saveat/dt)))
            m=Int64(saveat)
    		n=convert(Int64,ceil((tf-t0)/(m*dt)))
        end
    end



    U1 = Array{uType}(undef, s)
    U2 = Array{uType}(undef, s)
    U3 = Array{uType}(undef, s)
    U4 = Array{uType}(undef, s)
    U5 = Array{uType}(undef, s)
	U6 = Array{uType}(undef, s)
	for i in 1:s
		U1[i] = zero(u0)
		U2[i] = zero(u0)
		U3[i] = zero(u0)
		U4[i] = zero(u0)
		U5[i] = zero(u0)
		U6[i] = zero(u0)
	end

    cache=tcache{uType,uiType}(U1,U2,U3,U4,U5,U6,fill(true,s),[0],[0,0],[0.,0.])
	@unpack U,Uz,L,Lz,F,Dmin,Eval,rejects,nfcn,lambdas=cache

    ej=zero(u0)

	uu = Array{uType}[]
	tt = Array{tType2}[]
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
        k=0

        @inbounds begin
        for i in 1:m
		  j+=1
		  k+=1
   	      (status,it) = IRKStep!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
	                             maxtrials,initial_interp,abstol2s,reltol2s,adaptive,
								 mixed_precision,low_prec_type)

         if (status=="Failure")
			 println("Fail")
			 sol=DiffEqBase.build_solution(prob,alg,tt,uu,retcode= :Failure)
		     return(sol)
		 end
         tit+=it

		 if (dts[1]==0)
			 break
		 end

       end
	   end

	   cont = (sdt*(tj[1]+tj[2]) < sdt*tf) && (j<n*m)

       if (save_everystep==true) || (cont==false)
			push!(iters,convert(Int64,round(tit/k)))
			push!(uu,uj+ej)
			push!(tt,tj[1]+tj[2])
			push!(steps,dts[2])
       end

    end

	sol=DiffEqBase.build_solution(prob,alg,tt,uu,destats=destats,retcode= :Success)

	sol.destats.nf=nfcn[1]
	sol.destats.nf2=nfcn[2]
	sol.destats.nreject=rejects[1]
	sol.destats.naccept=j

	if (myoutputs==true)
    	return(sol,iters,steps)
	else
		return(sol)
	end

  end



function IRKStep!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,maxtrials,
  		               initial_interp,abstol,reltol,adaptive,
					   mixed_precision,low_prec_type)


   if (adaptive==true)
#	   println("IRKstep_adaptive. maxiter=", maxiter)
	  (status,it)= IRKstep_adaptive!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
 		     	             maxtrials,initial_interp,abstol,reltol,adaptive,
							 mixed_precision,low_prec_type)
   else
#	  println("IRKstep_fixed")
	  (status,it)= IRKstep_fixed!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
 			                      initial_interp,abstol,reltol,adaptive,
								  mixed_precision,low_prec_type)
   end

   return(status,it)

end


function IRKstep_fixed!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiter,
		               initial_interp,abstol,reltol,adaptive,
					   mixed_precision,low_prec_type)

        @unpack mu,hc,hb,nu,beta,beta2 = coeffs
        @unpack f,u0,p,tspan,kwargs=prob
		@unpack U,Uz,L,Lz,F,Dmin,Eval,rejects,nfcn,lambdas=cache

        lpp=kwargs[:lpp]
		uiType = eltype(uj)

		dt=dts[1]
		dtprev=dts[2]
		tf=tspan[2]

        elems = s*length(uj)

        tj = ttj[1]
        te = ttj[2]

        nit=0

    	if (dt != dtprev)
			HCoefficients!(mu,hc,hb,nu,dt,dtprev,uiType)
			@unpack mu,hc,hb,nu,beta,beta2 = coeffs
		end

    	if initial_interp
			@inbounds begin
        	for is in 1:s
            	for k in eachindex(uj)
                	aux=zero(eltype(uj))
                	for js in 1:s
                    	aux+=nu[is,js]*L[js][k]
                	end
                	U[is][k]=(uj[k]+ej[k])+aux
            	end
        	end
		    end
    	else
			@inbounds begin
        	for is in 1:s
            	@. U[is] = uj + ej
        	end
			end
    	end


    	iter = true # Initialize iter outside the for loop
    	plusIt=true

    	nit=1
		@inbounds begin
		for is in 1:s Dmin[is] .= Inf end
	    end


##
##  Low-prec iterations
##

##        println("***************************************************")
##        println("urratsa=",j)

        if (mixed_precision==true)

#        setprecision(BigFloat,precision(low_prec_type)+11)   #+11?

		@inbounds begin
		Threads.@threads for is in 1:s
				nfcn[2]+=1
				f(F[is], convert.(low_prec_type,U[is]),lpp,
				convert(low_prec_type, tj + hc[is]))
				@. L[is] = hb[is]*F[is]
		end
		end


    	while (iter)

            nit+=1

        	iter=false
        	D0=0

            @inbounds begin
    		for is in 1:s
        		Uz[is] .= U[is]
        		DiffEqBase.@.. U[is] = uj +
				                      (ej+mu[is,1]*L[1] + mu[is,2]*L[2]+
			                           mu[is,3]*L[3] + mu[is,4]*L[4]+
                                       mu[is,5]*L[5] + mu[is,6]*L[6]+
								       mu[is,7]*L[7] + mu[is,8]*L[8])
    		end

        	Threads.@threads for is in 1:s
            	Eval[is]=false
            	for k in eachindex(uj)
                    DY=abs(U[is][k]-Uz[is][k])
                    if DY>0.
                       Eval[is]=true
                       if DY< Dmin[is][k]
                          Dmin[is][k]=DY
                          iter=true
                       end
                   else
                       D0+=1
                   end
            	end

           		if Eval[is]==true
					nfcn[2]+=1
   				    f(F[is], convert.(low_prec_type,U[is]),lpp,
				     		 convert(low_prec_type, tj + hc[is]))
					@. L[is] = hb[is]*F[is]
           		end

    		end
		    end

#			println("low-prec-iters=",nit-1, " D0=", D0," Norm=", norm(U.-Uz))

    	end # while iter low-prec

#       println(" ")

#       setprecision(BigFloat,108)
#	   println("Recover precision:",precision(BigFloat))

#	   println(" ")

       iter = true # Initialize iter outside the for loop
       plusIt=true

#	   @inbounds begin
#	   for is in 1:s Dmin[is] .= Inf end
#	   end



	end # if (mixed_precision==true)

	   ##
	   ##  High-prec iterations
	   ##

	   @inbounds begin
	   Threads.@threads for is in 1:s
		   nfcn[1]+=1
	       f(F[is], U[is], p, tj + hc[is])
		   @. L[is] = hb[is]*F[is]
	   end
       end


       while (nit<maxiter && iter)

      	  nit+=1

		  iter=false
	      D0=0

	      @inbounds begin
	      for is in 1:s
	    	Uz[is] .= U[is]
		    DiffEqBase.@.. U[is] = uj + (ej+mu[is,1]*L[1] + mu[is,2]*L[2]+
			  				   mu[is,3]*L[3] + mu[is,4]*L[4]+
							   mu[is,5]*L[5] + mu[is,6]*L[6]+
							   mu[is,7]*L[7] + mu[is,8]*L[8])
	     end

	     Threads.@threads for is in 1:s
	     	Eval[is]=false
		    for k in eachindex(uj)
			    DY=abs(U[is][k]-Uz[is][k])
			    if DY>0.
			       Eval[is]=true
			       if DY< Dmin[is][k]
				      Dmin[is][k]=DY
				      iter=true
			      end
		       else
			      D0+=1
		      end
		   end

		   if Eval[is]==true
	    		nfcn[1]+=1
		    	f(F[is], U[is],lpp, tj + hc[is])
			    @. L[is] = hb[is]*F[is]
		   end
	    end
	    end

	    if (iter==false && D0<elems && plusIt)
		    iter=true
	    	plusIt=false
	    else
		    plusIt=true
	    end

#		 println("high-prec-iters=",nit-1, " D0=", D0, " Norm=", norm(U.-Uz))

     end # while iter high-prec



		indices = eachindex(uj)

        @inbounds begin
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
	    end

		res = Base.TwicePrecision(tj, te) + dt
		ttj[1] = res.hi
		ttj[2] = res.lo

    	dts[1]=min(dt,tf-(ttj[1]+ttj[2]))
		dts[2]=dt

        return("Success",nit)


end



function IRKstep_adaptive!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiter,maxtrials,
		                      initial_interp,abstol,reltol,adaptive,
							  mixed_precision,low_prec_type)


		@unpack mu,hc,hb,nu,beta,beta2 = coeffs
        @unpack f,u0,p,tspan=prob
		@unpack U,Uz,L,Lz,F,Dmin,Eval,rejects,nfcn,lambdas=cache

		uiType = eltype(uj)

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
			maxtrialsj=4*maxtrials
			maxiterj=2*maxiter
		else
			maxtrialsj=maxtrials
			maxiterj=maxiter
		end

		for is in 1:s Lz[is].=L[is] end

		while (!accept && ntrials<maxtrialsj)

			if (dt != dtprev)
				HCoefficients!(mu,hc,hb,nu,dt,dtprev,uiType)
				@unpack mu,hc,hb,nu,beta,beta2 = coeffs
			end

            @inbounds begin
        	for is in 1:s
            	for k in eachindex(uj)
                	aux=zero(eltype(uj))
                	for js in 1:s
                		aux+=nu[is,js]*Lz[js][k]
                	end
                	U[is][k]=(uj[k]+ej[k])+aux
            	end
        	end

         	Threads.@threads for is in 1:s
				nfcn[1]+=1
#                atomic_add!(nfcn[1],1)
        		f(F[is], U[is], p, tj + hc[is])
        		@. L[is] = hb[is]*F[is]
        	end
     		end

        	nit=1
			for is in 1:s Dmin[is] .= Inf end

			while (nit<maxiterj)
            	nit+=1

				@inbounds begin
        		for is in 1:s
            		Uz[is] .= U[is]
            		DiffEqBase.@.. U[is] = uj+(ej+mu[is,1]*L[1] + mu[is,2]*L[2]+
				                           mu[is,3]*L[3] + mu[is,4]*L[4]+
                                           mu[is,5]*L[5] + mu[is,6]*L[6]+
									       mu[is,7]*L[7] + mu[is,8]*L[8])
        		end

				Threads.@threads for is in 1:s
                    Eval[is]=false
					for k in eachindex(uj)
						if (abs(U[is][k]-Uz[is][k])>0.)
                            Eval[is]=true
						end
					end
					if (Eval[is]==true)
						nfcn[1]+=1
#		                atomic_add!(nfcn[1],1)
	            		f(F[is], U[is], p, tj + hc[is])
	            		@. L[is] = hb[is]*F[is]
					end
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

        @inbounds begin
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
	    end

		res = Base.TwicePrecision(tj, te) + dt
		ttj[1] = res.hi
		ttj[2] = res.lo

		if (j==1)
            dts[1]=min(max(dt/2,min(2*dt,dt/lambda)),tf-(ttj[1]+ttj[2]))
		else
            hath1=dt/lambda
			hath2=dtprev/lambdaprev
			tildeh=hath1*(hath1/hath2)^(lambda/lambdaprev)
			barlamb1=(dt+tildeh)/(hath1+tildeh)
			barlamb2=(dtprev+dt)/(hath2+hath1)
			barh=hath1*(hath1/hath2)^(barlamb1/barlamb2)
			dts[1]= min(max(dt/2,min(2*dt,barh)),tf-(ttj[1]+ttj[2]))
		end
		dts[2]=dt
        lambdas[2]=lambda

        return("Success",nit)


end


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
