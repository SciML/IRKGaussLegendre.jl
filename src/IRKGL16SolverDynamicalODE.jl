#include("IRKCoefficients.jl")

abstract type IRKAlgorithm3  <: OrdinaryDiffEqAlgorithm end
struct IRKGL163 <: IRKAlgorithm3 end

#function DiffEqBase.__solve(prob::DiffEqBase.DynamicalODEProblem{isinplace},
function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
     alg::IRKAlgorithm3,args...;
#     dt::tType=0.,
     dt=0.,
     saveat=1,
     maxiter=100,              #10
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

#    println("IRKGL163....DynamicalProblem")


	s = 8
    destats = DiffEqBase.DEStats(0)

	@unpack tspan,p=prob
	f1=prob.f.f1
	f2=prob.f.f2
	r0=prob.u0.x[1]
	v0=prob.u0.x[2]

	u0 = ArrayPartition(r0,v0)
    t0=tspan[1]
	tf=tspan[2]
	tType2=eltype(tspan)
	uiType = eltype(u0)

	uLowType=typeof(convert.(low_prec_type,u0))

    reltol2s=uiType(sqrt(reltol))
	abstol2s=uiType(sqrt(abstol))

    coeffs=tcoeffs{uiType}(zeros(s,s),zeros(s),zeros(s),
	                       zeros(s,s),zeros(s,s))
	@unpack mu,hc,hb,nu,alpha = coeffs


    if (dt==0)
		d0=MyNorm(u0,abstol2s,reltol2s)
		du0=similar(u0)
		f1(du0.x[1], u0.x[1],u0.x[2], p, t0)
		f2(du0.x[2], u0.x[1],u0.x[2], p, t0)
    	d1=MyNorm(du0,abstol2s,reltol2s)
		if (d0<1e-5 || d1<1e-5)
			dt=convert(tType2,1e-6)
		else
			dt=convert(tType2,0.01*(d0/d1))
		end
	end

	dt=min(dt,tf-t0)

	EstimateCoeffs!(alpha,uiType)
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
	   if (save_everystep==false)
    		m=1
		    n=Inf
       else
			m=Int64(saveat)
		    n=Inf
      end
    end

    if (adaptive==false)
	    if (save_everystep==false)
#			m=convert(Int64,ceil((tf-t0)/(dt)))
            m=1
#    		n=convert(Int64,ceil((tf-t0)/(m*dt)))
			n=Inf
    	else
#			m=convert(Int64,(round(saveat/dt)))
            m=Int64(saveat)
#    		n=convert(Int64,ceil((tf-t0)/(m*dt)))
			n=Inf
        end
    end

	U1 = Array{uType}(undef, s)
    U2 = Array{uType}(undef, s)
    U3 = Array{uType}(undef, s)
    U4 = Array{uType}(undef, s)
    U5 = Array{uType}(undef, s)
	U6 = Array{uType}(undef, s)
	U11 = Array{uLowType}(undef, s)
	U12 = Array{uLowType}(undef, s)
	U13 = Array{uLowType}(undef, s)
	U14 = Array{uLowType}(undef, s)
	U15 = Array{uLowType}(undef, s)
	U16 = Array{uLowType}(undef, s)
	U17 = Array{uLowType}(undef, s)
	for i in 1:s
		U1[i] = zero(u0)
		U2[i] = zero(u0)
		U3[i] = zero(u0)
		U4[i] = zero(u0)
		U5[i] = zero(u0)
		U6[i] = zero(u0)
		U11[i] = zero(convert.(low_prec_type,u0))
		U12[i] = zero(convert.(low_prec_type,u0))
		U13[i] = zero(convert.(low_prec_type,u0))
		U14[i] = zero(convert.(low_prec_type,u0))
		U15[i] = zero(convert.(low_prec_type,u0))
		U16[i] = zero(convert.(low_prec_type,u0))
		U17[i] = zero(convert.(low_prec_type,u0))
	end

	lmu=convert.(low_prec_type,mu)
	lhb=convert.(low_prec_type,hb)

    cache=tcache{uType,uiType,uLowType,low_prec_type}(U1,U2,U3,U4,U5,U6,
	             fill(true,s),[0],[0,0],fill(zero(uiType),2),
	             U11,U12,U13,U14,U15,U16,U17,
				 fill(zero(low_prec_type),s),lhb,lmu)
	@unpack U,Uz,L,Lz,F,Dmin,Eval,rejects,nfcn,lambdas,
	        Ulow,DU,DF,DL,delta,Fa,Fb,normU,lhb,lmu=cache

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
   	      (status,it) = IRKStepDynODE!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
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



function IRKStepDynODE!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,maxtrials,
  		               initial_interp,abstol,reltol,adaptive,
					   mixed_precision,low_prec_type)


   if (adaptive==true)
#	   println("IRKstep_adaptive. maxiter=", maxiter)
	  (status,it)= IRKstepDynODE_adaptive!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
 		     	             maxtrials,initial_interp,abstol,reltol,adaptive,
							 mixed_precision,low_prec_type)
   else
#	  println("IRKstep_fixed")
     (status,it)= IRKstepDynODE_fixed!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiter,
 									 initial_interp,abstol,reltol,adaptive)
   end

   return(status,it)

end


function IRKstepDynODE_fixed!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiter,
		               initial_interp,abstol,reltol,adaptive)

		@unpack mu,hc,hb,nu,alpha = coeffs
		@unpack tspan,p=prob
		f1=prob.f.f1
		f2=prob.f.f2
#		r0=prob.v0
#		v0=prob.u0
		@unpack U,Uz,L,Lz,F,Dmin,Eval,rejects,nfcn,lambdas=cache

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
			@unpack mu,hc,hb,nu,alpha = coeffs
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
		for is in 1:s Dmin[is] .= Inf end

        @inbounds begin
    	Threads.@threads for is in 1:s
			nfcn[1]+=1
        	f1(F[is].x[1], U[is].x[1],U[is].x[2], p, tj + hc[is])
			f2(F[is].x[2], U[is].x[1],U[is].x[2], p, tj + hc[is])
        	@. L[is] = hb[is]*F[is]
    	end
	    end

#		println("***************************************************")
#       println("urratsa=",j)

    	while (nit<maxiter && iter)

        	nit+=1
        	iter=false
        	D0=0

#           First part

            @inbounds begin
    		for is in 1:s
        		Uz[is].x[1] .= U[is].x[1]
        		DiffEqBase.@.. U[is].x[1] = uj.x[1] + (ej.x[1]+mu[is,1]*L[1].x[1] + mu[is,2]*L[2].x[1]+
			                           mu[is,3]*L[3].x[1] + mu[is,4]*L[4].x[1]+
                                       mu[is,5]*L[5].x[1] + mu[is,6]*L[6].x[1]+
								       mu[is,7]*L[7].x[1] + mu[is,8]*L[8].x[1])
    		end

        	Threads.@threads for is in 1:s
            	Eval[is]=false
            	for k in eachindex(uj.x[1])
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
              		f2(F[is].x[2], U[is].x[1],U[is].x[2], p,  tj + hc[is])
              		@. L[is].x[2] = hb[is]*F[is].x[2]
           		end

    		end


#           Second part

			for is in 1:s
        		Uz[is].x[2] .= U[is].x[2]
        		DiffEqBase.@.. U[is].x[2] = uj.x[2] + (ej.x[2]+mu[is,1]*L[1].x[2] + mu[is,2]*L[2].x[2]+
	                           mu[is,3]*L[3].x[2] + mu[is,4]*L[4].x[2]+
                               mu[is,5]*L[5].x[2] + mu[is,6]*L[6].x[2]+
    					       mu[is,7]*L[7].x[2] + mu[is,8]*L[8].x[2])
    		end

        	Threads.@threads for is in 1:s
            	Eval[is]=false
        		for k in eachindex(uj.x[2])
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
#					nfcn[1]+=1
              		f1(F[is].x[1], U[is].x[1], U[is].x[2], p,  tj + hc[is])
              		@. L[is].x[1] = hb[is]*F[is].x[1]
           		end

    		end
	    	end #inbounds


        	if (iter==false && D0<elems && plusIt)
            	iter=true
            	plusIt=false
        	else
            	plusIt=true
        	end

#	  	 println(j,",","high-=",nit-1, ",", D0, ",", norm(U.-Uz))

    	end # while iter


	    if (uiType<:CompiledFloats)

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

		else
		 	@. uj+=L[1]+L[2]+L[3]+L[4]+L[5]+L[6]+L[7]+L[8]
		 	ttj[1]=tj+dt
		end

    	dts[1]=min(dt,tf-(ttj[1]+ttj[2]))
		dts[2]=dt

        return("Success",nit)


end



function IRKstepDynODE_adaptive!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiter,maxtrials,
		                      initial_interp,abstol,reltol,adaptive,
							  mixed_precision,low_prec_type)


		@unpack mu,hc,hb,nu,alpha = coeffs
		@unpack tspan,p=prob
		f1=prob.f.f1
		f2=prob.f.f2
#		r0=prob.v0
#		v0=prob.u0
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
		else
			maxtrialsj=maxtrials
		end

		for is in 1:s Lz[is].=L[is] end

		while (!accept && ntrials<maxtrialsj)

			if (dt != dtprev)
				HCoefficients!(mu,hc,hb,nu,dt,dtprev,uiType)
				@unpack mu,hc,hb,nu,alpha = coeffs
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

			iter = true # Initialize iter outside the for loop
			plusIt=true
        	nit=1
			for is in 1:s Dmin[is] .= Inf end

        	Threads.@threads for is in 1:s
				nfcn[1]+=1
				f1(F[is].x[1], U[is].x[1],U[is].x[2], p, tj + hc[is])
				f2(F[is].x[2], U[is].x[1],U[is].x[2], p, tj + hc[is])
        		@. L[is] = hb[is]*F[is]
        	end
		    end

			while (nit<maxiter && iter)

            	nit+=1
                iter=false
				D0=0
#               First part

                @inbounds begin
        		for is in 1:s
            		Uz[is].x[1] .= U[is].x[1]
            		DiffEqBase.@.. U[is].x[1] = uj.x[1]+(ej.x[1]+mu[is,1]*L[1].x[1] + mu[is,2]*L[2].x[1]+
				                           mu[is,3]*L[3].x[1] + mu[is,4]*L[4].x[1]+
                                           mu[is,5]*L[5].x[1] + mu[is,6]*L[6].x[1]+
									       mu[is,7]*L[7].x[1] + mu[is,8]*L[8].x[1])
        		end

				Threads.@threads for is in 1:s
	                Eval[is]=false
					for k in eachindex(uj.x[1])
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

                   if (Eval[is]==true)
						nfcn[1]+=1
	            		f2(F[is].x[2], U[is].x[1],U[is].x[2], p,  tj + hc[is])
	            		@. L[is].x[2] = hb[is]*F[is].x[2]
				   end
	       	    end
		       end #inbound

#               Second part

                @inbounds begin
        		for is in 1:s
            		Uz[is].x[2] .= U[is].x[2]
            		DiffEqBase.@.. U[is].x[2] = uj.x[2]+(ej.x[2]+mu[is,1]*L[1].x[2] + mu[is,2]*L[2].x[2]+
				                           mu[is,3]*L[3].x[2] + mu[is,4]*L[4].x[2]+
                                           mu[is,5]*L[5].x[2] + mu[is,6]*L[6].x[2]+
									       mu[is,7]*L[7].x[2] + mu[is,8]*L[8].x[2])
        		end

				Threads.@threads for is in 1:s
					Eval[is]=false
					for k in eachindex(uj.x[2])
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

                   if (Eval[is]==true)
#						nfcn[1]+=1
	            		f1(F[is].x[1], U[is].x[1],U[is].x[2], p,  tj + hc[is])
	            		@. L[is].x[1] = hb[is]*F[is].x[1]
				   end
	       	   end
		      end #inbound

			  if (iter==false && D0<elems && plusIt)
					iter=true
					plusIt=false
			  else
					plusIt=true
			  end


        	end # while iter

            ntrials+=1

     		estimate=ErrorEst(U,F,dt,alpha,abstol,reltol)
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


		if (uiType<:CompiledFloats)

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

	    else
		 	@. uj+=L[1]+L[2]+L[3]+L[4]+L[5]+L[6]+L[7]+L[8]
		 	ttj[1]=tj+dt
		end


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
