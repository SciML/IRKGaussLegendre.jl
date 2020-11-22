
mutable struct tcoeffs{T}
       mu::Array{T,2}
       hc::Array{T,1}
	   hb::Array{T,1}
	   nu::Array{T,2}
	   alpha::Array{T,2}
end

mutable struct tcache{uType,elTypeu,uLowType,low_prec_type}
	U::Array{uType,1}
	Uz::Array{uType,1}
	L::Array{uType,1}
	Lz::Array{uType,1}
	F::Array{uType,1}
	Dmin::Array{uType,1}
	Eval::Array{Bool,1}
	DY::Array{elTypeu,1}
	rejects::Array{Int64,1}
	nfcn::Array{Int64, 1}
	lambdas::Array{elTypeu,1}
end


mutable struct tcacheMix{uType,elTypeu,uLowType,low_prec_type}
	U::Array{uType,1}
	Uz::Array{uType,1}
	L::Array{uType,1}
	Lz::Array{uType,1}
	F::Array{uType,1}
	Dmin::Array{uType,1}
	Eval::Array{Bool,1}
	DY::Array{elTypeu,1}
	rejects::Array{Int64,1}
	nfcn::Array{Int64, 1}
	lambdas::Array{elTypeu,1}
	Ulow:: Array{uLowType,1}
	DU::Array{uLowType,1}
	DF::Array{uLowType,1}
	DL::Array{uLowType,1}
	delta::Array{uLowType,1}
	Fa::Array{uLowType,1}
	Fb::Array{uLowType,1}
	normU::Array{low_prec_type,1}
	lhb::Array{low_prec_type,1}
	lmu::Array{low_prec_type,2}
end



abstract type IRKAlgorithm  <: OrdinaryDiffEqAlgorithm end
struct IRKGL16 <: IRKAlgorithm end


function DiffEqBase.__solve(prob::DiffEqBase.AbstractODEProblem{uType,tType,isinplace},
     alg::IRKAlgorithm,args...;
#     dt::tType=0.,
     dt=0.,
     saveat=1,
     maxiters=100,
	 maxtrials=3,            # maximum of unsuccessful trials
     save_everystep=true,
     initial_interp=true,
	 adaptive=true,
	 reltol=1e-6,
	 abstol=1e-6,
	 myoutputs=false,
	 threading=false,
     mixed_precision=false,
	 low_prec_type = Float64,
     kwargs...) where{uType,tType,isinplace}

	 s = 8
     destats = DiffEqBase.DEStats(0)

    if (typeof(prob.f)<:DynamicalODEFunction)
	    @unpack tspan,p=prob
        f1=prob.f.f1
        f2=prob.f.f2
        r0=prob.u0.x[1]
        v0=prob.u0.x[2]
	    u0 = ArrayPartition(r0,v0)
    elseif (typeof(prob.f)<:ODEFunction)
	    @unpack f,u0,tspan,p,kwargs=prob
    else
     println("Error: incorrect ODEFunction")
	 sol=DiffEqBase.build_solution(prob,alg,[],[],retcode= :Failure)
     return(sol)
    end

    t0=tspan[1]
    tf=tspan[2]
    tType2=eltype(tspan)
    uiType = eltype(u0)

    lu0 = convert.(low_prec_type,u0)
	uLowType=typeof(lu0)

	coeffs=tcoeffs{uiType}(zeros(s,s),zeros(s),zeros(s),
	               zeros(s,s),zeros(s,s))

	@unpack mu,hc,hb,nu,alpha = coeffs

    Treltol=convert(uiType,reltol)
    Tabstol=convert(uiType,abstol)

    if (dt==0)
		d0=MyNorm(u0,Tabstol,Treltol)
		du0=similar(u0)
		if (typeof(prob.f)<:DynamicalODEFunction)
			f1(du0.x[1], u0.x[1],u0.x[2], p, t0)
            f2(du0.x[2], u0.x[1],u0.x[2], p, t0)
        else # ODEFunction
		    f(du0, u0, p, t0)
		end
    	d1=MyNorm(du0,Tabstol,Treltol)
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

   if (save_everystep==false)
		m=1
	    n=Inf
   else
		m=Int64(saveat)
	    n=Inf
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

    if (mixed_precision==true && typeof(prob.f)<:ODEFunction)

		U11 = Array{uLowType}(undef, s)
		U12 = Array{uLowType}(undef, s)
		U13 = Array{uLowType}(undef, s)
		U14 = Array{uLowType}(undef, s)
		U15 = Array{uLowType}(undef, s)
		U16 = Array{uLowType}(undef, s)
		U17 = Array{uLowType}(undef, s)
		for i in 1:s
			U11[i] = zero(lu0)
			U12[i] = zero(lu0)
			U13[i] = zero(lu0)
			U14[i] = zero(lu0)
			U15[i] = zero(lu0)
			U16[i] = zero(lu0)
			U17[i] = zero(lu0)
		end

		lmu=convert.(low_prec_type,mu)
		lhb=convert.(low_prec_type,hb)

		cache=tcacheMix{uType,uiType,uLowType,low_prec_type}(U1,U2,U3,U4,U5,U6,
		             fill(true,s),fill(zero(uiType),s),
					 [0],[0,0],fill(zero(uiType),2),
		             U11,U12,U13,U14,U15,U16,U17,
					 fill(zero(low_prec_type),s),lhb,lmu)

		@unpack U,Uz,L,Lz,F,Dmin,Eval,DY,rejects,nfcn,lambdas,
		        Ulow,DU,DF,DL,delta,Fa,Fb,normU,lhb,lmu=cache

	 else
		cache=tcache{uType,uiType,uLowType,low_prec_type}(U1,U2,U3,U4,U5,U6,
	 				 fill(true,s),fill(zero(uiType),s),
					 [0],[0,0],fill(zero(uiType),2))
		@unpack U,Uz,L,Lz,F,Dmin,Eval,DY,rejects,nfcn,lambdas=cache

    end

#   initialization output variables
	uu = uType[]
	tt = tType2[]
    iters = Int[]
	steps = Int[]

	push!(uu,u0)
	push!(tt,t0)
	push!(iters,0)
	push!(steps,0)
    tj = [t0, zero(t0)]
    uj = copy(u0)
	ej=zero(u0)

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
   	      (status,it) = IRKStep!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiters,
	                             maxtrials,initial_interp,Tabstol,Treltol,adaptive,
								 threading,mixed_precision,low_prec_type)

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
		    push!(tt,tj[1]+tj[2])
			push!(uu,uj+ej)

			if (myoutputs==true)
			    push!(iters,convert(Int64,round(tit/k)))
			    push!(steps,dts[2])
			end

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



function IRKStep!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiters,maxtrials,
  		               initial_interp,abstol,reltol,adaptive,threading,
					   mixed_precision,low_prec_type)

   if (typeof(prob.f)<:ODEFunction)

   		if (adaptive==true)
      		if (mixed_precision==true)
	    		(status,it)= IRKstep_adaptive_Mix!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiters,
 		     	             maxtrials,initial_interp,abstol,reltol,adaptive,
						     threading, mixed_precision,low_prec_type)
	  		else
	    		(status,it)= IRKstep_adaptive!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiters,
   		     	             maxtrials,initial_interp,abstol,reltol,adaptive,
							 threading)
      		end
       else
      	   if (mixed_precision==true)
	            (status,it)= IRKstep_fixed_Mix!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiters,
 			                 initial_interp,abstol,reltol,adaptive,threading,
				   			 mixed_precision,low_prec_type)
	       else
		        (status,it)= IRKstep_fixed!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiters,
				     		 initial_interp,abstol,reltol,adaptive,
							 threading)
	      end
      end

   else  # (typeof(prob.f<:DynamicalODEFunction))

	   if (adaptive==true)
	       (status,it)= IRKstepDynODE_adaptive!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiters,
				        maxtrials,initial_interp,abstol,reltol,adaptive,threading)
	   else
		 (status,it)= IRKstepDynODE_fixed!(s,j,tj,uj,ej,prob,dts,coeffs,cache,maxiters,
					   initial_interp,abstol,reltol,adaptive, threading)
	   end

   end

   return(status,it)

end


function IRKstep_fixed!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiters,
		               initial_interp,abstol,reltol,adaptive,threading)

        @unpack mu,hc,hb,nu,alpha = coeffs
        @unpack f,u0,p,tspan=prob
		@unpack U,Uz,L,Lz,F,Dmin,Eval,DY,rejects,nfcn,lambdas=cache

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


    	iter = true
    	plusIt=true

    	nit=1
		@inbounds begin
		for is in 1:s Dmin[is] .= Inf end
	    end

       if threading==true
	   	@inbounds begin
	   		Threads.@threads for is in 1:s
		   		nfcn[1]+=1
	       		f(F[is], U[is], p, tj + hc[is])
		   		@. L[is] = hb[is]*F[is]
	   		end
       	end
       else
	     @inbounds begin
   	   		for is in 1:s
   		   		nfcn[1]+=1
   	       		f(F[is], U[is], p, tj + hc[is])
   		   		@. L[is] = hb[is]*F[is]
   	   		end
         end
	   end


       while (nit<maxiters && iter)

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
		  end  #inbounds

          if threading==true
	         Threads.@threads for is in 1:s
	     	    Eval[is]=false
		     	for k in eachindex(uj)
			    	DY[is]=abs(U[is][k]-Uz[is][k])
			    	if DY[is]>0.
			       		Eval[is]=true
			       		if DY[is]< Dmin[is][k]
				      		Dmin[is][k]=DY[is]
				      		iter=true
			      		end
		       		else
			      		D0+=1
		     		end
		    	end

		    	if Eval[is]==true
	    			nfcn[1]+=1
		    		f(F[is], U[is],p, tj + hc[is])
			    	@. L[is] = hb[is]*F[is]
		   		end
	      	end

		 else  #threading=false
			for is in 1:s
			   Eval[is]=false
			   for k in eachindex(uj)
			      DY[is]=abs(U[is][k]-Uz[is][k])
			      if DY[is]>0.
				     Eval[is]=true
				     if DY[is]< Dmin[is][k]
					    Dmin[is][k]=DY[is]
					    iter=true
				     end
			      else
				     D0+=1
			     end
		      end

		      if Eval[is]==true
			     nfcn[1]+=1
			     f(F[is], U[is],p, tj + hc[is])
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

    end  # ehile

	if (uiType<:CompiledFloats)

        #~Compensated summation
    	indices = eachindex(uj)
    	@inbounds begin
		for k in indices
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


function IRKstep_fixed_Mix!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiters,
		               initial_interp,abstol,reltol,adaptive,threading,
					   mixed_precision,low_prec_type)

		@unpack mu,hc,hb,nu,alpha = coeffs
		@unpack f,u0,p,tspan,kwargs=prob

		@unpack U,Uz,L,Lz,F,Dmin,Eval,DY,rejects,nfcn,lambdas,
				Ulow,DU,DF,DL,delta,Fa,Fb,
				normU,lhb,lmu=cache

        if !isempty(kwargs) lpp=kwargs[:lpp]
		else
			lpp=[]
		end
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
			lhb.=hb
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

    	iter = true
    	plusIt=true

    	nit=1
		@inbounds begin
		for is in 1:s Dmin[is] .= Inf end
	    end


       if threading==true
	       @inbounds begin
	       Threads.@threads for is in 1:s
		      nfcn[1]+=1
	          f(F[is], U[is], p, tj + hc[is])
		      @. L[is] = hb[is]*F[is]
	       end
           end
	   else
		   @inbounds begin
		   for is in 1:s
			  nfcn[1]+=1
			  f(F[is], U[is], p, tj + hc[is])
			  @. L[is] = hb[is]*F[is]
		   end
		   end
       end

	   lmax=1


       while (nit<maxiters && iter)

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
			Ulow[is].=U[is]
			normU[is]=copy(norm(Ulow[is]))
	      end
      	  end # inbound

          if threading==true
		  	Threads.@threads for is in 1:s
  				Eval[is]=false
  				for k in eachindex(uj)
					DY[is]=abs(Rdigits(U[is][k],10)-Rdigits(Uz[is][k],10))
	  				if DY[is]>0.
						Eval[is]=true
						if DY[is]< Dmin[is][k]
							Dmin[is][k]=DY[is]
							iter=true
						end
	 	   	   		else
				   		D0+=1
			   		end
 		    	end

 	        	if Eval[is]==true
	  				nfcn[1]+=1
	  				f(F[is], U[is],p, tj + hc[is])
					@. delta[is] = muladd(F[is],hb[is],-L[is])
		    	else
			    	delta[is].=0
 		    	end

		    	DL[is].=delta[is]
	       end

	     else

			for is in 1:s
   				Eval[is]=false
   				for k in eachindex(uj)
 					DY[is]=abs(Rdigits(U[is][k],10)-Rdigits(Uz[is][k],10))
 	  				if DY[is]>0.
 						Eval[is]=true
 						if DY[is]< Dmin[is][k]
 							Dmin[is][k]=DY[is]
 							iter=true
 						end
 	 	   	   		else
 				   		D0+=1
 			   		end
  		    	end

  	        	if Eval[is]==true
 	  				nfcn[1]+=1
 	  				f(F[is], U[is],p, tj + hc[is])
 					@. delta[is] = muladd(F[is],hb[is],-L[is])
 		    	else
 			    	delta[is].=0
  		    	end

 		    	DL[is].=delta[is]
    	     end
	     end


		 lmax=min(lmax*2,6)

         for l in 1:lmax

	       for is in 1:s
			   if (Eval[is]==true)
				DiffEqBase.@.. DU[is] =  lmu[is,1]*DL[1]+lmu[is,2]*DL[2]+
		   					        	 lmu[is,3]*DL[3]+lmu[is,4]*DL[4]+
		   								 lmu[is,5]*DL[5]+lmu[is,6]*DL[6]+
		   								 lmu[is,7]*DL[7]+lmu[is,8]*DL[8]
			  end
           end

           if threading==true

	       	Threads.@threads for is in 1:s

                if (Eval[is]==true && norm(DU[is])!=0)

        			beta=1e-6*norm(normU[is])/norm(DU[is])
					nfcn[2]+=2
					tjci=convert(low_prec_type, tj + hc[is])
					f(Fa[is], muladd.(beta,DU[is],Ulow[is]),lpp,tjci)
					f(Fb[is], muladd.(-beta,DU[is],Ulow[is]),lpp,tjci)
        			@. DF[is]=1/(2*beta)*(Fa[is]-Fb[is])
					@. DL[is] = muladd(DF[is],lhb[is],delta[is])

		    	end
          	 end

		  else

			  for is in 1:s

                  if (Eval[is]==true && norm(DU[is])!=0)

          			beta=1e-6*norm(normU[is])/norm(DU[is])
  					nfcn[2]+=2
  					tjci=convert(low_prec_type, tj + hc[is])
  					f(Fa[is], muladd.(beta,DU[is],Ulow[is]),lpp,tjci)
  					f(Fb[is], muladd.(-beta,DU[is],Ulow[is]),lpp,tjci)
          			@. DF[is]=1/(2*beta)*(Fa[is]-Fb[is])
  					@. DL[is] = muladd(DF[is],lhb[is],delta[is])

  		    	end
             end

          end


		end # end  for l

        for is in 1:s
			if (Eval[is]==true)
			  @. L[is] +=DL[is]
		   end
	    end

	    if (iter==false && D0<elems && plusIt)
		    iter=true
	    	plusIt=false
	    else
		    plusIt=true
	    end

     end # while

    if (uiType<:CompiledFloats)

#     ~ Compensated summation

		indices = eachindex(uj)
    	@inbounds begin
		for k in indices
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




function IRKstep_adaptive!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiters,maxtrials,
		                   initial_interp,abstol,reltol,adaptive,threading)


		@unpack mu,hc,hb,nu,alpha = coeffs
        @unpack f,u0,p,tspan=prob
		@unpack U,Uz,L,Lz,F,Dmin,Eval,DY,rejects,nfcn,lambdas=cache

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

        	if initial_interp
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
			    end
        	else
				@inbounds begin
            	for is in 1:s
                	@. U[is] = uj + ej
            	end
			    end
        	end

            if threading==true
            	@inbounds begin
        		Threads.@threads for is in 1:s
					nfcn[1]+=1
            		f(F[is], U[is], p, tj + hc[is])
            		@. L[is] = hb[is]*F[is]
        		end
		    	end
			else
				@inbounds begin
        		for is in 1:s
					nfcn[1]+=1
            		f(F[is], U[is], p, tj + hc[is])
            		@. L[is] = hb[is]*F[is]
        		end
		    	end
			end

        	iter = true
        	plusIt=true

        	nit=1
			for is in 1:s Dmin[is] .= Inf end

        	while (nit<maxiters && iter)

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
				end #inbound


               if threading==true

            		Threads.@threads for is in 1:s
						Eval[is]=false
                		for k in eachindex(uj)
                        	DY[is]=abs(U[is][k]-Uz[is][k])
                        	if DY[is]>0.
						   		Eval[is]=true
                           		if DY[is]< Dmin[is][k]
                              		Dmin[is][k]=DY[is]
                              		iter=true
                           		end
                       		else
                           		D0+=1
                       		end
                		end

               	  		if Eval[is]==true
				   	 		nfcn[1]+=1
                  	 		f(F[is], U[is], p,  tj + hc[is])
                  	 		@. L[is] = hb[is]*F[is]
                  		end
	       			end

				else

					for is in 1:s
						Eval[is]=false
                		for k in eachindex(uj)
                        	DY[is]=abs(U[is][k]-Uz[is][k])
                        	if DY[is]>0.
						   		Eval[is]=true
                           		if DY[is]< Dmin[is][k]
                              		Dmin[is][k]=DY[is]
                              		iter=true
                           		end
                       		else
                           		D0+=1
                       		end
                		end

               	  		if Eval[is]==true
				   	 		nfcn[1]+=1
                  	 		f(F[is], U[is], p,  tj + hc[is])
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

#			~Compensated summation

			indices = eachindex(uj)
        	@inbounds begin
			for k in indices
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



function IRKstep_adaptive_Mix!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiters,maxtrials,
		                   initial_interp,abstol,reltol,adaptive,threading,
						   mixed_precision,low_prec_type)


		@unpack mu,hc,hb,nu,alpha = coeffs
        @unpack f,u0,p,tspan,kwargs=prob

		@unpack U,Uz,L,Lz,F,Dmin,Eval,DY,rejects,nfcn,lambdas,
	         	Ulow,DU,DF,DL,delta,Fa,Fb,
	         	normU,lhb,lmu=cache

        if !isempty(kwargs) lpp=kwargs[:lpp]
		else
			lpp=[]
		end
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
				lhb.=hb
			end

        	if initial_interp
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
			    end
        	else
				@inbounds begin
            	for is in 1:s
                	@. U[is] = uj + ej
            	end
			    end
        	end

            if threading==true
            	@inbounds begin
        		Threads.@threads for is in 1:s
					nfcn[1]+=1
            		f(F[is], U[is], p, tj + hc[is])
            		@. L[is] = hb[is]*F[is]
        		end
		    	end
			else
				@inbounds begin
        		for is in 1:s
					nfcn[1]+=1
            		f(F[is], U[is], p, tj + hc[is])
            		@. L[is] = hb[is]*F[is]
        		end
		    	end
			end

            lmax=1
        	iter = true
        	plusIt=true

        	nit=1
			for is in 1:s Dmin[is] .= Inf end

        	while (nit<maxiters && iter)

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
					Ulow[is].=U[is]
					normU[is]=copy(norm(Ulow[is]))
        		end
				end #inbound

				if threading==true
            		Threads.@threads for is in 1:s
						Eval[is]=false
                		for k in eachindex(uj)
							DY[is]=abs(Rdigits(U[is][k],10)-Rdigits(Uz[is][k],10))
                        	if DY[is]>0.
						   		Eval[is]=true
                           		if DY[is]< Dmin[is][k]
                              		Dmin[is][k]=DY[is]
                              	iter=true
                           		end
                       		else
                           		D0+=1
                       		end
                		end

               			if Eval[is]==true
							nfcn[1]+=1
                  			f(F[is], U[is], p,  tj + hc[is])
							@. delta[is] = muladd(F[is],hb[is],-L[is])
						else
							delta[is].=0
               			end

						DL[is].=delta[is]
           			end

			    else

					for is in 1:s
						Eval[is]=false
                		for k in eachindex(uj)
							DY[is]=abs(Rdigits(U[is][k],10)-Rdigits(Uz[is][k],10))
                        	if DY[is]>0.
						   		Eval[is]=true
                           		if DY[is]< Dmin[is][k]
                              		Dmin[is][k]=DY[is]
                              	iter=true
                           		end
                       		else
                           		D0+=1
                       		end
                		end

               			if Eval[is]==true
							nfcn[1]+=1
                  			f(F[is], U[is], p,  tj + hc[is])
							@. delta[is] = muladd(F[is],hb[is],-L[is])
						else
							delta[is].=0
               			end

						DL[is].=delta[is]
           			end

			    end

				lmax=min(lmax*2,6)

				for l in 1:lmax
					for is in 1:s
						if (Eval[is]==true)
							DiffEqBase.@.. DU[is] =lmu[is,1]*DL[1]+lmu[is,2]*DL[2]+
												   lmu[is,3]*DL[3]+lmu[is,4]*DL[4]+
												   lmu[is,5]*DL[5]+lmu[is,6]*DL[6]+
												   lmu[is,7]*DL[7]+lmu[is,8]*DL[8]
						end
					end

					if threading==true
						Threads.@threads for is in 1:s
							if (Eval[is]==true && norm(DU[is])!=0)
								beta=1e-6*norm(normU[is])/norm(DU[is])
								nfcn[2]+=2
								tjci=convert(low_prec_type, tj + hc[is])
								f(Fa[is], muladd.(beta,DU[is],Ulow[is]),lpp,tjci)
								f(Fb[is], muladd.(-beta,DU[is],Ulow[is]),lpp,tjci)
								@. DF[is]=1/(2*beta)*(Fa[is]-Fb[is])
								@. DL[is] = muladd(DF[is],lhb[is],delta[is])
							end
						end
					else

						for is in 1:s
							if (Eval[is]==true && norm(DU[is])!=0)
								beta=1e-6*norm(normU[is])/norm(DU[is])
								nfcn[2]+=2
								tjci=convert(low_prec_type, tj + hc[is])
								f(Fa[is], muladd.(beta,DU[is],Ulow[is]),lpp,tjci)
								f(Fb[is], muladd.(-beta,DU[is],Ulow[is]),lpp,tjci)
								@. DF[is]=1/(2*beta)*(Fa[is]-Fb[is])
								@. DL[is] = muladd(DF[is],lhb[is],delta[is])
							end
						end
					end

			     end # end for l

			     for is in 1:s
				   if (Eval[is]==true)
					   @. L[is] +=DL[is]
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

#           ~Compensated summation

			indices = eachindex(uj)

        	@inbounds begin
				for k in indices
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


#
# DynamicalODEFunction
#

function IRKstepDynODE_fixed!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiters,
		               initial_interp,abstol,reltol,adaptive,threading)

		@unpack mu,hc,hb,nu,alpha = coeffs
		@unpack tspan,p=prob
		f1=prob.f.f1
		f2=prob.f.f2
		@unpack U,Uz,L,Lz,F,Dmin,Eval,DY,rejects,nfcn,lambdas=cache

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

		iter = true
		plusIt=true

    	nit=1
		for is in 1:s Dmin[is] .= Inf end

		if threading==true
        	@inbounds begin
    		Threads.@threads for is in 1:s
				nfcn[1]+=1
        		f1(F[is].x[1], U[is].x[1],U[is].x[2], p, tj + hc[is])
				f2(F[is].x[2], U[is].x[1],U[is].x[2], p, tj + hc[is])
        		@. L[is] = hb[is]*F[is]
    		end
	    	end

		else

			@inbounds begin
    		for is in 1:s
				nfcn[1]+=1
        		f1(F[is].x[1], U[is].x[1],U[is].x[2], p, tj + hc[is])
				f2(F[is].x[2], U[is].x[1],U[is].x[2], p, tj + hc[is])
        		@. L[is] = hb[is]*F[is]
    		end
	    	end

		end


    	while (nit<maxiters && iter)

        	nit+=1
        	iter=false
        	D0=0

#           First part

            @inbounds begin
    		for is in 1:s
        		Uz[is].x[1] .= U[is].x[1]
        		DiffEqBase.@.. U[is].x[1] = uj.x[1] +
				              (ej.x[1]+mu[is,1]*L[1].x[1] + mu[is,2]*L[2].x[1]+
			                           mu[is,3]*L[3].x[1] + mu[is,4]*L[4].x[1]+
                                       mu[is,5]*L[5].x[1] + mu[is,6]*L[6].x[1]+
								       mu[is,7]*L[7].x[1] + mu[is,8]*L[8].x[1])
    		end
			end #inbounds

            if threading==true

        		Threads.@threads for is in 1:s
            		Eval[is]=false
            		for k in eachindex(U[is].x[1])
                    	DY[is]=abs(U[is].x[1][k]-Uz[is].x[1][k])
                    	if DY[is]>0.
                       		Eval[is]=true
                       		if DY[is]< Dmin[is].x[1][k]
                          		Dmin[is].x[1][k]=DY[is]
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

			else

				for is in 1:s
            		Eval[is]=false
            		for k in eachindex(U[is].x[1])
                    	DY[is]=abs(U[is].x[1][k]-Uz[is].x[1][k])
                    	if DY[is]>0.
                       		Eval[is]=true
                       		if DY[is]< Dmin[is].x[1][k]
                          		Dmin[is].x[1][k]=DY[is]
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
			end

#           Second part

			for is in 1:s
        		Uz[is].x[2] .= U[is].x[2]
        		DiffEqBase.@.. U[is].x[2] = uj.x[2] +
				              (ej.x[2]+mu[is,1]*L[1].x[2] + mu[is,2]*L[2].x[2]+
	                                   mu[is,3]*L[3].x[2] + mu[is,4]*L[4].x[2]+
                                       mu[is,5]*L[5].x[2] + mu[is,6]*L[6].x[2]+
    					               mu[is,7]*L[7].x[2] + mu[is,8]*L[8].x[2])
    		end

			if threading==true
        		Threads.@threads for is in 1:s
            		Eval[is]=false
        			for k in eachindex(U[is].x[2])
	                  	DY[is]=abs(U[is].x[2][k]-Uz[is].x[2][k])
	                  	if DY[is]>0.
	                       Eval[is]=true
                       	   if DY[is]< Dmin[is].x[2][k]
			                  Dmin[is].x[2][k]=DY[is]
	                          iter=true
	                       end
	                   else
     	                    D0+=1
		               end
	           		end

           			if Eval[is]==true
#						nfcn[1]+=1
              			f1(F[is].x[1], U[is].x[1], U[is].x[2], p,  tj + hc[is])
              			@. L[is].x[1] = hb[is]*F[is].x[1]
           			end
    			end

			else
				for is in 1:s
            		Eval[is]=false
        			for k in eachindex(U[is].x[2])
	                  	DY[is]=abs(U[is].x[2][k]-Uz[is].x[2][k])
	                  	if DY[is]>0.
	                       Eval[is]=true
                       	   if DY[is]< Dmin[is].x[2][k]
			                  Dmin[is].x[2][k]=DY[is]
	                          iter=true
	                       end
	                   else
     	                    D0+=1
		               end
	           		end

           			if Eval[is]==true
#						nfcn[1]+=1
              			f1(F[is].x[1], U[is].x[1], U[is].x[2], p,  tj + hc[is])
              			@. L[is].x[1] = hb[is]*F[is].x[1]
           			end
    			end
			end

        	if (iter==false && D0<elems && plusIt)
            	iter=true
            	plusIt=false
        	else
            	plusIt=true
        	end

    	end # while iter


	    if (uiType<:CompiledFloats)

#           ~Compensated summation

			indices = eachindex(uj)
        	@inbounds begin
			for k in indices
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



function IRKstepDynODE_adaptive!(s,j,ttj,uj,ej,prob,dts,coeffs,cache,maxiters,maxtrials,
		                      initial_interp,abstol,reltol,adaptive,threading)

		@unpack mu,hc,hb,nu,alpha = coeffs
		@unpack tspan,p=prob
		f1=prob.f.f1
		f2=prob.f.f2
		@unpack U,Uz,L,Lz,F,Dmin,Eval,DY,rejects,nfcn,lambdas=cache

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
		    end

			iter = true
			plusIt=true
        	nit=1
			for is in 1:s Dmin[is] .= Inf end

			if threading==true
        		Threads.@threads for is in 1:s
					nfcn[1]+=1
					f1(F[is].x[1], U[is].x[1],U[is].x[2], p, tj + hc[is])
					f2(F[is].x[2], U[is].x[1],U[is].x[2], p, tj + hc[is])
        			@. L[is] = hb[is]*F[is]
        		end
			else
				for is in 1:s
					nfcn[1]+=1
					f1(F[is].x[1], U[is].x[1],U[is].x[2], p, tj + hc[is])
					f2(F[is].x[2], U[is].x[1],U[is].x[2], p, tj + hc[is])
	        		@. L[is] = hb[is]*F[is]
	        	end
     		end

			while (nit<maxiters && iter)

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
			    end #inbound


				if threading==true

					Threads.@threads for is in 1:s
	                	Eval[is]=false
						for k in eachindex(U[is].x[1])
							DY[is]=abs(U[is].x[1][k]-Uz[is].x[1][k])
							if DY[is]>0.
                            	Eval[is]=true
								if DY[is]< Dmin[is].x[1][k]
									Dmin[is].x[1][k]=DY[is]
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

				else

					for is in 1:s
	                	Eval[is]=false
						for k in eachindex(U[is].x[1])
							DY[is]=abs(U[is].x[1][k]-Uz[is].x[1][k])
							if DY[is]>0.
                            	Eval[is]=true
								if DY[is]< Dmin[is].x[1][k]
									Dmin[is].x[1][k]=DY[is]
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
				end


#               Second part

                @inbounds begin
        		for is in 1:s
            		Uz[is].x[2] .= U[is].x[2]
            		DiffEqBase.@.. U[is].x[2] = uj.x[2]+(ej.x[2]+mu[is,1]*L[1].x[2] + mu[is,2]*L[2].x[2]+
				                           mu[is,3]*L[3].x[2] + mu[is,4]*L[4].x[2]+
                                           mu[is,5]*L[5].x[2] + mu[is,6]*L[6].x[2]+
									       mu[is,7]*L[7].x[2] + mu[is,8]*L[8].x[2])
        		end
    		    end #inbound

				if threading==true

					Threads.@threads for is in 1:s
						Eval[is]=false
						for k in eachindex(U[is].x[2])
							DY[is]=abs(U[is].x[2][k]-Uz[is].x[2][k])
							if DY[is]>0.
								Eval[is]=true
								if DY[is]< Dmin[is].x[2][k]
									Dmin[is].x[2][k]=DY[is]
									iter=true
								end
							else
								D0+=1
							end
				   		end

                   		if (Eval[is]==true)
#							nfcn[1]+=1
	            			f1(F[is].x[1], U[is].x[1],U[is].x[2], p,  tj + hc[is])
	            			@. L[is].x[1] = hb[is]*F[is].x[1]
				   		end
	       	   		end

		   	   else

				  for is in 1:s
   					Eval[is]=false
   					for k in eachindex(U[is].x[2])
   						DY[is]=abs(U[is].x[2][k]-Uz[is].x[2][k])
   						if DY[is]>0.
   							Eval[is]=true
   							if DY[is]< Dmin[is].x[2][k]
   								Dmin[is].x[2][k]=DY[is]
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
		   	   end

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

#			~ Compensated summation
			indices = eachindex(uj)
        	@inbounds begin
			for k in indices
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
