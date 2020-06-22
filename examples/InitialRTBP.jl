function InitialRTBP(T=Float64)

   μ = parse(BigFloat,"0.012277471")

   q=[parse(BigFloat,"0.994"),parse(BigFloat,"0")]
   p=[parse(BigFloat,"0"),parse(BigFloat,"-2.00158510637908252240537862224")+0.994]

   u0=convert.(T,[q;p])

   return u0,convert(T,μ)

end
