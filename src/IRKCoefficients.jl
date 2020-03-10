function PolInterp(X::AbstractVector{ctype}, Y::AbstractMatrix{ctype}, Z::AbstractVector{ctype}) where {ctype}
    N = length(X)
    M = length(Z)
    K = size(Y,1)
    if size(Y,2)!=N
        error("columns(Y) != length(X)")
    end
    pz = zeros(ctype,K,M)
    for i=1:N
        lag = 1.
        for j=1:N
            if (j!=i)
                lag *= X[i]-X[j]
            end
        end
        lag = 1/lag
        for m=1:M
            liz = lag
            for j=1:N
                if (j!=i)
                    liz *= Z[m]-X[j]
                end
            end
            for k=1:K
                pz[k,m] += Y[k,i]*liz
            end
        end
    end
    return pz
end


function IRK8Coefficients(h::AbstractFloat, hprev::AbstractFloat)

      mu = convert.(typeof(h), [
             # s=1
             parse(BigFloat,"0.5") parse(BigFloat,"-8.1894963105581497136508164735930892e-02") parse(BigFloat,"+4.0042703777945052339969045601553216e-02") parse(BigFloat,"-2.4721345803200374868044581645082913e-02") parse(BigFloat,"+1.6976173236371093026881708761116103e-02 ") parse(BigFloat,"-1.2225914113298206053942531721943548e-02") parse(BigFloat,"+8.7485667691973688117766530139122287e-03") parse(BigFloat,"-5.4828082532158826793409353214950813e-03")
             # s=2
             parse(BigFloat,"+1.0818949631055814971365081647359309e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.6958924300832723329070964616248016e-02") parse(BigFloat,"+4.4941126302625688139830943466131237e-02") parse(BigFloat,"-2.8759775474749310978230557041068536e-02")  parse(BigFloat,"+2.0017127636408709173710417097426705e-02") parse(BigFloat,"-1.4074355889166929145973516652599415e-02") parse(BigFloat,"+8.7485667691973688117766530139122287e-03")
             # s=3
             parse(BigFloat,"+9.5995729622205494766003095439844678e-01") parse(BigFloat,"+1.0869589243008327233290709646162480e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.8093838732308313442213871391320316e-02") parse(BigFloat,"+4.6165464814800034116730885592456977e-02")  parse(BigFloat,"-2.9603873064977937463012598212122325e-02") parse(BigFloat,"+2.0017127636408709173710417097426705e-02") parse(BigFloat,"-1.2225914113298206053942531721943548e-02")
             # s=4
             parse(BigFloat,"+1.0247213458032003748680445816450829e+00") parse(BigFloat,"+9.5505887369737431186016905653386876e-01") parse(BigFloat,"+1.0880938387323083134422138713913203e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.8347161109827784250707380600804499e-02")  parse(BigFloat,"+4.6165464814800034116730885592456977e-02") parse(BigFloat,"-2.8759775474749310978230557041068536e-02") parse(BigFloat,"+1.6976173236371093026881708761116103e-02")
             # s=5
             parse(BigFloat,"+9.8302382676362890697311829123888390e-01") parse(BigFloat,"+1.0287597754747493109782305570410685e+00 ") parse(BigFloat,"+9.5383453518519996588326911440754302e-01") parse(BigFloat,"+1.0883471611098277842507073806008045e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.8093838732308313442213871391320316e-02") parse(BigFloat,"+4.4941126302625688139830943466131237e-02") parse(BigFloat,"-2.4721345803200374868044581645082913e-02")
             # s=6
             parse(BigFloat,"+1.0122259141132982060539425317219435e+00") parse(BigFloat,"+9.7998287236359129082628958290257329e-01")  parse(BigFloat,"+1.0296038730649779374630125982121223e+00")  parse(BigFloat,"+9.5383453518519996588326911440754302e-01") parse(BigFloat,"+1.0880938387323083134422138713913203e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.6958924300832723329070964616248016e-02") parse(BigFloat,"+4.0042703777945052339969045601553216e-02")
             # s=7
             parse(BigFloat,"+9.9125143323080263118822334698608777e-01") parse(BigFloat,"+1.0140743558891669291459735166525994e+00")  parse(BigFloat,"+9.7998287236359129082628958290257329e-01")  parse(BigFloat,"+1.0287597754747493109782305570410685e+00") parse(BigFloat,"+9.5505887369737431186016905653386876e-01 ") parse(BigFloat,"+1.0869589243008327233290709646162480e+00") parse(BigFloat,"0.5") parse(BigFloat,"-8.1894963105581497136508164735930892e-02")
             # s=8
             parse(BigFloat,"+1.0054828082532158826793409353214951e+00") parse(BigFloat,"+9.9125143323080263118822334698608777e-01")  parse(BigFloat,"+1.0122259141132982060539425317219435e+00")  parse(BigFloat,"+9.8302382676362890697311829123888390e-01") parse(BigFloat,"+1.0247213458032003748680445816450829e+00")  parse(BigFloat,"+9.5995729622205494766003095439844678e-01") parse(BigFloat,"+1.0818949631055814971365081647359309e+00") parse(BigFloat,"0.5")
          ])

          b = convert.(typeof(h),
#             h*[ parse(BigFloat,"+5.0614268145188129576265677154981094e-02"),
                [ parse(BigFloat,"+5.0614268145188129576265677154981094e-02"),
                 parse(BigFloat,"+1.1119051722668723527217799721312045e-01"),
                 parse(BigFloat,"+1.5685332293894364366898110099330067e-01"),
                 parse(BigFloat,"+1.8134189168918099148257522463859781e-01"),
                 parse(BigFloat,"+1.8134189168918099148257522463859781e-01"),
                 parse(BigFloat,"+1.5685332293894364366898110099330067e-01"),
                 parse(BigFloat,"+1.1119051722668723527217799721312045e-01"),
                 parse(BigFloat,"+5.0614268145188129576265677154981094e-02")
               ])


          c= convert.(typeof(h),
#             h*[ parse(BigFloat,"+1.9855071751231884158219565715263505e-02"),
                [ parse(BigFloat,"+1.9855071751231884158219565715263505e-02"),
                 parse(BigFloat,"+1.0166676129318663020422303176208480e-01"),
                 parse(BigFloat,"+2.3723379504183550709113047540537686e-01"),
                 parse(BigFloat,"+4.0828267875217509753026192881990801e-01"),
                 parse(BigFloat,"+5.9171732124782490246973807118009203e-01"),
                 parse(BigFloat,"+7.6276620495816449290886952459462321e-01"),
                 parse(BigFloat,"+8.9833323870681336979577696823791522e-01"),
                 parse(BigFloat,"+9.8014492824876811584178043428473653e-01")
          ])


  s = length(b)


""" \tilde b_i coefficients i=2,s """

  B=[1/(k+1) for k in 0:s-2]
  M=[ones(s-1)';[c[i]^j for i in 2:s, j in 1:s-2]']
  tb=M\B
  beta=b[1:s]-vcat(0,tb)

""" mu & hb """

  for i in 1:s
      for j in i+1:s
          mu[i,j] = 1 - mu[j,i]
      end
  end

  hb=h*b
  hb1 = (h-sum(hb[2:end-1]))/2
  hb[1] = hb1
  hb[end] = hb1

  hc=h*c


 """ Interpolate coefficients """

  if (hprev==0.)
      nu=zeros(s,s)
  else
      lambda=h/hprev
      X=vcat(-c[end:-1:1],[0])
      Y=hcat(mu,zeros(s))
      Z=lambda*c
      nu=-PolInterp(X,Y,Z)
  end


return (hb = hb, hc = hc, mu = mu, nu=nu', beta=beta)

end
