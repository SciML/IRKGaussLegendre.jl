function InitialNBody15(T=Float64)

   """
        Initial Value from DE430
        Given in the next order:
          Sun,Mercury,Venus,Earth+Moon,Mars,Jupiter,Saturn,Uranus,Neptune,Pluto,
          Ceres, Pallas, Vesta, Iris, Bamberga
   """

Gm = [0.295912208285591100e-3,
      0.491248045036476000e-10,
      0.724345233264412000e-9,
      0.888769244512563400e-9+0.109318945074237400e-10,
      0.954954869555077000e-10,
      0.282534584083387000e-6,
      0.845970607324503000e-7,
      0.129202482578296000e-7,
      0.152435734788511000e-7,
      0.217844105197418000e-11,
      # Ceres, Pallas, Vesta, Iris, Bamberga
      0.140047655617234400e-12,
      0.310444819893871300e-13,
      0.385475018780881000e-13,
      0.213643444257140700e-14,
      0.138862658985619900e-14
  ]

N=length(Gm)

q = [0.00450250878464055477, 0.00076707642709100705, 0.00026605791776697764,
     0.36176271656028195477, -0.09078197215676599295,-0.08571497256275117236,
     0.61275194083507215477, -0.34836536903362219295,-0.19527828667594382236,
     0.12051741410138465477, -0.92583847476914859295,-0.40154022645315222236,
    -0.11018607714879824523, -1.32759945030298299295,-0.60588914048429142236,
    -5.37970676855393644523, -0.83048132656339789295,-0.22482887442656542236,
    7.89439068290953155477, 4.59647805517127300705, 1.55869584283189997764,
   -18.26540225387235944523, -1.16195541867586999295,-0.25010605772133802236,
   -16.05503578023336944523, -23.94219155985470899295,-9.40015796880239402236,
   -30.48331376718383944523, -0.87240555684104999295, 8.91157617249954997764
]

v = [-0.00000035174953607552, 0.00000517762640983341, 0.00000222910217891203,
      0.00336749397200575848, 0.02489452055768343341, 0.01294630040970409203,
      0.01095206842352823448, 0.01561768426786768341, 0.00633110570297786403,
      0.01681126830978379448, 0.00174830923073434441, 0.00075820289738312913,
      0.01448165305704756448, 0.00024246307683646861, -0.00028152072792433877,
      0.00109201259423733748, -0.00651811661280738459,-0.00282078276229867897,
      -0.00321755651650091552, 0.00433581034174662541, 0.00192864631686015503,
      0.00022119039101561468, -0.00376247500810884459,-0.00165101502742994997,
      0.00264276984798005548, -0.00149831255054097759,-0.00067904196080291327,
      0.00032220737349778078, -0.00314357639364532859,-0.00107794975959731297
]

qq = copy(q)
vv = copy(v)


# Ceres, Pallas, Vesta, Iris, Bamberga
qsun=qq[1:3]
qCeres=qsun+[1.438681809676469747, -2.204373633189407045, -1.326397853361325874]
qPallas=qsun+[0.203832272462290465, -3.209619436062307152, 0.623843179079393351]
qVesta=qsun+[0.182371836377417107,  2.386628211277654010, 0.924596062836265498]
qIris=qsun+[1.892475267790300286, -0.848414748075139946, -0.157159319044464590]
qBamberga=qsun+[1.398759064223541682, -1.287476729008325105, -0.669098428660833799]

vsun=vv[1:3]
vCeres=vsun+[0.008465406136316316, 0.004684247977335608, 0.000466157738595739]
vPallas=vsun+[0.008534313855651248, -0.000860659210123161, -0.000392901992572746]
vVesta=vsun+[-0.010174496747119257, 0.000041478190529952, 0.001344157634155624]
vIris=vsun+[0.002786950314570632, 0.011314057384917047, 0.004975132577079665]
vBamberga=vsun+[0.007164363244556328, 0.009219958777618218, 0.006857861727407507]


qq = [qq; qCeres; qPallas; qVesta; qIris; qBamberga]
vv = [vv; vCeres; vPallas; vVesta; vIris; vBamberga]

q0 = reshape(qq,3,:)
v0 = reshape(vv,3,:)

q0bar = [sum(Gm .* q0[j,:])/sum(Gm) for j in 1:3]
v0bar = [sum(Gm .* v0[j,:])/sum(Gm) for j in 1:3]

q0 = q0 .- q0bar
v0 = v0 .- v0bar

u0 = Array{T}(undef,2,3,N)
u0[1,:,:] = v0
u0[2,:,:] = q0
#u0[2,:,:] = v0
#u0[1,:,:] = q0
return u0, Gm
end