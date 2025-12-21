
#
#  Schwarzschild problem
#  (cartesian coordinates)
#  2025-07-08

function Ham_Schwarzschildv2(u, parms)
    E=parms[1]
    L=parms[2]
    β=parms[3]
    x=u[1]
    y=u[2]
    px=u[3]
    py=u[4]

    H1 = H1v2(u, parms)
    H2 = H2v2(u, parms)
    H3 = H3v2(u, parms)

    return H1+H2+H3
end

function H1v2(u, parms)
    E=parms[1]
    L=parms[2]
    β=parms[3]
    x=u[1]
    y=u[2]
    px=u[3]
    py=u[4]

    x2=x*x
    y2=y*y
    r2=x2+y2
    E2=E*E
    r=sqrt(r2)

    w = L-β/2*y2

    H1=1/(2*y2)*w^2-1/2*(r/(r-2))*E2

    return H1
end

function H2v2(u, parms)
    E=parms[1]
    L=parms[2]
    β=parms[3]
    x=u[1]
    y=u[2]
    px=u[3]
    py=u[4]

    H2 = 1/2*(px*px+py*py)

    return H2
end

function H3v2(u, parms)
    E=parms[1]
    L=parms[2]
    β=parms[3]
    x=u[1]
    y=u[2]
    px=u[3]
    py=u[4]

    x2=x*x
    y2=y*y
    r2=x2+y2
    r=sqrt(r2)

    c=x/r
    s=y/r

    H3 = -(px*c+py*s)^2/r
end

function flowH1Schwarzschildv2!(uj, ej, h, parms)
    E=parms[1]
    L=parms[2]
    β=parms[3]

    x=uj[1]
    y=uj[2]
    px=uj[3]
    py=uj[4]

    z = y*y
    r2 = x*x+z
    r=sqrt(r2)
    w = L-β/2*z
    u = 1-2/r
    W = -w/z
    U = -E*E/(2*u*u)
    R = U/(r*r2)
    Z = R+(W-β)*W/2

    gradpx=-2*x*R
    gradpy=-2*y*Z

    uj[3]=px-h*gradpx
    uj[4]=py-h*gradpy

    return nothing
end

function flowH2Schwarzschildv2!(uj, ej, h, parms)
    x=uj[1]
    y=uj[2]
    px=uj[3]
    py=uj[4]

    uj[1]=x+h*px
    uj[2]=y+h*py

    return nothing
end

function flowH3Schwarzschildv2!(uj, ej, h, parms)
    x=uj[1]
    y=uj[2]
    px=uj[3]
    py=uj[4]

    r2=x^2+y^2
    r=sqrt(r2)

    c=x/r
    s=y/r
    ptheta=-y*px+x*py
    nu = x*px+y*py
    mu=nu*nu/(r*r2)

    nu_=nu-3*h*mu
    pr_=cbrt(mu*nu_)
    r_=nu_/pr_

    uj[1]=r_*c
    uj[2]=r_*s
    ptr=ptheta/r_
    uj[3]=c*pr_-s*ptr
    uj[4]=s*pr_+c*ptr

    return nothing
end
