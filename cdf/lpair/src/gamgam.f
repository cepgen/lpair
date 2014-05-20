       SUBROUTINE gamgam(ebeam,v1,v2,v3,v5,v6,v7,vmin,vmax,dj,nopt,x,nm)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       COMMON/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5,
     a st5,cp3,sp3,cp5,sp5
       COMMON/variac/al3,al4,be4,be5,de3,de5,pp3,pp4,pp5
       COMMON/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
       COMMON/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
       COMMON/dotps/q1dq,q1dq2,w6
       COMMON /extra/s1,s2,t1,t2
       COMMON/dotp/p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,p1k2,p2k1
       COMMON/civita/epsi,g5,g6,a5,a6,bb
       COMMON/ext/ctg,stg,cpg,spg
       COMMON/angu/ctcm6,stcm6
       COMMON/qvec/qve(4)
       DIMENSION x(7)
       data pi/3.14159265358979d+00/,const/2.1868465d+10/
       w6=v6*v6
       w7=v7*v7
       wmin=v6+v7
       IF(wmin.LT.vmin)wmin=vmin
       wmin=wmin*wmin
       e=2.*ebeam
       s=e*e
       wmax=e-v3-v5
       IF(wmax.GT.vmax)wmax=vmax
       wmax=wmax*wmax
       xw=x(5)
       CALL mapw2(w4,xw,wmin,wmax,dw)
       v4=dsqrt(w4)
       w=v4
       CALL orient(s,v1,v2,v3,v4,v5,dj,nopt,x)
       IF(t1.GT.0.OR.t2.GT.0)dj=0.
       IF(dj.EQ.0) RETURN
       ecm6=(w4+w6-w7)/(2.*v4)
       pcm6=dsqrt(ecm6*ecm6-w6)
       dj=dj*dw*pcm6/(v4*const*s)
       e3mp3=w3/(e3+p3)
       e1mp1=w1/(e1+p)
       eg=(w4+t1-t2)/(2.*v4)
       pg=dsqrt(eg*eg-t1)
       pgx=-pp3*cp3*ct4-st4*(de3-e1mp1+e3mp3+p3*al3)
       pgy=-pp3*sp3
       pgz=v4*de3/(e4+p4)-e4*de3*al4/v4-pp3*cp3
     a      *e4*st4/v4+e4*ct4/v4*(p3*al3+e3mp3-e1mp1)
       pgp=dsqrt(pgx*pgx+pgy*pgy)
       pgg=dsqrt(pgp*pgp+pgz*pgz)
       IF(pgg.GT.pgp*0.9.AND.pgg.GT.pg)pg=pgg
       stg=pgp/pg
       cpg=pgx/pgp
       spg=pgy/pgp
       ctg=dsqrt(1.-stg*stg)
       IF(pgz.LT.0)ctg=-ctg
       xx6=x(6)
       IF(nm.EQ.0) GO TO 1
       amap=.5*(w4-t1-t2)
       bmap=.5*dsqrt(((w4-t1-t2)**2-4.*t1*t2)*(1.-4.*w6/w4))
       ymap=(amap+bmap)/(amap-bmap)
       beta=ymap**(2.*xx6-1.)
       xx6=(amap/bmap*(beta-1.)/(beta+1.)+1.)*0.5
       IF(xx6.GT.1.)xx6=1.
       IF(xx6.LT.0.)xx6=0.
       ctcm6=1.-2.*xx6
       ddd=(amap+bmap*ctcm6)*(amap-bmap*ctcm6)/amap/bmap*dlog(ymap)
       dj=dj*ddd*0.5
 1     ctcm6=1.-2.*xx6
       stcm6=2.*dsqrt(xx6*(1.-xx6))
       phicm6=2.*pi*x(7)
       cpcm6=dcos(phicm6)
       spcm6=dsin(phicm6)
       pcm6x=pcm6*stcm6*cpcm6
       pcm6y=pcm6*stcm6*spcm6
       pcm6z=pcm6*ctcm6
       pc6z=ctg*pcm6z-stg*pcm6x
       h1=stg*pcm6z+ctg*pcm6x
       pc6x=cpg*h1-spg*pcm6y
       qcx=2.*pc6x
       qcz=2.*pc6z
       p6y=cpg*pcm6y+spg*h1
       e6=(e4*ecm6+p4*pc6z)/v4
       h2=(e4*pc6z+p4*ecm6)/v4
       p6x=ct4*pc6x+st4*h2
       p6z=ct4*h2-st4*pc6x
       qve(1)=p4*qcz/v4
       qve(3)=2.*p6y
       hq=e4*qcz/v4
       qve(2)=ct4*qcx+st4*hq
       qve(4)=ct4*hq-st4*qcx
       p6=dsqrt(e6*e6-w6)
       e7=e4-e6
       p7=dsqrt(e7*e7-w7)
       p7x=pp4-p6x
       p7y=-p6y
       p7z=p4*ct4-p6z
       pp6=dsqrt(p6x*p6x+p6y*p6y)
       pp7=dsqrt(p7x*p7x+p7y*p7y)
       ct6=p6z/p6
       st6=pp6/p6
       ct7=p7z/p7
       st7=pp7/p7
       cp6=p6x/pp6
       sp6=p6y/pp6
       cp7=p7x/pp7
       sp7=p7y/pp7
       q1dq=eg*(2.*ecm6-v4)-2.*pg*pcm6*ctcm6
       q1dq2=0.5*(w4-t1-t2)
       bb=t1*t2+(w4*stcm6*stcm6+4.*w6*ctcm6*ctcm6)*pg*pg
       q0=qve(1)
       qx=qve(2)
       qy=qve(3)
       qz=qve(4)
       c1=(qx*sp3-qy*cp3)*pp3
       c2=(qz*e1-q0*p)*pp3
       c3=(w31*e1*e1+2.*w1*de3*e1-w1*de3*de3+pp3*pp3*e1*e1)
     a      /(e3*p+p3*ct3*e1)
       b1=(qx*sp5-qy*cp5)*pp5
       b2=(qz*e2+q0*p)*pp5
       b3=(w52*e2*e2+2.*w2*de5*e2-w2*de5*de5+pp5*pp5*e2*e2)
     a      /(e2*p5*ct5-e5*p)
       r12=c2*sp3+qy*c3
       r13=-c2*cp3-qx*c3
       r22=b2*sp5+qy*b3
       r23=-b2*cp5-qx*b3
       epsi=p12*c1*b1+r12*r22+r13*r23
       g5=w1*c1*c1+r12*r12+r13*r13
       g6=w2*b1*b1+r22*r22+r23*r23
       a5=-(qx*cp3+qy*sp3)*pp3*p1k2-(e1*q0-p*qz)*(cp3*cp5+sp3*sp5)
     a      *pp3*pp5+(de5*qz+q0*(p+p5*ct5))*c3
       a6=-(qx*cp5+qy*sp5)*pp5*p2k1-(e2*q0+p*qz)*(cp3*cp5+sp3*sp5)
     a      *pp3*pp5+(de3*qz-q0*(p-p3*ct3))*b3
C     LF FIXME... introduced to filter out these negative-energy internal photon 1 events
c       IF(E3.GT.E1) dj=0.
c       IF(E5.Gt.E2) dj=0.
C     LF FIXME
       RETURN
       END
