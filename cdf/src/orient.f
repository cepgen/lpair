      SUBROUTINE orient(s,v1,v2,v3,v4,v5,dj,nopt,y)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,
     a     ct5,st5,cp3,sp3,cp5,sp5
      COMMON/variac/al3,al4,be4,be5,de3,de5,pp3,pp4,pp5
      COMMON/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
      COMMON/extra/s1,s2,t1,t2
      COMMON/levi/gram,dd1,dd2,dd3,dd4,dd5,delta,g4,sa1,sa2
      COMMON/dotp/p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,q1,q2
      DIMENSION y(4)
      CALL pickin(s,v1,v2,v3,v4,v5,dj,nopt,y)
      IF(dj.EQ.0)GO TO 10
      e=dsqrt(s)
      re=0.5/e
      e1=re*(s+w12)
      e2=re*(s-w12)
      p=re*sl1
      de3=re*(s2-w3+w12)
      de5=re*(s1-w5-w12)
      e3=e1-de3
      IF(e3.gt.e1) GOTO 10      ! Laurent workaround
      e4=de3+de5
      e5=e2-de5
      IF(e5.gt.e2) GOTO 10      ! Laurent workaround
      IF(e4.LT.v4) GO TO 10
      p3=dsqrt(e3*e3-w3)
      p4=dsqrt((e4-v4)*(e4+v4))
      IF(p4.EQ.0) GO TO 10
      p5=dsqrt(e5*e5-w5)
      pp3=dsqrt(dd1/s)/p
      pp5=dsqrt(dd3/s)/p
      st3=pp3/p3
      st5=pp5/p5
      IF(st3.GT.1..OR.st5.GT.1.) GO TO 10
      ct3=dsqrt(1.-st3*st3)
      ct5=dsqrt(1.-st5*st5)
      IF(e1*e3.LT.p13)ct3=-ct3
      IF(e2*e5.GT.p25)ct5=-ct5
      al3=st3*st3/(1.+ct3)
      be5=st5*st5/(1.-ct5)
      IF(dd5.LT.0) GO TO 10
      pp4=dsqrt(dd5/s)/p
      st4=pp4/p4
      IF(st4.GT.1.) GO TO 10
      ct4=dsqrt(1.-st4*st4)
      IF(e1*e4.LT.p14)ct4=-ct4
      al4=1.-ct4
      be4=1.+ct4
      IF(ct4.LT.0)be4=st4*st4/al4
      IF(ct4.ge.0)al4=st4*st4/be4
      rr=dsqrt(-gram/s)/(p*pp4)
      sp3=rr/pp3
      sp5=-rr/pp5
      IF(dabs(sp3).GT.1..OR.dabs(sp5).GT.1.) GO TO 10
      cp3=-dsqrt(1.-sp3*sp3)
      cp5=-dsqrt(1.-sp5*sp5)
      a1=pp3*cp3-pp5*cp5
      IF(dabs(pp4+pp3*cp3+cp5*pp5).LT.dabs(dabs(a1)-pp4)) GO TO 1
      IF(a1.LT.0)cp5=-cp5
      IF(a1.ge.0)cp3=-cp3
 1    RETURN
 10   dj=0.
      RETURN
      END
