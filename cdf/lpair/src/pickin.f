      SUBROUTINE pickin(s,v1,v2,v3,v4,v5,dj,nopt,y)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION y(4)
      COMMON/pickzz/w1,w2,w3,w4,w5,d1,d2,d5,d7,sl1
      COMMON/extra/s1,s2,t1,t2
      COMMON/accura/acc3,acc4
      COMMON/levi/gram,dd1,dd2,dd3,dd4,dd5,delta,g4,sa1,sa2
      COMMON/dotp/p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,p1k2,p2k1
      data pi/3.14159265358979d+00/
      x1=y(1)
      x2=y(2)
      x3=y(3)
      w1=v1*v1
      w2=v2*v2
      w3=v3*v3
      w4=v4*v4
      w5=v5*v5
      sig=v4+v5
      sig1=sig*sig
      sig2=sig1
      d1=w3-w1
      d2=w5-w2
      d5=w1-w2
      d6=w4-w5
      ss=s+d5
      rl1=ss*ss-4.*w1*s
      IF(rl1.le.0) GO TO 20
      sl1=dsqrt(rl1)
      IF(nopt.ne.0) GO TO 1
      smax=s+w3-2.*v3*dsqrt(s)
      CALL maps2(s2,x3,sig1,smax,ds2)
      sig1=s2
 1    sp=s+w3-sig1
      d3=sig1-w2
      rl2=sp*sp-4.*s*w3
      IF(rl2.le.0) GO TO 20
      sl2=sqrt(rl2)
      t1max=w1+w3-(ss*sp+sl1*sl2)/(2.*s)
      t1min=(d1*d3+(d3-d1)*(d3*w1-d1*w2)/s)/t1max
      CALL mapt1(t1,x1,t1min,t1max,dt1)
      d4=w4-t1
      d8=t1-w2
      t13=t1-w1-w3
      sa1=-(t1-d1)*(t1-d1)*0.25+w1*t1
      IF(sa1.ge.0) GO TO 20
      sl3=dsqrt(-sa1)
      IF(w1.EQ.0) GO TO 3
      sb=(s*(t1-d1)+d5*t13)/(2.*w1)+w3
      sd=sl1*sl3/w1
      se=(s*(t1*(s+t13-w2)-w2*d1)+w3*(d5*d8+w2*w3))/w1
      IF(dabs((sb-sd)/sd).LT.1.0) GO TO 2
      splus=sb-sd
      s2max=se/splus
      GO TO 4
 2    s2max=sb+sd
      splus=se/s2max
      GO TO 4
 3    s2max=(s*(t1*(s+d8-w3)-w2*w3)+w2*w3*(w2+w3-t1))/ss/t13
      splus=sig2
 4    s2x=s2max
      IF(nopt)5,6,7
 5    IF(splus.GT.sig2)sig2=splus
      IF(nopt.LT.-1)CALL maps2(s2,x3,sig2,s2max,ds2)
      IF(nopt.EQ.-1)CALL mapla(s2,t1,w2,x3,sig2,s2max,ds2)
 6    s2x=s2
 7    r1=s2x-d8
      r2=s2x-d6
      rl4=(r1*r1-4.*w2*s2x)*(r2*r2-4.*w5*s2x)
      IF(rl4.le.0) GO TO 20
      sl4=dsqrt(rl4)
      t2max=w2+w5-(r1*r2+sl4)/(2.*s2x)
      t2min=(d2*d4+(d4-d2)*(d4*w2-d2*t1)/s2x)/t2max
      CALL mapt2(t2,x2,t2min,t2max,dt2)
      d7=t1-t2
      r3=d4-t2
      r4=d2-t2
      b=r3*r4-2.*(t1+w2)*t2
      c=t2*d6*d8+(d6-d8)*(d6*w2-d8*w5)
      t25=t2-w2-w5
      sa2=-r4*r4*0.25+w2*t2
      IF(sa2.ge.0) GO TO 20
      sl6=2.*dsqrt(-sa2)
      g4=-0.25*r3*r3+t1*t2
      IF(g4.ge.0) GO TO 20
      sl7=sqrt(-g4)*2.
      sl5=sl6*sl7
      IF(dabs((sl5-b)/sl5).LT.1.0) GO TO 8
      s2p=(sl5-b)/(2.*t2)
      s2min=c/(t2*s2p)
      GO TO 9
 8    s2min=(-sl5-b)/(2.*t2)
      s2p=c/(t2*s2min)
 9    IF(nopt.GT.1)CALL maps2(s2,x3,s2min,s2max,ds2)
      IF(nopt.EQ.1)CALL mapla(s2,t1,w2,x3,s2min,s2max,ds2)
      ap=-(s2+d8)*(s2+d8)*0.25+s2*t1
      IF(w1.EQ.0) GO TO 10
      dd1=-w1*(s2-s2max)*(s2-splus)*0.25
      GO TO 11
 10   dd1=ss*t13*(s2-s2max)*0.25
 11   dd2=-t2*(s2-s2p)*(s2-s2min)*0.25
      yy4=dcos(pi*y(4))
      dd=dd1*dd2
      p12=0.5*(s-w1-w2)
      st=s2-t1-w2
      delb=(2.*w2*r3+r4*st)*(4.*p12*t1-(t1-d1)*st)/(16.*ap)
      IF(dd.le.0) GO TO 20
      delta=delb-yy4*st*dsqrt(dd)/(2.*ap)
      s1=t2+w1+(2.*p12*r3-4.*delta)/st
      IF(ap.ge.0) GO TO 20
      dj=ds2*dt1*dt2*pi*pi/(8.*sl1*dsqrt(-ap))
      gram=(1.-yy4)*(1.+yy4)*dd/ap
      p13=-t13*0.5
      p14=(d7+s1-w3)*0.5
      p15=(s+t2-s1-w2)*0.5
      p23=(s+t1-s2-w1)*0.5
      p24=(s2-d7-w5)*0.5
      p25=-t25*0.5
      p34=(s1-w3-w4)*0.5
      p35=(s+w4-s1-s2)*0.5
      p45=(s2-w4-w5)*0.5
      p1k2=(s1-t2-w1)*0.5
      p2k1=st*0.5
      IF(w2.EQ.0) GO TO 14
      sbb=(s*(t2-d2)-d5*t25)/(2.*w2)+w5
      sdd=sl1*sl6/(2.*w2)
      see=(s*(t2*(s+t25-w1)-w1*d2)+w5*(w1*w5-d5*(t2-w1)))/w2
      IF(sbb/sdd.LT.0) GO TO 12
      s1p=sbb+sdd
      s1m=see/s1p
      GO TO 13
 12   s1m=sbb-sdd
      s1p=see/s1m
 13   dd3=-w2*(s1p-s1)*(s1m-s1)*0.25
      GO TO 15
 14   s1p=(s*(t2*(s-w5+t2-w1)-w1*w5)+w1*w5*(w1+w5-t2))/t25/(s-d5)
      dd3=-t25*(s-d5)*(s1p-s1)*0.25
 15   acc3=(s1p-s1)/(s1p+s1)
      ssb=t2+w1-r3*(d1-t1)*0.5/t1
      ssd=sl3*sl7/t1
      sse=(t2-w1)*(w4-w3)+(t2-w4+d1)*((t2-w1)*w3-(w4-w3)*w1)/t1
      IF(ssb/ssd.LT.0) GO TO 16
      s1pp=ssb+ssd
      s1pm=sse/s1pp
      GO TO 17
 16   s1pm=ssb-ssd
      s1pp=sse/s1pm
 17   dd4=-t1*(s1-s1pp)*(s1-s1pm)*0.25
      acc4=(s1-s1pm)/(s1+s1pm)
      dd5=dd1+dd3+((p12*(t1-d1)*0.5-w1*p2k1)*(p2k1*(t2-d2)-w2*r3)
     a     -delta*(2.*p12*p2k1-w2*(t1-d1)))/p2k1
      RETURN
 20   dj=0.
      RETURN
      END
