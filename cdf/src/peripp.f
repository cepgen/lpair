       DOUBLE PRECISION FUNCTION peripp(nup,ndown)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       COMMON/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
       COMMON/levi/gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
       COMMON/civita/epsi,g5,g6,a5,a6,bb
       COMMON/extra/s1,s2,t1,t2
       COMMON/dotps/q1dq,q1dq2,w6
       data rho/.585d+00/
       IF(nup.GT.0) GO TO 1
       u1=1.
       u2=1.
       GO TO 3
1      IF(nup.GT.1) GO TO 2
       xt=1.-t1/.71
       xt=xt*xt
       xu=2.79/xt
       u1=xu*xu
       tau=t1/(4.*w1)
       u2=(1./(xt*xt)-xu*xu*tau)/(1.-tau)
       GO TO 3
2      x=t1/(t1-w3)
       en=w31-t1
       tau=t1/(4.*w1)
       rhot=rho-t1
       u1=(-.86926*rho*rho*w31/rhot/rhot-2.23422*w1*(1.-x)**4
     1 /(x*(x*.96-1.26)+1.))/t1
       u2=(-tau*u1-.12549*w31*t1*rho/rhot/rhot*w31*w31/en/en/w1)
     1 /(1.-en*en/(4.*w1*t1))
3      IF(ndown.GT.0) GO TO 4
       v1=1.
       v2=1.
       GO TO 6
4      IF(ndown.GT.1) GO TO 5
       xt=1.-t2/.71
       xt=xt*xt
       xu=2.79/xt
       v1=xu*xu
       tau=t2/(4.*w2)
       v2=(1./(xt*xt)-xu*xu*tau)/(1.-tau)
       GO TO 6
5      x=t2/(t2-w5)
       en=w52-t2
       tau=t2/(4.*w2)
       rhot=rho-t2
       v1=(-.86926*rho*rho*w52/rhot/rhot-2.23422*w2*(1.-x)**4
     1 /(x*(x*.96-1.26)+1.))/t2
       v2=(-tau*v1-.12549*w52*t2*rho/rhot/rhot*w52*w52/en/en/w2)
     1 /(1.-en*en/(4.*w2*t2))
6      CONTINUE
       qqq=q1dq*q1dq
       qdq=4.*w6-w4
       t22=512.*(bb*(delta*delta-gram)-(epsi-delta*(qdq+q1dq2))**2
     a -a1*a6*a6-a2*a5*a5-a1*a2*qqq)
       t12=128.*(-bb*(d2+g6)-2.*(t1+2.*w6)*(a2*qqq+a6*a6))*t1
       t21=128.*(-bb*(d4+g5)-2.*(t2+2.*w6)*(a1*qqq+a5*a5))*t2
       t11=64.*(bb*(qqq-g4-qdq*(t1+t2+2.*w6))-2.*(t1+2.*w6)*(t2+2.*w6)
     a *qqq)*t1*t2
       peripp=(((u1*v1*t11+u2*v1*t21+u1*v2*t12+u2*v2*t22)/(t1*t2*bb))
     a /(t1*t2*bb))*.25
       RETURN
       END
