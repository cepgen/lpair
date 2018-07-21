      DOUBLE PRECISION function f(x)
      IMPLICIT NONE

*---- Local variables
      DOUBLE PRECISION x,dj,pi
      DOUBLE PRECISION peripp
      DOUBLE PRECISION st3,wnvmass,yy
     1     ,w31min,w31max,w13,w132,dw13
     1     ,w52min,w52max,w25,w252,dw25
      DOUBLE PRECISION pt6,pz6,pt7,pz7,eta6,eta7
      INTEGER i

*---- Common blocks variables
      DOUBLE PRECISION me,mu,ebeam,const,sq
      DOUBLE PRECISION e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st4,ct4
     1     ,ct5,st5,cp3,sp3,cp5,sp5
      DOUBLE PRECISION e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
      DOUBLE PRECISION xl,v1,v2,av
      DOUBLE PRECISION s1,s2,t1,t2
      DOUBLE PRECISION w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
      DOUBLE PRECISION gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
      DOUBLE PRECISION epsi,g5,g6,a5,a6,bb
      DOUBLE PRECISION q1dq,q1dq2,w6
      INTEGER nn
      DOUBLE PRECISION angcut,encut,etacut,mxcut

      COMMON/inpu/me,mu,ebeam,const,sq
      COMMON/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4
     1     ,ct5,st5,cp3,sp3,cp5,sp5
      COMMON/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
      COMMON/lplot/xl(10),v1(2),v2(2),av(10)
      COMMON/extra/s1,s2,t1,t2
      COMMON/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
      COMMON/levi/gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
      COMMON/civita/epsi,g5,g6,a5,a6,bb
      COMMON/dotps/q1dq,q1dq2,w6
      COMMON/tell/nn
      COMMON/cuts/angcut,encut,etacut,mxcut
*     
      INTEGER IPAR(20)
      REAL*8 LPAR(20)
      COMMON/DATAPAR/IPAR,LPAR
*     
      DIMENSION x(10)
      data pi/3.141592459d+00/
      nn=nn+1
      
*---- proton form factors
      
      IF((ipar(5).EQ.8).OR.(ipar(5).EQ.9)) THEN
         w31min = (me+0.135)**2
         w31max = (sq-me-2.*mu)**2
         yy = w31max/w31min
         w132 = w31min*yy**x(8)
         w13 = dsqrt(w132)
         IF(ipar(10).EQ.0) THEN
            dw13 = w132*dlog(yy) ! Vermaseren (paper)
         ELSEIF(ipar(10).EQ.1) THEN
            dw13 = w132*dlog(yy)/(me*me) ! Vermaseren (experimental)
         ENDIF
         w52min = (me+0.135)**2
         w52max = (sq-w13-2*mu)**2
         yy = w52max/w52min
         w252 = w52min*yy**x(9)
         w25 = dsqrt(w252)
         IF(ipar(10).EQ.0) THEN
            dw25 = w252*dlog(yy) ! Vermaseren (paper)
         ELSEIF(ipar(10).EQ.1) THEN
            dw25 = w252*dlog(yy)/(me*me) ! Vermaseren (experimental)
         ENDIF
      ENDIF
     
      IF(ipar(5).EQ.9) THEN
*     inelastic-inelastic
         IF((w13.gt.mxcut).or.(w25.gt.mxcut)) GOTO 30
         CALL gamgam(ebeam,me,me,w13,w25,mu,mu,0.D+00,sq,dj,0,x,1)
      ELSEIF(ipar(5).EQ.8) THEN
*     inelastic-elastic
         IF(w13.gt.mxcut**2) GOTO 30
         CALL gamgam(ebeam,me,me,w13,me,mu,mu,0.D+00,sq,dj,0,x,1)
      ELSEIF(ipar(5).EQ.7) THEN
*     elastic-elastic
         CALL gamgam(ebeam,me,me,me,me,mu,mu,0.0D+02,sq,dj,0,x,1)
      ENDIF
*     
      IF(dj.EQ.0)GO TO 20
*     
      IF(ipar(9).EQ.1) THEN
*     
*     * -------------------------------- *
*     * ---- Bryan's analysis cuts ----- *
*     * -------------------------------- *
*     
*     1 ----- cos(theta) cuts
*     IF ( dabs(ct6) .GT. angcut ) goto 30
*     IF ( dabs(ct7) .GT. angcut ) goto 30
*     
*     2 ----- rapidity cuts
*     
         pt6 = p6*st6
         pt7 = p7*st7
         pz6 = p6*ct6
         pz7 = p7*ct7
*     
         eta6=dsign(dlog((dsqrt(pt6**2+pz6**2)+dabs(pz6))/pt6),pz6)
         eta7=dsign(dlog((dsqrt(pt7**2+pz7**2)+dabs(pz7))/pt7),pz7)
*     
         IF (dabs(eta6).GT.etacut) goto 30
         IF (dabs(eta7).GT.etacut) goto 30
*     
*     3 ----- transverse momentum cuts
*     
         IF ( ( p6*st6.LT.lpar(7) ).OR.( p7*st7.LT.lpar(7) ) ) goto 30
*     
*     4 ----- invariant mass cuts
*     
         wnvmass = dsqrt( (e6+e7)**2 - (p6*st6*cp6 + p7*st7*cp7)**2
     &        - (p6*st6*sp6 + p7*st7*sp7)**2
     &        - (    p6*ct6 +     p7*ct7)**2 )
         IF( (wnvmass.LT.lpar(8) ).OR.( wnvmass.GT.lpar(9) ) ) goto 30
*     
      ELSEIF(ipar(9).EQ.0) THEN
*     
*     * ----------------------------------- *
*     * ---- Vermaseren analysis cuts ----- *
*     * ----------------------------------- *
*     
*     1 --- invariant mass cut (UNUSED)
*     IF ( w4 .LT. 9. ) goto 30
*     
*     2 --- energy cuts (UNUSED)
*     IF ( e6 .LT. encut ) goto 30
*     IF ( e7 .LT. encut ) goto 30
*     
*     3 --- cos(theta) cuts (USED)
         IF ( dabs(ct6) .GT. angcut ) goto 30
         IF ( dabs(ct7) .GT. angcut ) goto 30
*     
*     4 --- transverse momentum cuts (UNUSED)
*     IF ( p6*st6.LT.0.4*dsqrt(w4) ) goto 30
*     IF ( p7*st7.LT.0.4*dsqrt(w4) ) goto 30
*     
*     5 --- transverse momentum cuts and cos(theta) cuts (USED)
         IF ( ( p6*st6.LT.1.0 ).AND.( dabs(ct6).LT.0.75 ) ) goto 30
         IF ( ( p7*st7.LT.1.0 ).AND.( dabs(ct7).LT.0.75 ) ) goto 30
*     
*     6 --- longitudinal momentum cuts and cos(theta) cuts (USED)
         IF ( ( abs(p6*ct6).LT.1.0 ).AND.( dabs(ct6).GT.0.75 ) ) goto 30
         IF ( ( abs(p7*ct7).LT.1.0 ).AND.( dabs(ct7).GT.0.75 ) ) goto 30
*     
*     7 --- tagging cut (UNUSED)
*     ppcut = 0.5*(p3*st3)**2 + 0.5*(p5*st5)**2 - 0.25*(p4*st4)**2 
*     IF ( ppcut .GT. 0.01 ) goto 30
*     
      ENDIF

*     
*     ------------------------------------ *
*     ---- matrix element calculation ---- *
*     ------------------------------------ *
*     
      IF(ipar(5).EQ.9) THEN
         f=const*dj*peripp(2,2)*dw13*dw25 ! inelastic-inelastic p p
      ELSEIF(ipar(5).EQ.8) THEN
         f=const*dj*peripp(2,1)*dw13 ! inelastic-elastic p p
      ELSEIF(ipar(5).EQ.7) THEN 
         IF(ipar(8).EQ.2212) THEN
            f=const*dj*peripp(1,1) ! elastic-elastic p p
         ELSEIF(ipar(8).EQ.11) THEN
            f=const*dj*peripp(0,0) ! electron-positron
         ENDIF
      ENDIF
*     
      IF(f.LT.0)GO TO 20 
*     
*     ---- fill entries for histograms
      do 2 i=1,7
 2       xl(i)=x(i)
         xl(1) = dlog10(-t1)
         xl(2) = dlog10(-t2)
         xl(3) = dsqrt(w4)
         xl(4) = p6*st6
         xl(5) = p7*st7
         xl(6) = p4*st4
         v1(1) = xl(1)
         v2(1) = xl(2)
*     
*     ---- fill ntuple entries and genz (move to accept.f)
*     
c         CALL pawfil1
c         CALL genzfil
*     
         RETURN
 20      CONTINUE
c     PRINT *,'Matrix element is negative'
 30      f = 0.
         do 3 i=1,7
 3          xl(i)=-100.
         v1(1) = 0.
         v2(1) = 0.
         RETURN
      END
