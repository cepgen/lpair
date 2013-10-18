
       PROGRAM eemumu
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       DOUBLE PRECISION me,mu,f
       COMMON/inpu/me,mu,ebeam,const,sq
       COMMON/tell/nn
       COMMON/ini/xxx,yyy
       COMMON/outp/nout
       COMMON/cuts/angcut,encut,etacut
*
       INTEGER ndim,npoin,nprin,ntreat,nevent
*     
*     --- LPAIR data COMMON block
*     
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*     
*     --- HEPEVT COMMON block
       PARAMETER (NMXHEP=10000)
       COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &      JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &      VHEP(4,NMXHEP)
       REAL PHEP,VHEP
*     
       EXTERNAL f
*     
       CALL fileini(IPAR,LPAR)  ! initialize variables
*     CALL pawini  ! initialize histograms
       CALL genzini             ! initialize genz structures
*     
       pi = dacos(-1.d+00)
       nout=6
       
       IF(ipar(4).EQ.1) THEN
          OPEN (12,file='ineemm.data',status='old')
          REWIND 12
       ENDIF
       
       IF(ipar(17).EQ.1) THEN
          OPEN(20,file='events.ascii',status='new')
       ENDIF
       
       nn=0
*     
*     ---- particle masses (GeV)
       me = lpar(1)             ! incoming particle
       mu = lpar(2)             ! outgoing particle
*     
*     ---- beam energy (GeV)
       ebeam = lpar(3)
       sq=2.*ebeam
*     
*     ---- angle cuts
       angcut = dcos(pi*lpar(4))
*     
*     ---- rapidity cuts
       etacut =  lpar(5)
*     
*     ---- energy cuts
       encut = lpar(6)
*     
*     ---- constant (convert GeV**2 to picobarns)
       const=(19.732d+03)**2
*     
*     ---- VEGAS integration
*     
       xx = ranf(1211)
       OPEN (15,file='dl2.vegas.grid',status='new')
*     CALL vegas(fxn,bcc,ndim,nCALL,itmx,nprn,igraph)
       CALL vegas(f,0.1d-03,ipar(5),ipar(6),ipar(7),+1,ipar(4))
       CALL save(7,15)
       PRINT*,'Wrote VEGAS grid to dl2.vegas.grid'
       CLOSE (15)
*     
*     ---- cross-section calculation only: exit PROGRAM
*     
       IF(ipar(1).EQ.1) THEN
          PRINT*,'IPAR(1) = 1: cross-section calculation complete'
          stop
       ENDIF       
*     
*     ---- SETGEN: find local min/max
*     
       ndim   = ipar(5)
       npoin  = ipar(14)
       nbin   = ipar(15)
       nprin  = 1
       ntreat = ipar(13)
*     
       OPEN (15,file='dl2.vegas.grid',status='old')
       CALL restr(7,15)
       PRINT*,'Read dl2.vegas.grid'
       CLOSE (15)
       OPEN (16,file='dl2.lattice.1',status='new')
*     CALL setgen(f,ndim,npoin,nbin,nprin,ntreat)
       CALL setgen(f,ndim,npoin,nbin,nprin,ntreat)
       PRINT*,'setgen complete'
       CALL save2(7,16)
       PRINT*,'Wrote SETGEN maxima to dl2.lattice.1'
       CLOSE (16)
*     
*     ---- GENERA: generate some events
*     
       ndim   = ipar(5)
       nevent = ipar(12)
       nstrat = 0
       ntreat = ipar(13)
*     
       OPEN (15,file='dl2.vegas.grid',status='old')
       CALL restr(7,15)
       PRINT*,'Read dl2.vegas.grid'
       CLOSE (15)
       OPEN (16,file='dl2.lattice.1',status='old')
       CALL restr2(7,16)
       CLOSE (16)
       PRINT*,'Read maxima from dl2.lattice.1'
       OPEN (17,file='dl2.lattice.2',status='new')
       xxxx = ranf(1236785)
*     CALL genera(f,ndim,nevent,nstrat,ntreat)
       CALL genera(f,ndim,nevent,nstrat,ntreat)
       CALL save2(7,17)
       PRINT*,'Wrote new maxima to dl2.lattice.2'
       CLOSE (17)
*     
       CLOSE(20)                ! CLOSE events.ascii
*     
*     CALL pawend  ! CLOSE histograms
       CALL genzend             ! CLOSE genz
*     
       PRINT*,'GENZ events generated: NEVHEP = ',NEVHEP
*     
       stop
       END

      DOUBLE PRECISION function f(x)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION me,mu
      COMMON/inpu/me,mu,ebeam,const,sq
      COMMON/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5
     1     ,st5,cp3,sp3,cp5,sp5
      COMMON/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
      COMMON/lplot/xl(10),v1(2),v2(2),av(10)
      COMMON/extra/s1,s2,t1,t2
      COMMON/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
      COMMON/levi/gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
      COMMON/civita/epsi,g5,g6,a5,a6,bb
      COMMON/dotps/q1dq,q1dq2,w6
      COMMON/tell/nn
      COMMON/cuts/angcut,encut,etacut
*     
*     --- LPAIR data COMMON block
*     
      INTEGER IPAR(20)
      REAL*8 LPAR(20)
      COMMON/DATAPAR/IPAR,LPAR
*     
      DIMENSION x(10)
      data pi/3.141592459d+00/
      nn=nn+1
      
*     ---- proton form factors
      
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
         CALL gamgam(ebeam,me,me,w13,w25,mu,mu,0.D+00,sq,dj,0,x,1)
      ELSEIF(ipar(5).EQ.8) THEN
*     inelastic-elastic
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
*     CALL pawfil1
         CALL genzfil
*     
         RETURN
 20      PRINT *,'Matrix element is negative'
 30      f = 0.
         do 3 i=1,7
 3          xl(i)=-100.
         v1(1) = 0.
         v2(1) = 0.
         RETURN
      END
      
      
      
      SUBROUTINE orient(s,v1,v2,v3,v4,v5,dj,nopt,y)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      COMMON/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5,
     a     st5,cp3,sp3,cp5,sp5/variac/al3,al4,be4,be5,de3,de5,
     a     pp3,pp4,pp5
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
      e4=de3+de5
      e5=e2-de5
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
      
*----------------------------------------------------------------------
      
      SUBROUTINE pickin(s,v1,v2,v3,v4,v5,dj,nopt,y)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      DIMENSION y(4)
      COMMON/pickzz/w1,w2,w3,w4,w5,d1,d2,d5,d7,sl1
      COMMON/extra/s1,s2,t1,t2/accura/acc3,acc4
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

*----------------------------------------------------------------------

      SUBROUTINE maps2(s2,x,smin,smax,ds)
      IMPLICIT DOUBLE PRECISION (a-z)
      y=smax/smin
      s2=smin*y**x
      ds=s2*dlog(y)
      RETURN
      END

*----------------------------------------------------------------------
      
      SUBROUTINE mapla(x,y,z,u,xm,xp,d)
      IMPLICIT DOUBLE PRECISION (a-z)
      xmb=xm-y-z
      xpb=xp-y-z
      c=-4.*y*z
      alp=dsqrt(xpb*xpb+c)
      alm=dsqrt(xmb*xmb+c)
      am=xmb+alm
      ap=xpb+alp
      yy=ap/am
      zz=yy**u
      x=y+z+(am*zz-c/(am*zz))*0.5
      ax=dsqrt((x-y-z)**2+c)
      d=ax*dlog(yy)
      RETURN
      END
      
*----------------------------------------------------------------------

      SUBROUTINE inplot(now,ff,pdx)
      IMPLICIT DOUBLE PRECISION ( a-h,o-z)
      COMMON/lplot/xl(10),v1(2),v2(2),av(10)
      DIMENSION zav(10),yav(10),zsv(10),ysv(10),ztv(10)
      DIMENSION xlmax(10),xlmin(10),nlp(10),ltop(10),text(8,10),ll(10)
      DIMENSION numb(12)
      DIMENSION xls(42,10),yls(42,10),nlsn(42,10),mlsn(42,10),dls(10)
     1,xlav(10),xlsq(10),xlava(10),sxa(10),tlim(6),top(10),xltq(10)
      DIMENSION nbin(41),nlog(41),slog(41),tlog(41),hv(12)
      DIMENSION v1max(2),v1min(2),v2max(2),v2min(2),nv1(2)
     1,nv2(2),vtext(6,2)
      DIMENSION vm(12,12,2),nvm(12,12,2),bin1(2),bin2(2),vol(2)
     1,wm(12,12,2),mvm(12,12,2)
      COMMON/result/y,si,u,v
      character*1 hmin,hplus,hblank,hstar,char(40)
      save
      data tlim/1.6d+00,2.5d+00,4.0d+00,6.666666667d+00,10.d+00,
     +16.d+00/
      data hmin,hplus,hblank,hstar/'-','+',' ','*'/
      data mls,mav,ndmax/10,10,2/
      data ngraph/0/
      igraph=now
      now=0
      kk=0
      itt=0
      IF(igraph.le.0) GO TO 800
      IF(igraph.ne.ngraph)read (12,810)nls
      IF(igraph.ne.ngraph)PRINT 814,nls
      IF(nls.LT.0) nls=0
      IF(nls.EQ.0) GO TO 802
      IF(nls.GT.mls) GO TO 807
      IF(igraph.ne.ngraph)PRINT 815
      do 801 i=1,nls
      IF(igraph.ne.ngraph)
     1read (12,811)xlmin(i),xlmax(i),nlp(i),ltop(i),ll(i),
     2(text(j,i),j=1,8)
      IF(igraph.ne.ngraph) PRINT 816,i,xlmin(i),xlmax(i)
     1,nlp(i),ltop(i),ll(i),(text(j,i),j=1,8)
      IF(nlp(i).LT.1)nlp(i)=1
      IF(nlp(i).GT.40)nlp(i)=40
      dls(i)=(xlmax(i)-xlmin(i))/nlp(i)
      nlps=nlp(i)+2
      do 300 j=1,nlps
      yls(j,i)=0
300   mlsn(j,i)=0
801   CONTINUE
802   IF(igraph.ne.ngraph) read (12,810)ndd
      IF(igraph.ne.ngraph) PRINT 817,ndd
      IF(ndd.LT.0) ndd=0
      IF(ndd.EQ.0) GO TO 804
      IF(ndd.GT.ndmax) GO TO 807
      IF(igraph.ne.ngraph) PRINT 818
      do 803 i=1,ndd
      IF(igraph.ne.ngraph)
     1read (12,812)v1min(i),v1max(i),nv1(i),v2min(i),v2max(i),nv2(i)
     1,(vtext(j,i),j=1,6)
      IF(igraph.ne.ngraph) PRINT 819,i,v1min(i),v1max(i),nv1(i)
     1,v2min(i),v2max(i),nv2(i),(vtext(j,i),j=1,6)
      IF(nv1(i).LT.1)nv1(i)=1
      IF(nv2(i).LT.1)nv2(i)=1
      IF(nv1(i).GT.10)nv1(i)=10
      IF(nv2(i).GT.10)nv2(i)=10
      bin1(i)=(v1max(i)-v1min(i))/nv1(i)
      bin2(i)=(v2max(i)-v2min(i))/nv2(i)
      vol(i)=bin1(i)*bin2(i)
803   CONTINUE
      wtow=0.
      do 805 i=1,ndd
      do 805 j=1,12
      do 805 k=1,12
      wm(k,j,i)=0.
805   mvm(k,j,i)=0
804   CONTINUE
      IF(igraph.ne.ngraph)read (12,810)nave
      IF(igraph.ne.ngraph)PRINT 820,nave
      IF(nave.LT.0)nave=0
      IF(nave.GT.mav) GO TO 807
      do 11 i=1,mav
      yav(i)=0.
11    ysv(i)=0.
      kt=0
      GO TO 808
800   nave=0
      nls=0
      ndd=0
      GO TO 808
807   PRINT 813,mav,mls,ndmax
      stop
808   ngraph=igraph
      RETURN
      entry replot(now,ff,pdx)
      IF(nave.EQ.0) GO TO 49
      do 62 i=1,nave
      zav(i)=0.
      ztv(i)=0.
62    zsv(i)=0.
49    fsqa=0.
      kt=kt+1
      IF(nls.EQ.0) GO TO 303
      do 302 i=1,nls
      nlps=nlp(i)+2
      xlav(i)=0
      xltq(i)=0.
      xlsq(i)=0
      do 302 j=1,nlps
      xls(j,i)=0
302   nlsn(j,i)=0
303   CONTINUE
      IF(ndd.EQ.0) GO TO 403
      do 402 i=1,ndd
      n1=nv1(i)+2
      n2=nv2(i)+2
      do 402 i1=1,n1
      do 402 i2=1,n2
      vm(i1,i2,i)=0
402   nvm(i1,i2,i)=0
403   CONTINUE
      RETURN
      entry xplot(now,ff,pdx)
      fsqa=fsqa+ff*ff/pdx
      itt=itt+1
      IF(nls.EQ.0) GO TO 305
      do 304  i=1,nls
      nlps=(xl(i)-xlmin(i))/dls(i)+1.
      IF(nlps.LT.0)nlps=0
      IF(nlps.GT.nlp(i))nlps=nlp(i)+1
      nlps=nlps+1
      xls(nlps,i)=xls(nlps,i)+ff/dls(i)
      nlsn(nlps,i)=nlsn(nlps,i)+1
      xlav(i)=xlav(i)+ff*xl(i)
      xltq(i)=xltq(i)+ff*ff*xl(i)/pdx
304   xlsq(i)=xlsq(i)+(ff*xl(i))**2/pdx
305   CONTINUE
      IF(ndd.EQ.0) GO TO 405
      do 404 i=1,ndd
      i1=(v1(i)-v1min(i))/bin1(i)+2
      IF(i1.LT.1) i1=1
      IF(i1.GT.nv1(i)+2) i1=nv1(i)+2
      i2=(v2(i)-v2min(i))/bin2(i)+2
      IF(i2.LT.1) i2=1
      IF(i2.GT.nv2(i)+2) i2=nv2(i)+2
      vm(i1,i2,i)=vm(i1,i2,i)+ff/vol(i)
404   nvm(i1,i2,i)=nvm(i1,i2,i)+1
405   CONTINUE
      IF(nave.EQ.0) GO TO 99
      do 22 i=1,nave
      zav(i)=zav(i)+av(i)*ff
      ztv(i)=ztv(i)+ff*ff*av(i)/pdx
22    zsv(i)=zsv(i)+(av(i)*ff)**2/pdx
99    RETURN
      entry plotit(now,ff,pdx)
      IF(nls.EQ.0) GO TO 315
      IF(kk.GT.0) GO TO 307
      do 306 i=1,nls
      nlps=nlp(i)+2
      do 306 j=1,nlps
      mlsn(j,i)=nlsn(j,i)
306   yls(j,i)=xls(j,i)
      GO TO 310
307   vbef=vtot
      vu=(v/u)**2
      do 309 i=1,nls
      nlps=nlp(i)+2
      do 309 j=1,nlps
      IF(nlsn(j,i).EQ.0) GO TO 309
      IF(mlsn(j,i).EQ.0) GO TO 308
      al1=vu/nlsn(j,i)
      al2=vbef/mlsn(j,i)
      mlsn(j,i)=mlsn(j,i)+nlsn(j,i)
      yls(j,i)=(al2*xls(j,i)+al1*yls(j,i))/(al1+al2)
      GO TO 309
308   mlsn(j,i)=nlsn(j,i)
      yls(j,i)=xls(j,i)
309   CONTINUE
310   CONTINUE
      do 311 i=1,nls
      sxf=xlsq(i)-xlav(i)*xlav(i)
      sxt=xltq(i)-xlav(i)*u
      sx2=xlsq(i)/xlav(i)**2+fsqa/u**2-2.*xltq(i)/(xlav(i)*u)
      sx2=sx2*(xlav(i)/u)**2
      IF(kt.ne.1) GO TO 312
      xlava(i)=xlav(i)/u
      sxa(i)=sx2
      GO TO 311
312   xhelp=sx2+sxa(i)
      IF(xhelp.EQ.0) GO TO 311
      xlava(i)=(xlav(i)*sxa(i)/u+xlava(i)*sx2)/xhelp
      sxa(i)=sxa(i)*sx2/xhelp
311   CONTINUE
      vtot=(si/y)**2
      IF(now.ne.2) GO TO 315
      do 341 i=1,nls
      top(i)=0.
      nlps=nlp(i)+1
      do 341 j=2,nlps
      xls(j,i)=yls(j,i)/y
      IF(xls(j,i).GT.top(i))top(i)=xls(j,i)
341   CONTINUE
      do 342 i=1,nls
      IF(ltop(i).le.0)ltop(i)=i
      lto=ltop(i)
      IF(top(i).GT.top(lto))top(lto)=top(i)
342   CONTINUE
      ylog=0.5*dlog10(y*y)
      do 314 i=1,nls
      PRINT 321,i
      nlps=nlp(i)+1
      lto=ltop(i)
      top(i)=top(lto)
      IF(top(i).EQ.0)top(i)=1.
      an1=dlog10(top(i))
      n1=an1
      IF(n1.GT.an1)n1=n1-1
      z1=top(i)*10.**(-n1)
      do 343 l=1,4
      IF(z1.LT.tlim(l)) GO TO 344
343   CONTINUE
      l=5
344   IF(top(i).LT.1.6/(xlmax(i)-xlmin(i)))l=l+1
      topm=tlim(l)*10.**n1
      do 345 j=2,nlps
      nbin(j)=xls(j,i)*40./topm+1.5
      IF(ll(i).LT.0)nbin(j)=0
      IF(xls(j,i).le.0) GO TO 346
      tlog(j)=dlog10(xls(j,i))
      slog(j)=tlog(j)+ylog
      nlog(j)=(tlog(j)-n1)*8.+33.5
      IF(ll(i).GT.0)nlog(j)=0
      GO TO 345
346   slog(j)=0
      tlog(j)=0
      nlog(j)=0
345   CONTINUE
      PRINT 322,(text(j,i),j=1,8)
      n1p1=n1+1
      n1m4=n1-4
      PRINT 323,tlim(l),n1,n1p1,n1m4
      do 347 l=1,40
      char(l)=hmin
      IF(nlog(l+1).EQ.41)char(l)=hplus
      IF(nbin(l+1).EQ.41)char(l)=hstar
347   CONTINUE
      xmin=xlmin(i)
      xmax=xmin+dls(i)
      PRINT 324,xmin,xmax,yls(2,i),slog(2),xls(2,i),tlog(2),mlsn(2,i)
     1,char
***********************
*      CALL pawfil2(   i,real(xmin),real(xmax),real(xls(2,i)))
*      CALL pawfil2(10+i,real(xmin),real(xmax),real(yls(2,i)))
***********************
      do 348 j=3,nlps
      xmin=xmax
      xmax=xmin+dls(i)
      do 349 l=1,40
      char(l)=hblank
      IF(nlog(l+1).EQ.43-j)char(l)=hplus
      IF(nbin(l+1).EQ.43-j)char(l)=hstar
349   CONTINUE
      PRINT 324,xmin,xmax,yls(j,i),slog(j),xls(j,i),tlog(j),mlsn(j,i)
     1,char
***********************
*      CALL pawfil2(   i,real(xmin),real(xmax),real(xls(j,i)))
*      CALL pawfil2(10+i,real(xmin),real(xmax),real(yls(j,i)))
***********************
348   CONTINUE
      nlps1=nlps+1
      IF(nlps.EQ.41) GO TO 352
      do 351 j=nlps1,41
      do 350 l=1,40
      char(l)=hblank
      IF(nlog(l+1).EQ.43-j)char(l)=hplus
      IF(nbin(l+1).EQ.43-j)char(l)=hstar
350   CONTINUE
351   PRINT 325,char
352   do 353 l=1,40
      char(l)=hmin
      IF(nlog(l+1).EQ.1)char(l)=hplus
      IF(nbin(l+1).EQ.1)char(l)=hstar
353   CONTINUE
      PRINT 326,char
      el1=yls(1,i)*dls(i)
      el2=el1/y
      PRINT 327,el1,el2,mlsn(1,i)
      el1=yls(42,i)*dls(i)
      el2=el1/y
      PRINT 328,el1,el2,mlsn(nlps1,i)
      sxsq=dsqrt(sxa(i)/itt)
      PRINT 329,xlava(i),sxsq
314   CONTINUE
315   CONTINUE
      IF(ndd.EQ.0) GO TO 409
      wbef=wtot
      do 500 i=1,ndd
      nx=nv1(i)+2
      ny=nv2(i)+2
      IF(kk.GT.0) GO TO 502
      do 501 j=1,nx
      do 501 k=1,ny
      wm(j,k,i)=vm(j,k,i)
501   mvm(j,k,i)=nvm(j,k,i)
      GO TO 500
502   vu=(v/u)**2
      do 503 j=1,nx
      do 503 k=1,ny
      IF(nvm(j,k,i).EQ.0) GO TO 503
      IF(mvm(j,k,i).EQ.0) GO TO 504
      al1=vu/nvm(j,k,i)
      al2=vbef/mvm(j,k,i)
      mvm(j,k,i)=mvm(j,k,i)+nvm(j,k,i)
      wm(j,k,i)=(al2*vm(j,k,i)+al1*wm(j,k,i))/(al1+al2)
      GO TO 503
504   mvm(j,k,i)=nvm(j,k,i)
      wm(j,k,i)=vm(j,k,i)
503   CONTINUE
500   CONTINUE
      wtot=(si/y)**2
      IF(now.ne.2) GO TO 409
      do 408 i=1,ndd
      PRINT 481,i,(vtext(j,i),j=1,6)
      vvv=v2max(i)
      mvv=nv1(i)+2
      nvv=nv2(i)+1
      size=vol(i)/y
      do 406 i2=1,nvv
      j2=nvv+2-i2
      do 410 i1=1,mvv
410   numb(i1)=1000.*wm(i1,j2,i)*size+.5
      PRINT 486,(numb(i1),i1=1,mvv)
      PRINT 483,(wm(i1,j2,i),i1=1,mvv)
      PRINT 486,(mvm(i1,j2,i),i1=1,mvv)
      PRINT 484,vvv
      vvv=vvv-bin2(i)
      IF(dabs(vvv/bin2(i)).LT.1.e-10)vvv=0.
406   CONTINUE
      do 411 i1=1,mvv
411   numb(i1)=1000.*wm(i1,1,i)*size+.5
      PRINT 486,(numb(i1),i1=1,mvv)
      PRINT 483,(wm(i1,1,i),i1=1,mvv)
      PRINT 486,(mvm(i1,1,i),i1=1,mvv)
      PRINT 482
      mvv=mvv-1
      do 407 i1=1,mvv
      hv(i1)=v1min(i)+(i1-1)*bin1(i)
      IF(dabs(hv(i1)/bin1(i)).LT.1.d-10)hv(i1)=0.
407   CONTINUE
      PRINT 485,(hv(i1),i1=1,mvv)
408   CONTINUE
409   CONTINUE
      IF(nave.EQ.0) GO TO 23
      IF(now.EQ.2) PRINT 26
      do 24 i=1,nave
      sxf=zsv(i)-zav(i)*zav(i)
      sxt=zsv(i)/zav(i)**2+fsqa/u**2-2.*ztv(i)/(zav(i)*u)
      sx2=sxt*(zav(i)/u)**2
      IF(kt.ne.1) GO TO 21
      yav(i)=zav(i)/u
      ysv(i)=sx2
      GO TO 30
21    xhelp=sx2+ysv(i)
      IF(xhelp.EQ.0) GO TO 30
      yav(i)=(ysv(i)*zav(i)/u+yav(i)*sx2)/xhelp
      ysv(i)=ysv(i)*sx2/xhelp
30    yssq=dsqrt(ysv(i)/itt)
      IF(now.EQ.2)PRINT 27,i,yav(i),yssq
24    CONTINUE
23    now=1
      kk=kk+1
      RETURN
27    format(12x,i2,9x,d15.5,5x,d15.3)
26    format(1h1,10x,46hthe following are averages with error estimate/)
321   format(1h1,40x,40hsingle dIFferential cross-section number,i3///)
322   format(38h single dIFferential cross section of ,8a4/)
323   format(11x,6hlimits,9x,1hi,16x,19haccumulated results,15x,1hi,24x
     1,9hupper bin,6x,9hlower bin/26x,1hi,50x,20hi * linear      plot,
     2f8.2,5h*10**,i3,8x,1h0/5x,5hlower,7x,5hupper,4x,1hi,5x,5hds/dx,4x
     3,56halog10   (ds/dx)/s  alog10  points  i + logarithmic plot,6x
     4,4h10**,i3,8x,4h10**,i3/2x,24(1h-),1hi,50(1h-),1hi)
324   format(d12.4,d12.4,3h  i,2(d12.4,f8.2),i8,3h  i,4x,1hi,40a1,1hi)
325   format(26x,1hi,50x,1hi,4x,1hi,40a1,1hi)
326   format(2x,24(1h-),1hi,50(1h-),1hi,4x,1hi,40a1,1hi)
327   format(7x,15htotal underflow,4x,1hi,d12.4,d20.4,i16,2x,1hi)
328   format(7x,15htotal  overflow,4x,1hi,d12.4,d20.4,i16,2x,1hi)
329   format(//19x,21haccumulated average =,d12.5
     1/19x,21hestimated error     =,d12.5)
481   format(1h1,45x,40hdouble dIFferential cross-section number,i3/
     1/60x,7hx-axis ,3a4/60x,7hy-axis ,3a4/)
482   format(20x,11(1hi,9x))
483   format(11x,d9.3,11(1hi,d9.3))
484   format(1x,d10.3,  9h---------,11(10hi---------))
485   format(1h0,14x,11d10.3)
486   format(11x,i8,1x,11(1hi,i8,1x))
810   format(i2)
811   format(2d12.4,3i2,8a4)
812   format(2(2d10.3,i4),6a4)
813   format(12h1***error***,10x,24htoo many plots requested//
     122x,20hthe upper limits are //19x,i2,9h averages//19x,i2,22h one d
     2mensional plots//19x,i2,22h two DIMENSIONal plots////22x,25h***exe
     3cution is halted***   )
814   format(37h1number of single dIFferential cross
     1,20hsections requested =,i3/)
815   format(30h information on the data cards//
     13h  i,10x,5hxlmin,12x,5hxlmax,7x,
     224hbins  correllation  type,19x,4htext/)
816   format(i3,2e17.4,i8,i9,i10,5x,8a4)
817   format(37h0number of double dIFferential cross
     1,20hsections requested =,i3/)
818   format(30h information on the data cards//
     13h  i,8x,5hv1min,10x,5hv1max,4x,4hbins,8x,5hv2min
     2,10x,5hv2max,4x,4hbins,8x,6htext 1,8x,6htext 2/)
819   format(i3,2d15.3,i5,1x,2d15.3,i5,6x,3a4,2x,3a4)
820   format(31h0number of averages requested =,i3)
      END

*----------------------------------------------------------------------

       DOUBLE PRECISION FUNCTION peripp(nup,ndown)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       COMMON/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
       COMMON/levi/gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
       COMMON/civita/epsi,g5,g6,a5,a6,bb
       COMMON/extra/s1,s2,t1,t2/dotps/q1dq,q1dq2,w6
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

*----------------------------------------------------------------------

       SUBROUTINE gamgam(ebeam,v1,v2,v3,v5,v6,v7,vmin,vmax,dj,nopt,x,nm)
       IMPLICIT DOUBLE PRECISION (a-h,o-z)
       COMMON/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5,
     a st5,cp3,sp3,cp5,sp5/variac/al3,al4,be4,be5,de3,de5,pp3,pp4,pp5
       COMMON/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
       COMMON/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
       COMMON/dotps/q1dq,q1dq2,w6/extra/s1,s2,t1,t2
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
       RETURN
       END
      
      SUBROUTINE mapw2(w2,x,w2min,w2max,dw)
      IMPLICIT DOUBLE PRECISION (a-z)
*     wmin=1./w2min
*     wmax=1./w2max
*     dw=wmin-wmax
*     w2=1./(wmax+dw*x)
*     dw=dw*w2*w2
      y = w2max/w2min
      w2 = w2min*y**x
      dw = w2*dlog(y)
      RETURN
      END
      
      SUBROUTINE mapt1(t,x,tmin,tmax,dt)
      IMPLICIT DOUBLE PRECISION (a-z)
      y=tmax/tmin
      t=tmin*y**x
      dt=-t*dlog(y)
*     dt = tmin-tmax
*     t = tmin - x*dt
      RETURN
      END
      
      SUBROUTINE mapt2(t,x,tmin,tmax,dt)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      y=tmax/tmin
      t=tmin*y**x
      dt=-t*dlog(y)
*     dt = tmin-tmax
*     t = tmin - x*dt
      RETURN
      END
