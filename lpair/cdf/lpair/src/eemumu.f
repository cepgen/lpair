
       program test6
       implicit double precision (a-h,o-z)
       double precision me,mu,f
       common/inpu/me,mu,ebeam,const,sq
       common/tell/nn
       common/ini/xxx,yyy
       common/outp/nout
       common/cuts/angcut,encut,etacut
*
       integer ndim,npoin,nprin,ntreat,nevent
*
* --- LPAIR data common block
*
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*
* --- HEPEVT common block
       PARAMETER (NMXHEP=10000)
       COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     & JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
       REAL PHEP,VHEP
*
       external f
*
       call fileini ! initialize variables
       call pawini  ! initialize histograms
       call genzini ! initialize genz structures
*
       pi = dacos(-1.d+00)
       nout=6

       if(ipar(4).eq.1) then
         open (12,file='ineemm.data',status='old')
         rewind 12
       endif

       if(ipar(17).eq.1) then
          open(20,file='events.ascii',status='new')
       endif

       nn=0
*
* ---- particle masses (GeV)
       me = lpar(1)  ! incoming particle
       mu = lpar(2)  ! outgoing particle
*
* ---- beam energy (GeV)
       ebeam = lpar(3)
       sq=2.*ebeam
*
* ---- angle cuts
       angcut = dcos(pi*lpar(4))
*
* ---- rapidity cuts
       etacut =  lpar(5)
*
* ---- energy cuts
       encut = lpar(6)
*
* ---- constant (convert GeV**2 to picobarns)
       const=(19.732d+03)**2
*
* ---- VEGAS integration
*
       xx = ranf(1211)
       open (15,file='dl2.vegas.grid',status='new')
*      call vegas(fxn,bcc,ndim,ncall,itmx,nprn,igraph)
       call vegas(f,0.1d-03,ipar(5),ipar(6),ipar(7),+1,ipar(4))
       call save(7,15)
       print*,'Wrote VEGAS grid to dl2.vegas.grid'
       close (15)
*
* ---- cross-section calculation only: exit program
*
       if(ipar(1).eq.1) then
         print*,'IPAR(1) = 1: cross-section calculation complete'
         stop
       endif       
*
* ---- SETGEN: find local min/max
*
       ndim   = ipar(5)
       npoin  = ipar(14)
       nbin   = ipar(15)
       nprin  = 1
       ntreat = ipar(13)
*
       open (15,file='dl2.vegas.grid',status='old')
       call restr(7,15)
       print*,'Read dl2.vegas.grid'
       close (15)
       open (16,file='dl2.lattice.1',status='new')
*      call setgen(f,ndim,npoin,nbin,nprin,ntreat)
       call setgen(f,ndim,npoin,nbin,nprin,ntreat)
       print*,'setgen complete'
       call save2(7,16)
       print*,'Wrote SETGEN maxima to dl2.lattice.1'
       close (16)
*
* ---- GENERA: generate some events
*
       ndim   = ipar(5)
       nevent = ipar(12)
       nstrat = 0
       ntreat = ipar(13)
*
       open (15,file='dl2.vegas.grid',status='old')
       call restr(7,15)
       print*,'Read dl2.vegas.grid'
       close (15)
       open (16,file='dl2.lattice.1',status='old')
       call restr2(7,16)
       close (16)
       print*,'Read maxima from dl2.lattice.1'
       open (17,file='dl2.lattice.2',status='new')
       xxxx = ranf(1236785)
*      call genera(f,ndim,nevent,nstrat,ntreat)
       call genera(f,ndim,nevent,nstrat,ntreat)
       call save2(7,17)
       print*,'Wrote new maxima to dl2.lattice.2'
       close (17)
*
       close(20)    ! close events.ascii
*
       call pawend  ! close histograms
       call genzend ! close genz
*
       print*,'GENZ events generated: NEVHEP = ',NEVHEP
*
       stop
       end

       double precision function f(x)
       implicit double precision (a-h,o-z)
       double precision me,mu
       common/inpu/me,mu,ebeam,const,sq
       common/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5
     1 ,st5,cp3,sp3,cp5,sp5
       common/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
       common/lplot/xl(10),v1(2),v2(2),av(10)
       common/extra/s1,s2,t1,t2
       common/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
       common/levi/gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
       common/civita/epsi,g5,g6,a5,a6,bb
       common/dotps/q1dq,q1dq2,w6
       common/tell/nn
       common/cuts/angcut,encut,etacut
*
* --- LPAIR data common block
*
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*
       dimension x(10)
       data pi/3.141592459d+00/
       nn=nn+1

* ---- proton form factors

       if((ipar(5).eq.8).or.(ipar(5).eq.9)) then
         w31min = (me+0.135)**2
         w31max = (sq-me-2.*mu)**2
         yy = w31max/w31min
         w132 = w31min*yy**x(8)
         w13 = dsqrt(w132)
         if(ipar(10).eq.0) then
           dw13 = w132*dlog(yy)         ! Vermaseren (paper)
         elseif(ipar(10).eq.1) then
           dw13 = w132*dlog(yy)/(me*me) ! Vermaseren (experimental)
         endif
         w52min = (me+0.135)**2
         w52max = (sq-w13-2*mu)**2
         yy = w52max/w52min
         w252 = w52min*yy**x(9)
         w25 = dsqrt(w252)
         if(ipar(10).eq.0) then
           dw25 = w252*dlog(yy)         ! Vermaseren (paper)
         elseif(ipar(10).eq.1) then
           dw25 = w252*dlog(yy)/(me*me) ! Vermaseren (experimental)
         endif
       endif

       if(ipar(5).eq.9) then
*        inelastic-inelastic
         call gamgam(ebeam,me,me,w13,w25,mu,mu,0.D+00,sq,dj,0,x,1)
       elseif(ipar(5).eq.8) then
*        inelastic-elastic
         call gamgam(ebeam,me,me,w13,me,mu,mu,0.D+00,sq,dj,0,x,1)
       elseif(ipar(5).eq.7) then
*        elastic-elastic
         call gamgam(ebeam,me,me,me,me,mu,mu,0.0D+02,sq,dj,0,x,1)
       endif
*
       if(dj.eq.0)go to 20
*
       if(ipar(9).eq.1) then
*
*        * -------------------------------- *
*        * ---- Bryan's analysis cuts ----- *
*        * -------------------------------- *
*
*        1 ----- cos(theta) cuts
*         if ( dabs(ct6) .gt. angcut ) goto 30
*         if ( dabs(ct7) .gt. angcut ) goto 30
*
*        2 ----- rapidity cuts
*
         pt6 = p6*st6
         pt7 = p7*st7
         pz6 = p6*ct6
         pz7 = p7*ct7
*
         eta6=dsign(dlog((dsqrt(pt6**2+pz6**2)+dabs(pz6))/pt6),pz6)
         eta7=dsign(dlog((dsqrt(pt7**2+pz7**2)+dabs(pz7))/pt7),pz7)
*
         if (dabs(eta6).gt.etacut) goto 30
         if (dabs(eta7).gt.etacut) goto 30
*
*        3 ----- transverse momentum cuts
*
         if ( ( p6*st6.lt.lpar(7) ).or.( p7*st7.lt.lpar(7) ) ) goto 30
*
*        4 ----- invariant mass cuts
*
         wnvmass = dsqrt( (e6+e7)**2 - (p6*st6*cp6 + p7*st7*cp7)**2
     &                               - (p6*st6*sp6 + p7*st7*sp7)**2
     &                               - (    p6*ct6 +     p7*ct7)**2 )
         if( (wnvmass.lt.lpar(8) ).or.( wnvmass.gt.lpar(9) ) ) goto 30
*
       elseif(ipar(9).eq.0) then
*
*        * ----------------------------------- *
*        * ---- Vermaseren analysis cuts ----- *
*        * ----------------------------------- *
*
*        1 --- invariant mass cut (UNUSED)
*         if ( w4 .lt. 9. ) goto 30
*
*        2 --- energy cuts (UNUSED)
*         if ( e6 .lt. encut ) goto 30
*         if ( e7 .lt. encut ) goto 30
*
*        3 --- cos(theta) cuts (USED)
         if ( dabs(ct6) .gt. angcut ) goto 30
         if ( dabs(ct7) .gt. angcut ) goto 30
*
*        4 --- transverse momentum cuts (UNUSED)
*         if ( p6*st6.lt.0.4*dsqrt(w4) ) goto 30
*         if ( p7*st7.lt.0.4*dsqrt(w4) ) goto 30
*
*        5 --- transverse momentum cuts and cos(theta) cuts (USED)
         if ( ( p6*st6.lt.1.0 ).and.( dabs(ct6).lt.0.75 ) ) goto 30
         if ( ( p7*st7.lt.1.0 ).and.( dabs(ct7).lt.0.75 ) ) goto 30
*
*        6 --- longitudinal momentum cuts and cos(theta) cuts (USED)
         if ( ( abs(p6*ct6).lt.1.0 ).and.( dabs(ct6).gt.0.75 ) ) goto 30
         if ( ( abs(p7*ct7).lt.1.0 ).and.( dabs(ct7).gt.0.75 ) ) goto 30
*
*        7 --- tagging cut (UNUSED)
*         ppcut = 0.5*(p3*st3)**2 + 0.5*(p5*st5)**2 - 0.25*(p4*st4)**2 
*         if ( ppcut .gt. 0.01 ) goto 30
*
       endif
*
* ------------------------------------ *
* ---- matrix element calculation ---- *
* ------------------------------------ *
*
       if(ipar(5).eq.9) then
         f=const*dj*peripp(2,2)*dw13*dw25   ! inelastic-inelastic p p
       elseif(ipar(5).eq.8) then
         f=const*dj*peripp(2,1)*dw13        ! inelastic-elastic p p
       elseif(ipar(5).eq.7) then 
         if(ipar(8).eq.2212) then
           f=const*dj*peripp(1,1)           ! elastic-elastic p p
         elseif(ipar(8).eq.11) then
           f=const*dj*peripp(0,0)           ! electron-positron
         endif
       endif
*
       if(f.lt.0)go to 20 
*
* ---- fill entries for histograms
       do 2 i=1,7
2      xl(i)=x(i)
       xl(1) = dlog10(-t1)
       xl(2) = dlog10(-t2)
       xl(3) = dsqrt(w4)
       xl(4) = p6*st6
       xl(5) = p7*st7
       xl(6) = p4*st4
       v1(1) = xl(1)
       v2(1) = xl(2)
*
* ---- fill ntuple entries and genz (move to accept.f)
*
*       call pawfil1
*       call genzfil
*
       return
20     print *,'Matrix element is negative'
30     f = 0.
       do 3 i=1,7
3      xl(i)=-100.
       v1(1) = 0.
       v2(1) = 0.
       return
       end



       subroutine orient(s,v1,v2,v3,v4,v5,dj,nopt,y)
       implicit double precision (a-h,o-z)
       common/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5,
     a st5,cp3,sp3,cp5,sp5/variac/al3,al4,be4,be5,de3,de5,pp3,pp4,pp5
       common/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
       common/extra/s1,s2,t1,t2
       common/levi/gram,dd1,dd2,dd3,dd4,dd5,delta,g4,sa1,sa2
       common/dotp/p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,q1,q2
       dimension y(4)
       call pickin(s,v1,v2,v3,v4,v5,dj,nopt,y)
       if(dj.eq.0)go to 10
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
       if(e4.lt.v4)go to 10
       p3=dsqrt(e3*e3-w3)
       p4=dsqrt((e4-v4)*(e4+v4))
       if(p4.eq.0)go to 10
       p5=dsqrt(e5*e5-w5)
       pp3=dsqrt(dd1/s)/p
       pp5=dsqrt(dd3/s)/p
       st3=pp3/p3
       st5=pp5/p5
       if(st3.gt.1..or.st5.gt.1.)go to 10
       ct3=dsqrt(1.-st3*st3)
       ct5=dsqrt(1.-st5*st5)
       if(e1*e3.lt.p13)ct3=-ct3
       if(e2*e5.gt.p25)ct5=-ct5
       al3=st3*st3/(1.+ct3)
       be5=st5*st5/(1.-ct5)
       if(dd5.lt.0)go to 10
       pp4=dsqrt(dd5/s)/p
       st4=pp4/p4
       if(st4.gt.1.)go to 10
       ct4=dsqrt(1.-st4*st4)
       if(e1*e4.lt.p14)ct4=-ct4
       al4=1.-ct4
       be4=1.+ct4
       if(ct4.lt.0)be4=st4*st4/al4
       if(ct4.ge.0)al4=st4*st4/be4
       rr=dsqrt(-gram/s)/(p*pp4)
       sp3=rr/pp3
       sp5=-rr/pp5
       if(dabs(sp3).gt.1..or.dabs(sp5).gt.1.)go to 10
       cp3=-dsqrt(1.-sp3*sp3)
       cp5=-dsqrt(1.-sp5*sp5)
       a1=pp3*cp3-pp5*cp5
       if(dabs(pp4+pp3*cp3+cp5*pp5).lt.dabs(dabs(a1)-pp4))go to 1
       if(a1.lt.0)cp5=-cp5
       if(a1.ge.0)cp3=-cp3
1      return
10     dj=0.
       return
       end

       subroutine pickin(s,v1,v2,v3,v4,v5,dj,nopt,y)
       implicit double precision (a-h,o-z)
       dimension y(4)
       common/pickzz/w1,w2,w3,w4,w5,d1,d2,d5,d7,sl1
       common/extra/s1,s2,t1,t2/accura/acc3,acc4
       common/levi/gram,dd1,dd2,dd3,dd4,dd5,delta,g4,sa1,sa2
       common/dotp/p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,p1k2,p2k1
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
       if(rl1.le.0)go to 20
       sl1=dsqrt(rl1)
       if(nopt.ne.0)go to 1
       smax=s+w3-2.*v3*dsqrt(s)
       call maps2(s2,x3,sig1,smax,ds2)
       sig1=s2
1      sp=s+w3-sig1
       d3=sig1-w2
       rl2=sp*sp-4.*s*w3
       if(rl2.le.0)go to 20
       sl2=sqrt(rl2)
       t1max=w1+w3-(ss*sp+sl1*sl2)/(2.*s)
       t1min=(d1*d3+(d3-d1)*(d3*w1-d1*w2)/s)/t1max
       call mapt1(t1,x1,t1min,t1max,dt1)
       d4=w4-t1
       d8=t1-w2
       t13=t1-w1-w3
       sa1=-(t1-d1)*(t1-d1)*0.25+w1*t1
       if(sa1.ge.0)go to 20
       sl3=dsqrt(-sa1)
       if(w1.eq.0)go to 3
       sb=(s*(t1-d1)+d5*t13)/(2.*w1)+w3
       sd=sl1*sl3/w1
       se=(s*(t1*(s+t13-w2)-w2*d1)+w3*(d5*d8+w2*w3))/w1
       if(dabs((sb-sd)/sd).lt.1.0)go to 2
       splus=sb-sd
       s2max=se/splus
       go to 4
2      s2max=sb+sd
       splus=se/s2max
       go to 4
3      s2max=(s*(t1*(s+d8-w3)-w2*w3)+w2*w3*(w2+w3-t1))/ss/t13
       splus=sig2
4      s2x=s2max
       if(nopt)5,6,7
5      if(splus.gt.sig2)sig2=splus
       if(nopt.lt.-1)call maps2(s2,x3,sig2,s2max,ds2)
       if(nopt.eq.-1)call mapla(s2,t1,w2,x3,sig2,s2max,ds2)
6      s2x=s2
7      r1=s2x-d8
       r2=s2x-d6
       rl4=(r1*r1-4.*w2*s2x)*(r2*r2-4.*w5*s2x)
       if(rl4.le.0)go to 20
       sl4=dsqrt(rl4)
       t2max=w2+w5-(r1*r2+sl4)/(2.*s2x)
       t2min=(d2*d4+(d4-d2)*(d4*w2-d2*t1)/s2x)/t2max
       call mapt2(t2,x2,t2min,t2max,dt2)
       d7=t1-t2
       r3=d4-t2
       r4=d2-t2
       b=r3*r4-2.*(t1+w2)*t2
       c=t2*d6*d8+(d6-d8)*(d6*w2-d8*w5)
       t25=t2-w2-w5
       sa2=-r4*r4*0.25+w2*t2
       if(sa2.ge.0)go to 20
       sl6=2.*dsqrt(-sa2)
       g4=-0.25*r3*r3+t1*t2
       if(g4.ge.0)go to 20
       sl7=sqrt(-g4)*2.
       sl5=sl6*sl7
       if(dabs((sl5-b)/sl5).lt.1.0)go to 8
       s2p=(sl5-b)/(2.*t2)
       s2min=c/(t2*s2p)
       go to 9
8      s2min=(-sl5-b)/(2.*t2)
       s2p=c/(t2*s2min)
9      if(nopt.gt.1)call maps2(s2,x3,s2min,s2max,ds2)
       if(nopt.eq.1)call mapla(s2,t1,w2,x3,s2min,s2max,ds2)
       ap=-(s2+d8)*(s2+d8)*0.25+s2*t1
       if(w1.eq.0)go to 10
       dd1=-w1*(s2-s2max)*(s2-splus)*0.25
       go to 11
10     dd1=ss*t13*(s2-s2max)*0.25
11     dd2=-t2*(s2-s2p)*(s2-s2min)*0.25
       yy4=dcos(pi*y(4))
       dd=dd1*dd2
       p12=0.5*(s-w1-w2)
       st=s2-t1-w2
       delb=(2.*w2*r3+r4*st)*(4.*p12*t1-(t1-d1)*st)/(16.*ap)
       if(dd.le.0)go to 20
       delta=delb-yy4*st*dsqrt(dd)/(2.*ap)
       s1=t2+w1+(2.*p12*r3-4.*delta)/st
       if(ap.ge.0)go to 20
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
       if(w2.eq.0)go to 14
       sbb=(s*(t2-d2)-d5*t25)/(2.*w2)+w5
       sdd=sl1*sl6/(2.*w2)
       see=(s*(t2*(s+t25-w1)-w1*d2)+w5*(w1*w5-d5*(t2-w1)))/w2
       if(sbb/sdd.lt.0)go to 12
       s1p=sbb+sdd
       s1m=see/s1p
       go to 13
12     s1m=sbb-sdd
       s1p=see/s1m
13     dd3=-w2*(s1p-s1)*(s1m-s1)*0.25
       go to 15
14     s1p=(s*(t2*(s-w5+t2-w1)-w1*w5)+w1*w5*(w1+w5-t2))/t25/(s-d5)
       dd3=-t25*(s-d5)*(s1p-s1)*0.25
15     acc3=(s1p-s1)/(s1p+s1)
       ssb=t2+w1-r3*(d1-t1)*0.5/t1
       ssd=sl3*sl7/t1
       sse=(t2-w1)*(w4-w3)+(t2-w4+d1)*((t2-w1)*w3-(w4-w3)*w1)/t1
       if(ssb/ssd.lt.0)go to 16
       s1pp=ssb+ssd
       s1pm=sse/s1pp
       go to 17
16     s1pm=ssb-ssd
       s1pp=sse/s1pm
17     dd4=-t1*(s1-s1pp)*(s1-s1pm)*0.25
       acc4=(s1-s1pm)/(s1+s1pm)
       dd5=dd1+dd3+((p12*(t1-d1)*0.5-w1*p2k1)*(p2k1*(t2-d2)-w2*r3)
     a -delta*(2.*p12*p2k1-w2*(t1-d1)))/p2k1
       return
20     dj=0.
       return
       end

      subroutine maps2(s2,x,smin,smax,ds)
       implicit double precision (a-z)
      y=smax/smin
      s2=smin*y**x
      ds=s2*dlog(y)
      return
      end

       subroutine mapla(x,y,z,u,xm,xp,d)
       implicit double precision (a-z)
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
       return
       end


      subroutine inplot(now,ff,pdx)
       implicit double precision ( a-h,o-z)
      common/lplot/xl(10),v1(2),v2(2),av(10)
      dimension zav(10),yav(10),zsv(10),ysv(10),ztv(10)
      dimension xlmax(10),xlmin(10),nlp(10),ltop(10),text(8,10),ll(10)
      dimension numb(12)
      dimension xls(42,10),yls(42,10),nlsn(42,10),mlsn(42,10),dls(10)
     1,xlav(10),xlsq(10),xlava(10),sxa(10),tlim(6),top(10),xltq(10)
      dimension nbin(41),nlog(41),slog(41),tlog(41),hv(12)
      dimension v1max(2),v1min(2),v2max(2),v2min(2),nv1(2)
     1,nv2(2),vtext(6,2)
      dimension vm(12,12,2),nvm(12,12,2),bin1(2),bin2(2),vol(2)
     1,wm(12,12,2),mvm(12,12,2)
      common/result/y,si,u,v
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
      if(igraph.le.0) go to 800
      if(igraph.ne.ngraph)read (12,810)nls
      if(igraph.ne.ngraph)print 814,nls
      if(nls.lt.0) nls=0
      if(nls.eq.0) go to 802
      if(nls.gt.mls) go to 807
      if(igraph.ne.ngraph)print 815
      do 801 i=1,nls
      if(igraph.ne.ngraph)
     1read (12,811)xlmin(i),xlmax(i),nlp(i),ltop(i),ll(i),
     2(text(j,i),j=1,8)
      if(igraph.ne.ngraph) print 816,i,xlmin(i),xlmax(i)
     1,nlp(i),ltop(i),ll(i),(text(j,i),j=1,8)
      if(nlp(i).lt.1)nlp(i)=1
      if(nlp(i).gt.40)nlp(i)=40
      dls(i)=(xlmax(i)-xlmin(i))/nlp(i)
      nlps=nlp(i)+2
      do 300 j=1,nlps
      yls(j,i)=0
300   mlsn(j,i)=0
801   continue
802   if(igraph.ne.ngraph) read (12,810)ndd
      if(igraph.ne.ngraph) print 817,ndd
      if(ndd.lt.0) ndd=0
      if(ndd.eq.0) go to 804
      if(ndd.gt.ndmax) go to 807
      if(igraph.ne.ngraph) print 818
      do 803 i=1,ndd
      if(igraph.ne.ngraph)
     1read (12,812)v1min(i),v1max(i),nv1(i),v2min(i),v2max(i),nv2(i)
     1,(vtext(j,i),j=1,6)
      if(igraph.ne.ngraph) print 819,i,v1min(i),v1max(i),nv1(i)
     1,v2min(i),v2max(i),nv2(i),(vtext(j,i),j=1,6)
      if(nv1(i).lt.1)nv1(i)=1
      if(nv2(i).lt.1)nv2(i)=1
      if(nv1(i).gt.10)nv1(i)=10
      if(nv2(i).gt.10)nv2(i)=10
      bin1(i)=(v1max(i)-v1min(i))/nv1(i)
      bin2(i)=(v2max(i)-v2min(i))/nv2(i)
      vol(i)=bin1(i)*bin2(i)
803   continue
      wtow=0.
      do 805 i=1,ndd
      do 805 j=1,12
      do 805 k=1,12
      wm(k,j,i)=0.
805   mvm(k,j,i)=0
804   continue
      if(igraph.ne.ngraph)read (12,810)nave
      if(igraph.ne.ngraph)print 820,nave
      if(nave.lt.0)nave=0
      if(nave.gt.mav)go to 807
      do 11 i=1,mav
      yav(i)=0.
11    ysv(i)=0.
      kt=0
      go to 808
800   nave=0
      nls=0
      ndd=0
      go to 808
807   print 813,mav,mls,ndmax
      stop
808   ngraph=igraph
      return
      entry replot(now,ff,pdx)
      if(nave.eq.0) go to 49
      do 62 i=1,nave
      zav(i)=0.
      ztv(i)=0.
62    zsv(i)=0.
49    fsqa=0.
      kt=kt+1
      if(nls.eq.0) go to 303
      do 302 i=1,nls
      nlps=nlp(i)+2
      xlav(i)=0
      xltq(i)=0.
      xlsq(i)=0
      do 302 j=1,nlps
      xls(j,i)=0
302   nlsn(j,i)=0
303   continue
      if(ndd.eq.0) go to 403
      do 402 i=1,ndd
      n1=nv1(i)+2
      n2=nv2(i)+2
      do 402 i1=1,n1
      do 402 i2=1,n2
      vm(i1,i2,i)=0
402   nvm(i1,i2,i)=0
403   continue
      return
      entry xplot(now,ff,pdx)
      fsqa=fsqa+ff*ff/pdx
      itt=itt+1
      if(nls.eq.0) go to 305
      do 304  i=1,nls
      nlps=(xl(i)-xlmin(i))/dls(i)+1.
      if(nlps.lt.0)nlps=0
      if(nlps.gt.nlp(i))nlps=nlp(i)+1
      nlps=nlps+1
      xls(nlps,i)=xls(nlps,i)+ff/dls(i)
      nlsn(nlps,i)=nlsn(nlps,i)+1
      xlav(i)=xlav(i)+ff*xl(i)
      xltq(i)=xltq(i)+ff*ff*xl(i)/pdx
304   xlsq(i)=xlsq(i)+(ff*xl(i))**2/pdx
305   continue
      if(ndd.eq.0)go to 405
      do 404 i=1,ndd
      i1=(v1(i)-v1min(i))/bin1(i)+2
      if(i1.lt.1) i1=1
      if(i1.gt.nv1(i)+2) i1=nv1(i)+2
      i2=(v2(i)-v2min(i))/bin2(i)+2
      if(i2.lt.1) i2=1
      if(i2.gt.nv2(i)+2) i2=nv2(i)+2
      vm(i1,i2,i)=vm(i1,i2,i)+ff/vol(i)
404   nvm(i1,i2,i)=nvm(i1,i2,i)+1
405   continue
      if(nave.eq.0)go to 99
      do 22 i=1,nave
      zav(i)=zav(i)+av(i)*ff
      ztv(i)=ztv(i)+ff*ff*av(i)/pdx
22    zsv(i)=zsv(i)+(av(i)*ff)**2/pdx
99    return
      entry plotit(now,ff,pdx)
      if(nls.eq.0)go to 315
      if(kk.gt.0)go to 307
      do 306 i=1,nls
      nlps=nlp(i)+2
      do 306 j=1,nlps
      mlsn(j,i)=nlsn(j,i)
306   yls(j,i)=xls(j,i)
      go to 310
307   vbef=vtot
      vu=(v/u)**2
      do 309 i=1,nls
      nlps=nlp(i)+2
      do 309 j=1,nlps
      if(nlsn(j,i).eq.0)go to 309
      if(mlsn(j,i).eq.0)go to 308
      al1=vu/nlsn(j,i)
      al2=vbef/mlsn(j,i)
      mlsn(j,i)=mlsn(j,i)+nlsn(j,i)
      yls(j,i)=(al2*xls(j,i)+al1*yls(j,i))/(al1+al2)
      go to 309
308   mlsn(j,i)=nlsn(j,i)
      yls(j,i)=xls(j,i)
309   continue
310   continue
      do 311 i=1,nls
      sxf=xlsq(i)-xlav(i)*xlav(i)
      sxt=xltq(i)-xlav(i)*u
      sx2=xlsq(i)/xlav(i)**2+fsqa/u**2-2.*xltq(i)/(xlav(i)*u)
      sx2=sx2*(xlav(i)/u)**2
      if(kt.ne.1)go to 312
      xlava(i)=xlav(i)/u
      sxa(i)=sx2
      go to 311
312   xhelp=sx2+sxa(i)
      if(xhelp.eq.0)go to 311
      xlava(i)=(xlav(i)*sxa(i)/u+xlava(i)*sx2)/xhelp
      sxa(i)=sxa(i)*sx2/xhelp
311   continue
      vtot=(si/y)**2
      if(now.ne.2)go to 315
      do 341 i=1,nls
      top(i)=0.
      nlps=nlp(i)+1
      do 341 j=2,nlps
      xls(j,i)=yls(j,i)/y
      if(xls(j,i).gt.top(i))top(i)=xls(j,i)
341   continue
      do 342 i=1,nls
      if(ltop(i).le.0)ltop(i)=i
      lto=ltop(i)
      if(top(i).gt.top(lto))top(lto)=top(i)
342   continue
      ylog=0.5*dlog10(y*y)
      do 314 i=1,nls
      print 321,i
      nlps=nlp(i)+1
      lto=ltop(i)
      top(i)=top(lto)
      if(top(i).eq.0)top(i)=1.
      an1=dlog10(top(i))
      n1=an1
      if(n1.gt.an1)n1=n1-1
      z1=top(i)*10.**(-n1)
      do 343 l=1,4
      if(z1.lt.tlim(l))go to 344
343   continue
      l=5
344   if(top(i).lt.1.6/(xlmax(i)-xlmin(i)))l=l+1
      topm=tlim(l)*10.**n1
      do 345 j=2,nlps
      nbin(j)=xls(j,i)*40./topm+1.5
      if(ll(i).lt.0)nbin(j)=0
      if(xls(j,i).le.0) go to 346
      tlog(j)=dlog10(xls(j,i))
      slog(j)=tlog(j)+ylog
      nlog(j)=(tlog(j)-n1)*8.+33.5
      if(ll(i).gt.0)nlog(j)=0
      go to 345
346   slog(j)=0
      tlog(j)=0
      nlog(j)=0
345   continue
      print 322,(text(j,i),j=1,8)
      n1p1=n1+1
      n1m4=n1-4
      print 323,tlim(l),n1,n1p1,n1m4
      do 347 l=1,40
      char(l)=hmin
      if(nlog(l+1).eq.41)char(l)=hplus
      if(nbin(l+1).eq.41)char(l)=hstar
347   continue
      xmin=xlmin(i)
      xmax=xmin+dls(i)
      print 324,xmin,xmax,yls(2,i),slog(2),xls(2,i),tlog(2),mlsn(2,i)
     1,char
***********************
      call pawfil2(   i,real(xmin),real(xmax),real(xls(2,i)))
      call pawfil2(10+i,real(xmin),real(xmax),real(yls(2,i)))
***********************
      do 348 j=3,nlps
      xmin=xmax
      xmax=xmin+dls(i)
      do 349 l=1,40
      char(l)=hblank
      if(nlog(l+1).eq.43-j)char(l)=hplus
      if(nbin(l+1).eq.43-j)char(l)=hstar
349   continue
      print 324,xmin,xmax,yls(j,i),slog(j),xls(j,i),tlog(j),mlsn(j,i)
     1,char
***********************
      call pawfil2(   i,real(xmin),real(xmax),real(xls(j,i)))
      call pawfil2(10+i,real(xmin),real(xmax),real(yls(j,i)))
***********************
348   continue
      nlps1=nlps+1
      if(nlps.eq.41)go to 352
      do 351 j=nlps1,41
      do 350 l=1,40
      char(l)=hblank
      if(nlog(l+1).eq.43-j)char(l)=hplus
      if(nbin(l+1).eq.43-j)char(l)=hstar
350   continue
351   print 325,char
352   do 353 l=1,40
      char(l)=hmin
      if(nlog(l+1).eq.1)char(l)=hplus
      if(nbin(l+1).eq.1)char(l)=hstar
353   continue
      print 326,char
      el1=yls(1,i)*dls(i)
      el2=el1/y
      print 327,el1,el2,mlsn(1,i)
      el1=yls(42,i)*dls(i)
      el2=el1/y
      print 328,el1,el2,mlsn(nlps1,i)
      sxsq=dsqrt(sxa(i)/itt)
      print 329,xlava(i),sxsq
314   continue
315   continue
      if(ndd.eq.0)go to 409
      wbef=wtot
      do 500 i=1,ndd
      nx=nv1(i)+2
      ny=nv2(i)+2
      if(kk.gt.0)go to 502
      do 501 j=1,nx
      do 501 k=1,ny
      wm(j,k,i)=vm(j,k,i)
501   mvm(j,k,i)=nvm(j,k,i)
      go to 500
502   vu=(v/u)**2
      do 503 j=1,nx
      do 503 k=1,ny
      if(nvm(j,k,i).eq.0)go to 503
      if(mvm(j,k,i).eq.0)go to 504
      al1=vu/nvm(j,k,i)
      al2=vbef/mvm(j,k,i)
      mvm(j,k,i)=mvm(j,k,i)+nvm(j,k,i)
      wm(j,k,i)=(al2*vm(j,k,i)+al1*wm(j,k,i))/(al1+al2)
      go to 503
504   mvm(j,k,i)=nvm(j,k,i)
      wm(j,k,i)=vm(j,k,i)
503   continue
500   continue
      wtot=(si/y)**2
      if(now.ne.2) go to 409
      do 408 i=1,ndd
      print 481,i,(vtext(j,i),j=1,6)
      vvv=v2max(i)
      mvv=nv1(i)+2
      nvv=nv2(i)+1
      size=vol(i)/y
      do 406 i2=1,nvv
      j2=nvv+2-i2
      do 410 i1=1,mvv
410   numb(i1)=1000.*wm(i1,j2,i)*size+.5
      print 486,(numb(i1),i1=1,mvv)
      print 483,(wm(i1,j2,i),i1=1,mvv)
      print 486,(mvm(i1,j2,i),i1=1,mvv)
      print 484,vvv
      vvv=vvv-bin2(i)
      if(dabs(vvv/bin2(i)).lt.1.e-10)vvv=0.
406   continue
      do 411 i1=1,mvv
411   numb(i1)=1000.*wm(i1,1,i)*size+.5
      print 486,(numb(i1),i1=1,mvv)
      print 483,(wm(i1,1,i),i1=1,mvv)
      print 486,(mvm(i1,1,i),i1=1,mvv)
      print 482
      mvv=mvv-1
      do 407 i1=1,mvv
      hv(i1)=v1min(i)+(i1-1)*bin1(i)
      if(dabs(hv(i1)/bin1(i)).lt.1.d-10)hv(i1)=0.
407   continue
      print 485,(hv(i1),i1=1,mvv)
408   continue
409   continue
      if(nave.eq.0) go to 23
      if(now.eq.2) print 26
      do 24 i=1,nave
      sxf=zsv(i)-zav(i)*zav(i)
      sxt=zsv(i)/zav(i)**2+fsqa/u**2-2.*ztv(i)/(zav(i)*u)
      sx2=sxt*(zav(i)/u)**2
      if(kt.ne.1) go to 21
      yav(i)=zav(i)/u
      ysv(i)=sx2
      go to 30
21    xhelp=sx2+ysv(i)
      if(xhelp.eq.0)go to 30
      yav(i)=(ysv(i)*zav(i)/u+yav(i)*sx2)/xhelp
      ysv(i)=ysv(i)*sx2/xhelp
30    yssq=dsqrt(ysv(i)/itt)
      if(now.eq.2)print 27,i,yav(i),yssq
24    continue
23    now=1
      kk=kk+1
      return
27    format(12x,i2,9x,d15.5,5x,d15.3)
26    format(1h1,10x,46hthe following are averages with error estimate/)
321   format(1h1,40x,40hsingle differential cross-section number,i3///)
322   format(38h single differential cross section of ,8a4/)
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
481   format(1h1,45x,40hdouble differential cross-section number,i3/
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
     2mensional plots//19x,i2,22h two dimensional plots////22x,25h***exe
     3cution is halted***   )
814   format(37h1number of single differential cross
     1,20hsections requested =,i3/)
815   format(30h information on the data cards//
     13h  i,10x,5hxlmin,12x,5hxlmax,7x,
     224hbins  correllation  type,19x,4htext/)
816   format(i3,2e17.4,i8,i9,i10,5x,8a4)
817   format(37h0number of double differential cross
     1,20hsections requested =,i3/)
818   format(30h information on the data cards//
     13h  i,8x,5hv1min,10x,5hv1max,4x,4hbins,8x,5hv2min
     2,10x,5hv2max,4x,4hbins,8x,6htext 1,8x,6htext 2/)
819   format(i3,2d15.3,i5,1x,2d15.3,i5,6x,3a4,2x,3a4)
820   format(31h0number of averages requested =,i3)
      end

       double precision function peripp(nup,ndown)
       implicit double precision (a-h,o-z)
       common/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
       common/levi/gram,d1,d2,d3,d4,d5,delta,g4,a1,a2
       common/civita/epsi,g5,g6,a5,a6,bb
       common/extra/s1,s2,t1,t2/dotps/q1dq,q1dq2,w6
       data rho/.585d+00/
       if(nup.gt.0)go to 1
       u1=1.
       u2=1.
       go to 3
1      if(nup.gt.1)go to 2
       xt=1.-t1/.71
       xt=xt*xt
       xu=2.79/xt
       u1=xu*xu
       tau=t1/(4.*w1)
       u2=(1./(xt*xt)-xu*xu*tau)/(1.-tau)
       go to 3
2      x=t1/(t1-w3)
       en=w31-t1
       tau=t1/(4.*w1)
       rhot=rho-t1
       u1=(-.86926*rho*rho*w31/rhot/rhot-2.23422*w1*(1.-x)**4
     1 /(x*(x*.96-1.26)+1.))/t1
       u2=(-tau*u1-.12549*w31*t1*rho/rhot/rhot*w31*w31/en/en/w1)
     1 /(1.-en*en/(4.*w1*t1))
3      if(ndown.gt.0)go to 4
       v1=1.
       v2=1.
       go to 6
4      if(ndown.gt.1)go to 5
       xt=1.-t2/.71
       xt=xt*xt
       xu=2.79/xt
       v1=xu*xu
       tau=t2/(4.*w2)
       v2=(1./(xt*xt)-xu*xu*tau)/(1.-tau)
       go to 6
5      x=t2/(t2-w5)
       en=w52-t2
       tau=t2/(4.*w2)
       rhot=rho-t2
       v1=(-.86926*rho*rho*w52/rhot/rhot-2.23422*w2*(1.-x)**4
     1 /(x*(x*.96-1.26)+1.))/t2
       v2=(-tau*v1-.12549*w52*t2*rho/rhot/rhot*w52*w52/en/en/w2)
     1 /(1.-en*en/(4.*w2*t2))
6      continue
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
       return
       end

       subroutine gamgam(ebeam,v1,v2,v3,v5,v6,v7,vmin,vmax,dj,nopt,x,nm)
       implicit double precision (a-h,o-z)
       common/variab/e,e1,e2,e3,e4,e5,p,p3,p4,p5,ct3,st3,ct4,st4,ct5,
     a st5,cp3,sp3,cp5,sp5/variac/al3,al4,be4,be5,de3,de5,pp3,pp4,pp5
       common/variad/e6,e7,p6,p7,ct6,st6,ct7,st7,cp6,sp6,cp7,sp7,w
       common/pickzz/w1,w2,w3,w4,w5,w31,w52,w12,tau,sl1
       common/dotps/q1dq,q1dq2,w6/extra/s1,s2,t1,t2
       common/dotp/p12,p13,p14,p15,p23,p24,p25,p34,p35,p45,p1k2,p2k1
       common/civita/epsi,g5,g6,a5,a6,bb
       common/ext/ctg,stg,cpg,spg
       common/angu/ctcm6,stcm6
       common/qvec/qve(4)
       dimension x(7)
       data pi/3.14159265358979d+00/,const/2.1868465d+10/
       w6=v6*v6
       w7=v7*v7
       wmin=v6+v7
       if(wmin.lt.vmin)wmin=vmin
       wmin=wmin*wmin
       e=2.*ebeam
       s=e*e
       wmax=e-v3-v5
       if(wmax.gt.vmax)wmax=vmax
       wmax=wmax*wmax
       xw=x(5)
       call mapw2(w4,xw,wmin,wmax,dw)
       v4=dsqrt(w4)
       w=v4
       call orient(s,v1,v2,v3,v4,v5,dj,nopt,x)
       if(t1.gt.0.or.t2.gt.0)dj=0.
       if(dj.eq.0)return
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
     a *e4*st4/v4+e4*ct4/v4*(p3*al3+e3mp3-e1mp1)
       pgp=dsqrt(pgx*pgx+pgy*pgy)
       pgg=dsqrt(pgp*pgp+pgz*pgz)
       if(pgg.gt.pgp*0.9.and.pgg.gt.pg)pg=pgg
       stg=pgp/pg
       cpg=pgx/pgp
       spg=pgy/pgp
       ctg=dsqrt(1.-stg*stg)
       if(pgz.lt.0)ctg=-ctg
       xx6=x(6)
       if(nm.eq.0)go to 1
       amap=.5*(w4-t1-t2)
       bmap=.5*dsqrt(((w4-t1-t2)**2-4.*t1*t2)*(1.-4.*w6/w4))
       ymap=(amap+bmap)/(amap-bmap)
       beta=ymap**(2.*xx6-1.)
       xx6=(amap/bmap*(beta-1.)/(beta+1.)+1.)*0.5
       if(xx6.gt.1.)xx6=1.
       if(xx6.lt.0.)xx6=0.
       ctcm6=1.-2.*xx6
       ddd=(amap+bmap*ctcm6)*(amap-bmap*ctcm6)/amap/bmap*dlog(ymap)
       dj=dj*ddd*0.5
1      ctcm6=1.-2.*xx6
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
     a /(e3*p+p3*ct3*e1)
       b1=(qx*sp5-qy*cp5)*pp5
       b2=(qz*e2+q0*p)*pp5
       b3=(w52*e2*e2+2.*w2*de5*e2-w2*de5*de5+pp5*pp5*e2*e2)
     a /(e2*p5*ct5-e5*p)
       r12=c2*sp3+qy*c3
       r13=-c2*cp3-qx*c3
       r22=b2*sp5+qy*b3
       r23=-b2*cp5-qx*b3
       epsi=p12*c1*b1+r12*r22+r13*r23
       g5=w1*c1*c1+r12*r12+r13*r13
       g6=w2*b1*b1+r22*r22+r23*r23
       a5=-(qx*cp3+qy*sp3)*pp3*p1k2-(e1*q0-p*qz)*(cp3*cp5+sp3*sp5)
     a *pp3*pp5+(de5*qz+q0*(p+p5*ct5))*c3
       a6=-(qx*cp5+qy*sp5)*pp5*p2k1-(e2*q0+p*qz)*(cp3*cp5+sp3*sp5)
     a *pp3*pp5+(de3*qz-q0*(p-p3*ct3))*b3
       return
       end

      subroutine mapw2(w2,x,w2min,w2max,dw)
       implicit double precision (a-z)
*      wmin=1./w2min
*      wmax=1./w2max
*      dw=wmin-wmax
*      w2=1./(wmax+dw*x)
*      dw=dw*w2*w2
       y = w2max/w2min
       w2 = w2min*y**x
       dw = w2*dlog(y)
      return
      end

      subroutine mapt1(t,x,tmin,tmax,dt)
       implicit double precision (a-z)
      y=tmax/tmin
      t=tmin*y**x
      dt=-t*dlog(y)
*      dt = tmin-tmax
*      t = tmin - x*dt
      return
      end

      subroutine mapt2(t,x,tmin,tmax,dt)
       implicit double precision (a-h,o-z)
      y=tmax/tmin
      t=tmin*y**x
      dt=-t*dlog(y)
*      dt = tmin-tmax
*      t = tmin - x*dt
      return
      end
