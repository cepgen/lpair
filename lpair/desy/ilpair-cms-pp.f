**********************************************************
*                ILPAIR ver. 4.0 - KRAKOW                *
* Program for analyzing e+e- pair in pp->pp+e+e- process *
* Created 03 January 2005, last update 14.01.2005        *
* Author: Dariusz Bocian@cern.ch                         *
*
* JC: correct for the use of the GEANT random numbe generator
*
**********************************************************
*--   
      EXTERNAL LUDATA,PYDATA
*--   
      COMMON /lujets/ N,K(4000,5),P(4000,5),V(4000,5)
*--   COMMON with invariant e+e- mass calculated with double precision in gmufil.f routine
      DOUBLE PRECISION pairm,PM
      COMMON /gmulpm/ pairm,PM(5,9)
*     -
*-----------------------------------------------------------------
*--   0 section (D.B.)
*     - my general
      INTEGER NEV,Nprt
      PARAMETER (NEV=1E5)             ! number of events NEV
      PARAMETER (Nprt=NEV/100)        ! printing period Nprt
*     - general constants
      DOUBLE PRECISION pi,SQ2,clight,qe,me,mekg
      PARAMETER (SQ2=1.4142135624)
      PARAMETER (pi=3.14159265)
      PARAMETER (clight=2.99792458E8) ! c=299792458 m/s
      PARAMETER (qe=1.6E-19)          ! q=1.6E-19 C
      PARAMETER (me=5.10998903E-4)    ! me=5.10998903E-4 GeV
      PARAMETER (mekg=9.10938188-31)  ! me=9.10938188-31 kg
*     - my constants
      DOUBLE PRECISION B_FIELD,l_solenoid
      PARAMETER (B_FIELD=2.0)         ! B=2.0 T
      PARAMETER (l_solenoid=7.0)      ! l=5.3 m in ATLAS
*     - my for calculation
*     
      INTEGER ITOTHIT,ITOTHIL,ITOTHITR
      INTEGER IHITL,IHITR,goodEVT
      INTEGER IHITNR,IHITNL,IHITCR,IHITCL
      DOUBLE PRECISION ETall
      DOUBLE PRECISION RIHITR,RIHITL
*     
      INTEGER itmax,FLAG,flagv
      DOUBLE PRECISION r_B,l_B,theta_B,tL,tR,dtR,dtmax
      DOUBLE PRECISION Rp,Rpn,phi,dphi,phiout,charge
      DOUBLE PRECISION etap,Pp,PTp,thp
      REAL phiv
*     
*     - ntuple creation
*     
      INTEGER ip,icode
      REAL px,py,pz,E,m,Ptot,PT,Eta
      REAL fin,fhit,Rhit,dthit,thB,iz,ihit
      INTEGER Lpawc,maxp,NNTUP
      PARAMETER (NNTUP=4444)
      PARAMETER (maxp=100)	! maxp=max number of particles in event
      PARAMETER (Lpawc=1500000)	! HBOOK PARAMETER
*     
*     - ntuple declaration
*     
      COMMON /QUEST/ IQUEST(100)
      COMMON /KINE/ ip,
     &     icode(maxp),                ! icode - particle code
     &     px(maxp),py(maxp),pz(maxp), ! px-z - momentum
     &     E(maxp),m(maxp),            ! E=Energy, m=mass
     &     Ptot(maxp),                 ! Ptot - momentum of particle	
     &     PT(maxp),                   ! PT - transversal momentum of particle
     &     fin(maxp),                  ! fin - phi from generator
     &     fhit(maxp),                 ! fhit - phi with B-field effect
     &     Rhit(maxp),                 ! Rhit - particle impact point with B-field effect
     &     dthit(maxp),thB(maxp),      ! dthit - difference of time of flight for "B-field" and "straight" particle
     &     iz(maxp),                   ! iz - the side indicator Left or Right
     &     ihit(maxp), Eta(maxp)       ! ihit - indicate the "slowest" particle i.e. with max. dthit

      INTEGER ie,Ntot,Nch,Nn,NhitL,NhitP
      INTEGER nh,nhl,nhr
      REAL ET
      COMMON /EVENT/ ie,
     &     Ntot,                ! Ntot - total number particles in event
     &     Nch,                 ! Nch - number of charged particles in event
     &     Nn,                  ! Nn - number of neutral particles in event
     &     NhitP,               ! NhitP - number of particles in "signal" side (0<eta<3)   
     &     NhitL,               ! NhitL - number of particles in "reference" side (-3<eta<0)
     &     nh,                  ! nh - total number of accepted hits in accepted events
     &     nhl,                 ! nhl- number of accepted hits in accepted events - left side
     &     nhr,                 ! nhr- number of accepted hits in accepted events - right side
     &     ET                   ! ET - transversal energy
                                     
      INTEGER iv
      REAL vx,vy,vz,dthetax,dthetay
      COMMON /PVERTEX/ iv,
     &     vx,vy,vz,            ! vx-z - production vertex coordinates
     &     thetax,dthetay       ! dthetax, dthetay - vertex smearing due to the beam divergence

*     -
      COMMON /PAWC/ hbook(Lpawc)
*     -
*     random number seeds initialisation
*     
      CALL get_seeds(1)
*     
      IQUEST(10)=128000
      
      CALL HLIMIT(Lpawc)
      CALL HCDIR('//PAWC',' ')
      CALL HROPEN(99,'CLBR','lpair.hbook','N',4096,IER)
      CALL HBNT(NNTUP,'CLBR',' ')
      CALL HBNAME(NNTUP,'KINE',ip,'ip[0,100]:I,icode(ip):I
     &     ,px(ip):R,py(ip):R,pz(ip):R,E(ip):R,m(ip):R
     &     ,Ptot(ip):R,PT(ip):R,fin(ip):R,fhit(ip):R,Rhit(ip):R
     &     ,dthit(ip):R,thB(ip):R,iz(ip):R,ihit(ip):R,Eta(ip):R')
      CALL HBNAME(NNTUP,'EVENT',ie,'ie:I,Ntot:I,Nch:I
     &     ,Nn:I,NhitP:I,NhitL:I,nh,nhl,nhr:I,ET:R')
      CALL HBNAME(NNTUP,'PVERTEX',iv,'iv:I
     &     ,vx:R,vy:R,vz:R,dthetax:R,dthetay:R')
*     
C...  First section: initialization, booking etc.
*     
*     
*     - control histograms evertex.f
*     
      CALL hbook1(101,'x?gen!',200,-5.E-5,5.E-5,0.)
      CALL hbook1(102,'y?gen!',200,-5.E-5,5.E-5,0.)
      CALL hbook1(103,'z?gen!',200,-0.2,0.2,0.)
*     
      CALL hbook1(105,'dtheta?x!',200,-5.E-5,5.E-5,0.)
      CALL hbook1(106,'dtheta?y!',200,-5.E-5,5.E-5,0.)
*     
*     - global histogram definition
*     
      CALL hbook1(1101,'N?p!/event L=0',5,1.,6.,0.) !Number of right events if L=0 
      CALL hbook1(1102,'N?p!/event L=1',5,1.,6.,0.) !Number of right events if L=1
      CALL hbook1(1103,'N?p!/event L=2',5,1.,6.,0.) !Number of right events if L=2
      CALL hbook1(1201,'N?p!/event R=2',5,-1.,4.,0.) !Number of left events if R=2
      CALL hbook1(1202,'N?p!/event R=3',5,-1.,4.,0.) !Number of left events if R=3
      CALL hbook1(1203,'N?p!/event R=4',5,-1.,4.,0.) !Number of left events if R=4
*     
c     CALL HBOOK2(2001,'N?hitR! vs N?hitL!',300,0.,3.,300,0.,3.,0.)
c     CALL HBOOK2(2002,'N?hitR! vs N?hitL!',300,0.,3.,300,0.,3.,0.)
c     CALL HBOOK2(2003,'N?hitR! vs N?hitL!',2,2.,4.,2,0.,2.,0.)
*     
c     CALL HBOOK1(3001,'N?tmax!',12,-1.,5.,0.)    !Number of hits with max. time
c     CALL HBOOK1(3002,'N?tmax!',12,-1.,5.,0.)    !Number of hits with max. time and 0.5 ns delay
c     CALL HBOOK1(3003,'N?tmax!',12,-1.,5.,0.)    !Number of hits with min. time
c     CALL HBOOK1(3004,'N?tmax!',12,-1.,5.,0.)    !Number of hits with min. time and 0.5 ns delay
c     CALL HBOOK2(3005,'N?tmax1!vsN?tmax1!'
c     & ,12,-1.,5.,12,-1.,5.,0.)    !Number of hits with both 0.5 ns delay
*     
c     CALL HBOOK2(4001,'t?1!-t?2! vs t?1!+t?2!'
c     &                      ,20,-1.,3.,30,-1.,5.,0.)
*     
*     
C...  Second section: event loop.
      ie=0
      iv=0
      ETall=0.
*     
*     
      write(*,*) ' Generator init'
*     
      CALL zduini

      DO IEV=1,NEV
         CALL zduevt(ireturn)
         IF(MOD(IEV,Nprt).EQ.0) print *,' Event nr = ',IEV
*     
C...  Extract and fill event properties.
         goodEVT = 0
         Npe=N
         IPN=0
         ip=0
         ITOTHIT=0
         ITOTHITL=0
         ITOTHITR=0
         ITOTHITN=0
         ITOTHITC=0
*     
         IRTHIT=0
         IRTHITL=0
         IRTHITR=0
*     - OK wazne
         IHITNR=0
         IHITNL=0
         IHITCR=0
         IHITCL=0
         IHITR=0
         IHITL=0
*     
         itmax=0
         FLAG=0
         flagv=0
         dtmax=0.
         ETall=0.
*     -
         DO I=1,N
            IF(K(I,1).EQ.1) THEN ! all stable particles
               IPN=IPN+1
*     
               etap=sign(log((sqrt(P(I,1)**2+P(I,2)**2+P(I,3)**2) !pseudorapidity calculation
     &              +abs(P(I,3)))/SQRT(P(I,1)**2+P(I,2)**2)),P(I,3))
               
*     
*     
               if (flagv.eq.0) then ! vertex calculation
                  phiv=atan2(P(I,2),P(I,1))
                  CALL EVERTEX(phiv)  
                  flagv=1
                  
               endif
*     
*     
               thp=ATAN(SQRT(P(I,1)**2+P(I,2)**2)/ABS(P(I,3))) !particle production angle
               phi=atan2(P(I,2),P(I,1)) !particle azimuthal angle
               if(phi.lt.0) phi = phi+2*pi
               Rp=0.5*l_solenoid*(SQRT(P(I,1)**2+P(I,2)**2)/ABS(P(I,3))) !particle hypothetic impact point
               Pp=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2) !particle momentum
               PTp=SQRT(P(I,1)**2+P(I,2)**2) !particle transverse momentum
               tR=((0.5*l_solenoid)/clight)*ABS((P(I,4)/P(I,3))) !time of the hit
*     -counters for ntuple
               ITOTHIT=ITOTHIT+1 !counts all particles in -3<eta<3
               
               
*------------------------------------
*     -ntuple filling - particles in event
*     
*     print *,' First pz = ',P(I,3)
*     IF((IHITNL+IHITCL).gt.2.and.(IHITNR+IHITCR).gt.4)GOTO 99
*     IF(K(I,2).eq.22.or.abs(charge).eq.1)THEN
*     IF(IHITNR.ge.1.or.IHITNL.ge.1.or.IHITCR.ge.1.or.IHITCL.ge.1) THEN
*     IF(etap.ge.-3.0.and.etap.le.3.0)THEN
*     IF(Rpn.ge.0.4.and.Rpn.le.1.4.and.theta_B.le.2.*pi) THEN
*     
*     -counters for ntuple
*     print *,' Second pz = ',P(I,3)
               IRTHIT=IRTHIT+1  !counts all particles in 0.4<R<1.4
               IF(etap.lt.0) IRTHITL=IRTHITL+1 !counts all particles in 0.4<R<1.4 and left side
               IF(etap.ge.0) IRTHITR=IRTHITR+1 !counts all particles in 0.4<R<1.4 and right side
*     -
*     
*     print *,' Third pz = ',P(I,3)
               ip=ip+1
               icode(ip)= K(I,2) ! icode - particle code
               px(ip)   = P(I,1) ! px - x-component of momentu
               py(ip)   = P(I,2) ! py - y-component of momentum
               pz(ip)   = P(I,3) ! pz - z-component of momentum
               E(ip)    = P(I,4) ! E - Energy of particle
               m(ip)    = P(I,5) ! m - mass of particle
               Ptot(ip) = Pp    ! Ptot - momentum of particle
               PT(ip)   = PTp   ! PT - transversal momentum of particle
               fin(ip)  = phi   ! fin - phi from generator
               fhit(ip) = mod(phiout,2.*pi) ! fhit - phi with B-field effect
               Rhit(ip) = Rpn   ! Rhit - partilce impact point with B-field effect
               dthit(ip)= dtR*1.E9 ! dthit - difference of time of flight for "B-field" and "stright" particle
               thB(ip)  = theta_B ! angle of curvarure in B-field
               iz(ip)   = etap/ABS(etap) ! iz - the side indicator Left or Right
               Eta(ip)  = etap  ! Eta - pseudorapidity 
*     print *,' Fourth pz = ',pz(ip)
               IF(ABS(charge).eq.0) ihit(ip) = 0 ! ihit - indicate the "slowest" particle i.e. with max. dthit
               IF(ABS(charge).eq.1) THEN
                  ihit(ip) = itmax  
                  IF(itmax.eq.1.and.dtmax.eq.dtR) ihit(ip) = 1
                  IF(itmax.eq.1.and.dtmax.gt.dtR) ihit(ip) = 0
                  IF(itmax.eq.1.and.dtmax.lt.dtR) THEN
                     dtmax=dtR
                     ihit(ip) = 1
                     IF(ip.gt.1) THEN
                        DO in=1,ip-1
                           ihit(in) = 0
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
*     print *,' Fourth pz = ',pz(ip)
*     
*     ENDIF ! IF(Rpn.ge.0.4.and.Rpn.le.1.4.and.theta_B.le.2*pi) THEN
*     ENDIF !IF(etap.ge.-3.0.and.etap.le.3.0)THEN - accept only particles in eta = (-3,3)
*     ENDIF !count eq. at least 1
*     ENDIF !IF(K(I,2).eq.22.or.abs(charge).eq.1)THEN
            ENDIF               ! IF(K(I,1).eq.1) - stable particles
 97         CONTINUE
         ENDDO
*     
*     - GENERAL HISTOGRAMS
*     
c     CALL HF2(2001,REAL(ITOTHITR),REAL(ITOTHITL),1.)
c     CALL HF2(2002,REAL(ITOTHITR),REAL(ITOTHITL),1.)
c     CALL HF2(2003,1.000000000000,0.000000000000,1.)
*     
*     - Second part of ntuple filling
*     
*     
         goodEVT = 1            ! filling ntuple
*     
*     - GENERAL HISTOGRAMS
*     
c     CALL HF2(4001,ABS(dthit(ip)-dthit(ip-1))
c     &                  ,dthit(ip)+dthit(ip-1),1.)
c     DO in=1,ip
c     IF(ihit(in).eq.1) CALL HF1(3001,dthit(in),1.)
c     IF(ihit(in).eq.1.and.dthit(in).ge.0.5)
c     &     CALL HF1(3002,dthit(in),1.)
c     IF(ihit(in).eq.0) CALL HF1(3003,dthit(in),1.)
c     IF(ihit(in).eq.0.and.dthit(in).ge.0.5)
c     &     CALL HF1(3004,dthit(in),1.)
c     ENDDO
c     IF(dthit(ip).ge.0.5.and.dthit(ip-1).ge.0.5)
c     &     CALL HF2(3005,dthit(ip),dthit(ip-1),1.)
c     
*     
*     -ntuple filling c.d. - event properties
*     
*     
*     -end of the loop
*-----------------------------------------------------------------
*     
         IF(goodEVT.eq.1) CALL HFNT(NNTUP)
C...  End event loop.
 99      CONTINUE
*     
*     store random number seeds every 500 events
*     
         if(mod(iev,500).eq.0) CALL put_seeds(0)
*     
      ENDDO
*-----------------------------------------------------------------
      CALL zduend
      CALL hrout(0,ICYCLE,' ')
      CALL hrend ('CLBR')
      CLOSE(99)
*     
*     store random number seeds
*     
      CALL put_seeds(1)
*     
      PRINT *,'                       '
      PRINT *,'      THE END          '
      PRINT *,'                       '
      END
      
      
      
      
