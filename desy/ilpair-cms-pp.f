**********************************************************
*                                                        *
*       LPAIR ver. 4.1 - KRAKOW--LOUVAIN-LA-NEUVE        *
*                                                        *
* Monte Carlo generator to simulate                      *
*  (p/e)(p/e)->(p/p*/e)l+l-(p/e)                         *
* processes                                              *
*                                                        *
* Created 03 January 2005, last update November 2013     *
* Authors/Collaborators :                                *
*   Jos Vermaseren <>                                    *
*   Dariusz Bocian <dariusz.bocian@cern.ch>              *
*   Janusz Chwastowski <>                                *
*   Laurent Forthonne <laurent.forthomme@cern.ch>        *
*                                                        *
**********************************************************
*--   
      implicit none
*--   
      integer ireturn,iev,i,in,ipn

      external ludata,pydata
      integer n
      double precision k,p,v
      common /lujets/ N,K(4000,5),P(4000,5),V(4000,5)

      double precision tmx
      common /mykin/ tmx

      double precision pairm,pm
      common /gmulpm/ pairm,pm(5,9)
*     -
*-----------------------------------------------------------------
*--   0 section (D.B.)
*     - my general
      INTEGER NEV,Nprt
      PARAMETER (NEV=1E5)             ! number of events NEV
      PARAMETER (Nprt=NEV/100)        ! printing period Nprt
*     - general constants
      DOUBLE PRECISION pi,me
      PARAMETER (pi=3.14159265)
      PARAMETER (me=5.10998903E-4)    ! me=5.10998903E-4 GeV
*     - my for calculation
*     
      INTEGER goodEVT
      DOUBLE PRECISION ETall
*     
      INTEGER itmax,FLAG
      DOUBLE PRECISION tL,tR,dtmax
      DOUBLE PRECISION phi,dphi,charge
      DOUBLE PRECISION etap,Pp,PTp,thp
*     
*     - ntuple creation
*     
      INTEGER ip,icode
      REAL px,py,pz,E,m,Ptot,PT,Eta
      REAL fin,fhit,Rhit,dthit,thB,iz
      INTEGER Lpawc,maxp,NNTUP
      PARAMETER (NNTUP=4444)
      PARAMETER (maxp=100)	! maxp=max number of particles in event
      PARAMETER (Lpawc=30000)	! HBOOK PARAMETER
*     
*     - ntuple declaration
*     
      integer iquest
      COMMON /QUEST/ IQUEST(100)
      COMMON /KINE/ ip,
     &     icode(maxp),                ! icode - particle code
     &     px(maxp),py(maxp),pz(maxp), ! px-z - momentum
     &     E(maxp),m(maxp),            ! E=Energy, m=mass
     &     Ptot(maxp),                 ! Ptot - momentum of particle	
     &     PT(maxp),                   ! PT - transversal momentum of particle
     &     fin(maxp),                  ! fin - phi from generator
     &     fhit(maxp),                 ! fhit - phi with B-field effect
     &     iz(maxp),                   ! iz - the side indicator Left or Right
     &     Eta(maxp)                   ! Eta - indicate the pseudo-rapidity of the particle

      INTEGER ie,Ntot,Nch,Nn
      REAL ET,MX
      data ntot,nn,nch/0,0,0/
      COMMON /EVENT/ ie,
     &     Ntot,                ! Ntot - total number particles in event
     &     Nch,                 ! Nch - number of charged particles in event
     &     Nn,                  ! Nn - number of neutral particles in event
     &     ET,                  ! ET - transversal energy
     &     MX                   ! MX - proton remnant mass
                                     
*     -
      double precision hbook
      integer icycle,ier
      COMMON /PAWC/ hbook(Lpawc)
*     -
      IQUEST(10)=128000
      
c      CALL HLIMIT(Lpawc)
c      CALL HCDIR('//PAWC',' ')
c      CALL HROPEN(99,'CLBR','lpair.hbook','N',4096,IER)
c      CALL HBNT(NNTUP,'CLBR',' ')
c      CALL HBNAME(NNTUP,'KINE',ip,'ip[0,100]:I,icode(ip):I
c     &     ,px(ip):R,py(ip):R,pz(ip):R,E(ip):R,m(ip):R,Ptot(ip):R,PT(ip):R
c     &     ,fin(ip):R,fhit(ip):R,iz(ip):R,Eta(ip):R')
c      CALL HBNAME(NNTUP,'EVENT',ie,'ie:I,Ntot:I,Nch:I,Nn:I,ET:R,MX:R')
*     
*     
C...  Event loop.
      ie=0
      ETall=0.
*     
*     
      write(*,*) ' Generator init'
*     
      CALL zduini

      DO IEV=1,NEV
c         print *,'Calling ZDUEVT for the',IEV,'th time'
         CALL zduevt(ireturn)
         IF(MOD(IEV,Nprt).EQ.0) print *,' Event nr = ',IEV
*     
C...  Extract and fill event properties.
         goodEVT = 0
         IPN=0
         ip=0
*     
*     
         itmax=0
         FLAG=0
         dtmax=0.
         ETall=0.
C LF: Outgoing proton-like remnants invariant mass
         MX=TMX
*     -
         DO I=1,N
            IF(K(I,1).EQ.1) THEN ! all stable particles
               IPN=IPN+1
*     
               etap=sign(log((sqrt(P(I,1)**2+P(I,2)**2+P(I,3)**2) !pseudorapidity calculation
     &              +abs(P(I,3)))/SQRT(P(I,1)**2+P(I,2)**2)),P(I,3))
               
*     
!     particle production angle
               thp=ATAN(SQRT(P(I,1)**2+P(I,2)**2)/ABS(P(I,3)))
!     particle azimuthal angle
               phi=atan2(P(I,2),P(I,1))
               IF(phi.LT.0) phi = phi+2*pi
!     particle momentum
               Pp=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
!     particle transverse momentum
               PTp=SQRT(P(I,1)**2+P(I,2)**2)
               
*------------------------------------
*     ntuple filling - particles in event
*     
*     print *,' First pz = ',P(I,3)
*     
*     counters for ntuple
*     print *,' Second pz = ',P(I,3)
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
               iz(ip)   = etap/ABS(etap) ! iz - the side indicator Left or Right
               Eta(ip)  = etap  ! Eta - pseudorapidity 
*     print *,' Fourth pz = ',pz(ip)
               if(abs(charge).eq.0) then
                  nn=nn+1
               else
                  nch=nch+1
               ENDIF
*     print *,' Fourth pz = ',pz(ip)
*     
            ENDIF
 97         CONTINUE
         ENDDO
*     
         goodEVT = 1            ! filling ntuple
*     
*     -end of the loop
*-----------------------------------------------------------------
*     
c         IF(goodEVT.eq.1) CALL HFNT(NNTUP)
C...  End event loop.
 99      CONTINUE
*     
      ENDDO
*-----------------------------------------------------------------
      CALL zduend
c      CALL hrout(0,ICYCLE,' ')
c      CALL hrend ('CLBR')
      CLOSE(99)
*     
      PRINT *,'                       '
      PRINT *,'      THE END          '
      PRINT *,'                       '
      END
      
      
      
      
