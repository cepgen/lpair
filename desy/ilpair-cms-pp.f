**********************************************************
*                                                        *
*       LPAIR ver. 4.1 - KRAKOW - LOUVAIN-LA-NEUVE       *
*                                                        *
* Monte Carlo generator to simulate                      *
*  (p/e)(p/e)->(p/p*/e)l+l-(p/e)                         *
* processes in High Energy Physics                       *
*                                                        *
* Created 03 January 2005, last update November 2013     *
* Authors/Collaborators :                                *
*   Jos Vermaseren <>                                    *
*   Dariusz Bocian <dariusz.bocian@cern.ch>              *
*   Janusz Chwastowski <>                                *
*   Laurent Forthomme <laurent.forthomme@cern.ch>        *
*                                                        *
**********************************************************
*   
      implicit none
*   
      integer ireturn,iev,i,in
*
      external ludata,pydata
      integer n
      double precision k,p,v
      common /lujets/ N,K(4000,5),P(4000,5),V(4000,5)
*
      double precision tmx
      common /mykin/ tmx
*
      double precision pairm,pm
      common /gmulpm/ pairm
*
*-----------------------------------------------------------------
*
      integer NEV,Nprt
      parameter (NEV=1E5)             ! number of events NEV
      parameter (Nprt=NEV/100)        ! printing period Nprt
*
      double precision pi
      parameter (pi=3.14159265)
*     
      double precision ETall
*     
      double precision phi,charge
      double precision etap,Pp,PTp,thp
*     
      integer ip,icode
      real px,py,pz,E,m,Ptot,PT,Eta,fin,iz
      integer Lpawc,maxp,NNTUP
      parameter (maxp=100)	! maxp=max number of particles in event
*     
      common /kine/ ip,
     &     icode(maxp),                ! icode - particle code
     &     px(maxp),py(maxp),pz(maxp), ! px-z  - 3-momentum
     &     E(maxp),                    ! E     - energy
     &     m(maxp),                    ! m     - mass
     &     Ptot(maxp),                 ! Ptot  - momentum of particle	
     &     PT(maxp),                   ! PT    - transversal momentum of particle
     &     fin(maxp),                  ! fin   - phi from generator
     &     iz(maxp),                   ! iz    - the side indicator Left or Right
     &     Eta(maxp)                   ! Eta   - pseudo-rapidity of the particle
*
      integer ie,Ntot,Nch,Nn
      real ET,MX
      data ntot,nn,nch/0,0,0/
      common /event/ ie,
     &     Ntot,                       ! Ntot  - total number particles in event
     &     Nch,                        ! Nch   - number of charged particles in event
     &     Nn,                         ! Nn    - number of neutral particles in event
     &     ET,                         ! ET    - transversal energy
     &     MX                          ! MX    - proton remnant mass

      real kchg,pmas,parf,vckm
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
*     
C...  Event loop.
      ie=0
      ETall=0.
*     
*     
      call zduini
*
      DO IEV=1,NEV
         call zduevt(ireturn)
c         IF(MOD(IEV,Nprt).EQ.0) print *,' Event nr = ',IEV
*     
         ip=0
*     
         ETall=0.
!     outgoing proton-like remnants invariant mass
         MX=TMX
         print *,pairm,tmx
*
         DO I=1,N
            IF(K(I,1).EQ.1) THEN ! all stable particles
*     
!     pseudorapidity calculation
               etap=sign(log((sqrt(P(I,1)**2+P(I,2)**2+P(I,3)**2)
     &              +abs(P(I,3)))/SQRT(P(I,1)**2+P(I,2)**2)),P(I,3))
               
!     particle production angle
               thp=ATAN(SQRT(P(I,1)**2+P(I,2)**2)/ABS(P(I,3)))
!     particle azimuthal angle
               phi=atan2(P(I,2),P(I,1))
               IF(phi.LT.0) phi = phi+2*pi
!     particle momentum
               Pp=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
!     particle transverse momentum
               PTp=SQRT(P(I,1)**2+P(I,2)**2)
*     
*------------------------------------
*     
               ip=ip+1
               icode(ip)= K(I,2)         ! icode - particle code
               px(ip)   = P(I,1)         ! px    - x-component of momentum
               py(ip)   = P(I,2)         ! py    - y-component of momentum
               pz(ip)   = P(I,3)         ! pz    - z-component of momentum
               E(ip)    = P(I,4)         ! E     - Energy of particle
               m(ip)    = P(I,5)         ! m     - mass of particle
               Ptot(ip) = Pp             ! Ptot  - momentum of particle
               PT(ip)   = PTp            ! PT    - transversal momentum of particle
               fin(ip)  = phi            ! fin   - phi from generator
               iz(ip)   = etap/abs(etap) ! iz    - the side indicator Left or Right
               Eta(ip)  = etap           ! Eta   - pseudorapidity 
               if(abs(charge).eq.0) then
                  nn=nn+1
               else
                  nch=nch+1
               endif
               write(*,1000) px(ip),py(ip),pz(ip)
*     
            endif
 97         continue
         enddo
*     
*-----------------------------------------------------------------
*     
 99      continue
*     
      enddo
*     
*-----------------------------------------------------------------
*     
      call zduend
 1000 format(d20.8,d20.8,d20.8)
*     
      end
