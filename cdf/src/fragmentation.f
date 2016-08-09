c     Fragmentation algorithm
c     Connects LPAIR output to Jetset for string fragmentation
c
c     Authors :
c     ...
c     N. Schul (UCL, Louvain-la-Neuve)
c     L. Forthomme (UCL, Louvain-la-Neuve), Feb 2014

      SUBROUTINE FRAGMENTATION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c     Pair produced code
      INTEGER IPAIR,NLHE
      PARAMETER (IPAIR=24) 
      INTEGER L,NLEP 

      INTEGER NPART,NLIMAX
      PARAMETER (NLIMAX=13)
      REAL*4  PL(4,NLIMAX),PMXDA(4),PMXDB(4)
      REAL*4  I2MASS(NLIMAX)
      DATA    I2MASS/NLIMAX*-9999.9/

c     LUND common
      REAL*4        P(4000,5)
      INTEGER       N,K(4000,5)
      COMMON/LUJETS/N,K,P
      save /lujets/

c     Parameters pi
      REAL*4  PI2,PI
      PARAMETER (PI2=2.0*3.14159265)
      PARAMETER (PI=3.14159265)
c     For LUJOIN      
      INTEGER JLPSF1(2),JLPSF2(2)
      DATA    JLPSF1/10,11/
      DATA    JLPSF2/12,13/

c     for MX
      REAL*8  MX1,MX2,WX1,WX2,W1,W6,W3,W8
      INTEGER I2STAT(13),I2PART(13)
      INTEGER I2MO1(13),I2DA1(13),I2DA2(13)
      REAL*4 RANUDQ3,RANUDQ4
      COMMON/REMNTS/MX1,MX2

c     Allow W+ -> l+ v only
      COMMON/LUDAT3/MDME(8750,2)
      SAVE  /LUDAT3/
      
c     HEPEVT common
      PARAMETER (NMXHEP=10000)
      COMMON/HEPEVT/NEVHEP,NHEP,
     &     ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &     JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),
     &     PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      REAL PHEP,VHEP

c     Event common
      LOGICAL ACCEPTED
      INTEGER NDIM
      DOUBLE PRECISION x(10)
      COMMON/event/accepted,ndim,x

      INTEGER I,J
      LOGICAL FRAGMENT(2)
!
! INCOMING p+_1    == 1
! photon 1         == 2
! DISSOCIATED p+_1 == 3
! Quark/Diquark    == 4&5

!
! INCOMING p+_2    == 6
! photon 2         == 7
! DISSOCIATED p+_2 == 8
! Quark/Diquark    == 9&10

***** PHEP common block content
*   1  proton       (incoming, left)
*   2  proton       (incoming,right)
*   3  proton       (outgoing, left)
*   4  muon-pair (resonance product)
*   5  proton       (outgoing,right)
*   6  muon         (outgoing, left)
*   7  muon         (outgoing,right)
*
*
*    1 ------x------------ 3
*        p    \    p/q
*          gam \
*               \_________ 6
*               |  ell+-
*               | 4
*               |_________ 7
*               /  ell-+
*          gam /
*             /
*    2 ------x------------ 5
*        p         p/q
*
*****

      IF(NDIM.EQ.7) THEN
         FRAGMENT(1)=.FALSE.
         FRAGMENT(2)=.FALSE.
      ELSEIF(NDIM.EQ.8) THEN
         FRAGMENT(1)=.TRUE.
         FRAGMENT(2)=.FALSE.
      ELSEIF(NDIM.EQ.9) THEN
         FRAGMENT(1)=.TRUE.
         FRAGMENT(2)=.TRUE.
      ELSE
         GOTO 22
      ENDIF

c      PRINT *,'To fragment :',(FRAGMENT(I),I=1,2)

c...  READ the kinematics from the event
      NFAIL=0
 1    J=rand(0)*6561+1.
      NFAIL=NFAIL+1
      IF(NFAIL.GT.9999) THEN
         WRITE(*,*) 'SKIP Event'
         RETURN
      ENDIF

***** PL (Lund) kinematic quantities mapping
*
*    1 ------x------------ 5
*        p    \    p/q
*            3 \ gam
*               \_________ 8
*               |  ell+-
*               | 6
*               |_________ 9
*               /  ell-+
*            4 / gam
*             /
*    2 ------x------------ 7
*        p         p/q
*
*****

      NPART=9
c...  Fill the Lund common block
C     PARTICLE 1 = "PROTON1" <===================================
      PL(1,1)=PHEP(1,1)
      PL(2,1)=PHEP(2,1)
      PL(3,1)=PHEP(3,1)
      PL(4,1)=PHEP(4,1)
C     PARTICLE 2 = "PROTON2" <===================================
      PL(1,2)=PHEP(1,2)
      PL(2,2)=PHEP(2,2)
      PL(3,2)=PHEP(3,2)
      PL(4,2)=PHEP(4,2)
C     PARTICLE 5 = QUARK1 OUT <==================================
      PL(1,5)=PHEP(1,3)
      PL(2,5)=PHEP(2,3)
      PL(3,5)=PHEP(3,3)
      PL(4,5)=PHEP(4,3)
C     PARTICLE 7 = QUARK2 OUT <==================================
      PL(1,7)=PHEP(1,5)
      PL(2,7)=PHEP(2,5)
      PL(3,7)=PHEP(3,5)
      PL(4,7)=PHEP(4,5)
C     PARTICLE 3 = GAMMA1
      PL(1,3)=PHEP(1,1)-PHEP(1,3)
      PL(2,3)=PHEP(2,1)-PHEP(2,3)
      PL(3,3)=PHEP(3,1)-PHEP(3,3)
      PL(4,3)=PHEP(4,1)-PHEP(4,3)
C     PARTICLE 4 = GAMMA2
      PL(1,4)=PHEP(1,2)-PHEP(1,5)
      PL(2,4)=PHEP(2,2)-PHEP(2,5)
      PL(3,4)=PHEP(3,2)-PHEP(3,5)
      PL(4,4)=PHEP(4,2)-PHEP(4,5)
C     Particle 6 --> middle particle
      PL(1,6)=PHEP(1,4)
      PL(2,6)=PHEP(2,4)
      PL(3,6)=PHEP(3,4)
      PL(4,6)=PHEP(4,4)
CLF   PARTICLE 8 = LEPTON1
      PL(1,8)=PHEP(1,6)
      PL(2,8)=PHEP(2,6)
      PL(3,8)=PHEP(3,6)
      PL(4,8)=PHEP(4,6)
CLF   PARTICLE 9 = LEPTON2
      PL(1,9)=PHEP(1,7)
      PL(2,9)=PHEP(2,7)
      PL(3,9)=PHEP(3,7)
      PL(4,9)=PHEP(4,7)

c      DO 2000 I=1,7
c         PRINT *,'PHEP:',I,':',(PHEP(J,I),J=1,4)
c 2000 CONTINUE
c      PRINT *,(PHEP(I,1)+PHEP(I,2),I=1,4)
c      PRINT *,(PHEP(I,3)+PHEP(I,5)+PHEP(I,6)+PHEP(I,7),I=1,4)
c      PRINT *,(PHEP(I,1)+PHEP(I,2)-
c     &     (PHEP(I,3)+PHEP(I,5)+PHEP(I,6)+PHEP(I,7)),I=1,4)
c      PRINT *,''
c      PRINT *,'PHEP3=',(PHEP(I,3),I=1,4)
c      PRINT *,'PHEP5=',(PHEP(I,5),I=1,4)
c      PRINT *,(PHEP(I,6),I=1,4)
c      PRINT *,(PHEP(I,7),I=1,4)
c      PRINT *,'PL5=',(PL(I,7),I=1,4)
c      PRINT *,'PL7=',(PL(I,7),I=1,4)
c      PRINT *,PHEP(5,3),PHEP(5,5)
      
c     Default status for the outgoing protons is 1 ('stable')
      I2STAT(5)=1
      I2STAT(7)=1

c...  choose the quark content

      IF(FRAGMENT(1)) THEN
c...  compute the MX value
         WX1=(1.1449)*(102400/1.1449)**X(8)
         MX1=DSQRT(WX1)
c     smearing of the mass (N. Schul)
c     WRITE(*,*) 'mx1, mx2=',MX1,MX2
         VARYMX1=rand(0)*(0.1*MX1)
         MX1=MX1+VARYMX1
c     WRITE(*,*) 'corr mx1=',MX1

C====> INSERT THE MASS OF THE HADRONIC SYSTEM <==================
         I2MASS(5)=SNGL(MX1)
C===> RANDOM SELECTION OF U , D AND DI QUARKS <==================
         RANUDQ4=rand(0)
         IF (RANUDQ4 .LT. 1.0/9.0) THEN
            I2PART(10)=1
            I2PART(11)=2203
         ELSEIF (RANUDQ4 .LT. 5.0/9.0) THEN
            I2PART(10)=2
            I2PART(11)=2101
         ELSE
            I2PART(10)=2
            I2PART(11)=2103
         ENDIF
         ULMQ4=ULMASS(I2PART(10))
         ULMDQ4=ULMASS(I2PART(11))
         I2STAT(5)=21
      ENDIF
      
      IF(FRAGMENT(2)) THEN
c...  compute the MX value
      WX2=(1.1449)*(102400/1.1449)**X(9)
      MX2=DSQRT(WX2)
c     smearing of the mass (N. Schul)
c     WRITE(*,*) 'mx1, mx2=',MX1,MX2
      VARYMX2=rand(0)*(0.1*MX2)
      MX2=MX2+VARYMX2
c     WRITE(*,*) 'corr mx2=',MX2
      
C====> INSERT THE MASS OF THE HADRONIC SYSTEM <==================
         I2MASS(7)=SNGL(MX2)
C===> RANDOM SELECTION OF U , D AND DI QUARKS <==================
         RANUDQ3=rand(0)
         IF (RANUDQ3 .LT. 1.0/9.0) THEN
            I2PART(12)=1
            I2PART(13)=2203
         ELSEIF (RANUDQ3 .LT. 5.0/9.0) THEN
            I2PART(12)=2
            I2PART(13)=2101
         ELSE
            I2PART(12)=2
            I2PART(13)=2103
         ENDIF
         ULMQ3=ULMASS(I2PART(12))
         ULMDQ3=ULMASS(I2PART(13))
         I2STAT(7)=21
      ENDIF

C...  Set Lund code
*
*  1 : incoming proton 1
      I2STAT(1)=21
      I2PART(1)=2212
      I2MO1(1)=0
      I2DA1(1)=3
      I2DA2(1)=5
         
*  2 : incoming proton 2
      I2STAT(2)=21
      I2PART(2)=2212
      I2MO1(2)=0
      I2DA1(2)=4
      I2DA2(2)=7
         
*  3 : photon 1
      I2STAT(3)=11
      I2PART(3)=22
      I2MO1(3)=1
      I2DA1(3)=0
      I2DA2(3)=0
         
*  4 : photon 2
      I2STAT(4)=11
      I2PART(4)=22
      I2MO1(4)=2
      I2DA1(4)=6
      I2DA2(4)=0
         
*  5 : outgoing proton 1
      I2PART(5)=2212
      I2MO1(5)=1
      I2DA1(5)=0
      I2DA2(5)=0
         
*  6 : central system (mother : photon 2)
*     I2STAT(6)=1
*     I2PART(6)=14 ! muon neutrino ?!
      I2STAT(6)=11
      I2PART(6)=100
      I2MO1(6)=4
      I2DA1(6)=8
      I2DA2(6)=9
         
*  7 : outgoing proton 2
      I2PART(7)=2212
      I2MO1(7)=2
      I2DA1(7)=0
      I2DA2(7)=0
         
*  8 : muon 1 !!! need to check the charge
      I2STAT(8)=1
      I2PART(8)=-13
      I2MO1(8)=6
      I2DA1(8)=0
      I2DA2(8)=0

*  9 : muon 2 !!! need to check the charge
      I2STAT(9)=1
      I2PART(9)=13
      I2MO1(9)=6
      I2DA1(9)=0
      I2DA2(9)=0

      IF(FRAGMENT(1)) THEN
* 10 : quark 1
         I2STAT(10)=1
         I2MO1(10)=5
         I2DA1(10)=0
         I2DA2(10)=0
         
* 11 : diquark 1
         I2STAT(11)=1
         I2MO1(11)=5
         I2DA1(11)=0
         I2DA2(11)=0

C==== > CHOOSE RANDOM DIRECTION IN MX FRAME <====================
         RANMXP1=PI2*rand(0)
         RANMXT1=ACOS(2.0*rand(0)-1.0)
         
C==== > COMPUTE MOMENTUM OF DECAY PARTICLE FROM MX2 <============
         PMXP1=DSQRT((MX1**2-ULMDQ4**2+ULMQ4**2)**2/
     &        4.0/MX1/MX1-ULMQ4**2)
         
         PMXDB(1)=SIN(RANMXT1)*COS(RANMXP1)*PMXP1
         PMXDB(2)=SIN(RANMXT1)*SIN(RANMXP1)*PMXP1
         PMXDB(3)=COS(RANMXT1)*PMXP1
         PMXDB(4)=SQRT(PMXP1**2+ULMDQ4**2)      
         CALL LORENB(I2MASS(5),PL(1,5),PMXDB(1),PL(1,11))
         
         PMXDB(1)=-PMXDB(1)
         PMXDB(2)=-PMXDB(2)
         PMXDB(3)=-PMXDB(3)
         PMXDB(4)=SQRT(PMXP1**2+ULMQ4**2)
         CALL LORENB(I2MASS(5),PL(1,5),PMXDB(1),PL(1,10))

c     We added 2 'particles' in the event (quark and diquark)
         NPART=NPART+2

      ENDIF
         
      IF(FRAGMENT(2)) THEN
* 12 : quark 2
         I2STAT(12)=1
         I2MO1(12)=7
         I2DA1(12)=0
         I2DA2(12)=0
         
* 13 : diquark 2
         I2STAT(13)=1
         I2MO1(13)=7
         I2DA1(13)=0
         I2DA2(13)=0

C==== > CHOOSE RANDOM DIRECTION IN MX FRAME <====================
         RANMXP2=PI2*rand(0)
         RANMXT2=ACOS(2.0*rand(0)-1.0)
         
C==== > COMPUTE MOMENTUM OF DECAY PARTICLE FROM MX1 <============
         PMXP2=DSQRT((MX2**2-ULMDQ3**2+ULMQ3**2)**2/
     &        4.0/MX2/MX2-ULMQ3**2)
         
         PMXDA(1)=SIN(RANMXT2)*COS(RANMXP2)*PMXP2
         PMXDA(2)=SIN(RANMXT2)*SIN(RANMXP2)*PMXP2
         PMXDA(3)=COS(RANMXT2)*PMXP2
         PMXDA(4)=SQRT(PMXP2**2+ULMDQ3**2)
         CALL LORENB(I2MASS(7),PL(1,7),PMXDA(1),PL(1,13))
         
         PMXDA(1)=-PMXDA(1)
         PMXDA(2)=-PMXDA(2)
         PMXDA(3)=-PMXDA(3)
         PMXDA(4)=SQRT(PMXP2**2+ULMQ3**2)
         CALL LORENB(I2MASS(7),PL(1,7),PMXDA(1),PL(1,12))

c     We added 2 'particles' in the event (quark and diquark)
         NPART=NPART+2
      
      ENDIF
*
         
         
c      DO 2001 I=1,NPART
c         PRINT *,'PL:',I,':',(PL(J,I),J=1,4)
c 2001 CONTINUE
      
      CALL LUNSET(NPART)
C     ====> FILLING THE LUND COMMON <============================
      DO 201 I=1,NPART
C     SET MOTHER/DAUGHTER VALUES, MARKING PARTICLES AS DECAYED <=
         CALL LUKSET(I,I2STAT(I),I2PART(I),
     &        I2MO1(I),I2DA1(I),I2DA2(I),0)
C     SET PULS, ENERGY AND MASS OFF THE PARTICLES <==============
         CALL LUPSET(I,PL(1,I),PL(2,I),PL(3,I),PL(4,I),I2MASS(I))
 201  CONTINUE
      
      IF(FRAGMENT(1)) CALL LUJOIN(2,JLPSF1)
      IF(FRAGMENT(2)) CALL LUJOIN(2,JLPSF2)
      
      CALL LUEXEC
c     IF(K(16,2).EQ.91.AND.K(16,1).EQ.11) THEN
c     WRITE(*,*) '-> System non-inelastic'
c     GO TO 1
c     ENDIF
clf      IF(K(17,2).EQ.2212.AND.K(17,1).EQ.1) THEN
clf         WRITE(*,*) '-> System non-inelastic'
clf         NFAIL=NFAIL+1
clf         GO TO 1
clf      ENDIF
      
c      CALL LULIST(2)

      CALL LHEFIL
      RETURN

 22   PRINT *,'ERROR! The number of dimensions is incorrect'
      END
      
