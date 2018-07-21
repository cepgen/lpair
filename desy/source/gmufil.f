      SUBROUTINE GMUFIL
C
C   THIS SUBROUTINE SHOULD FILL THE GTR-BANK
C
      IMPLICIT none

      DOUBLE PRECISION ulmass
      DOUBLE PRECISION ulmq,ulmdq
      DOUBLE PRECISION pmxp,ranmxp,ranmxt
      DOUBLE PRECISION gmuselx
      INTEGER ncall,nextw,nfracs,nfrac3,nterm3,ninit,nfinal,npout
      INTEGER ipq,ipdq,i2line
      INTEGER i,j

      REAL*4  PI2
      PARAMETER (PI2=2.0*3.14159265)
*KEEP,INPU.
      REAL*8       ME,MU,MP,MX,S,SQ,PE,PP,EE,EP,CONST,PI
      COMMON /INPU/ME,MU,MP,MX,S,SQ,PE,PP,EE,EP,CONST,PI
      save /inpu/

*KEND.
      REAL*8         E,E1,E2,E3,E4,E5,P1,
     &                P3,P4,P5,CT3,ST3,CT4,ST4,CT5,
     &                ST5,CP3,SP3,CP5,SP5
      COMMON /VARIAB/E,E1,E2,E3,E4,E5,P1,
     &                P3,P4,P5,CT3,ST3,CT4,ST4,CT5,
     &                ST5,CP3,SP3,CP5,SP5
      REAL*8          E6,E7,P6,P7,CT6,ST6,CT7,ST7,CP6,SP6,CP7,SP7,W
      COMMON /VARIAD/ E6,E7,P6,P7,CT6,ST6,CT7,ST7,CP6,SP6,CP7,SP7,W
      REAL*8          S1,S2,T1,T2
      COMMON /EXTRA/  S1,S2,T1,T2
*KEEP,LTCOM.
      REAL*8          GAMMA,BETGAM
      COMMON /LTCOM/  GAMMA,BETGAM
      save /variab/,/variad/,/extra/

*KEEP,XQCOM.
      INTEGER IUSEDF
      DOUBLE PRECISION XQ,EQ,PQ,MQ,
     &                 QSCALE,XDENS(-6:2),
     &                 PSEA,PVALD
      COMMON /XQCOM/   XQ,EQ,PQ,MQ,
     &                 QSCALE,XDENS,
     &                 PSEA,PVALD,IUSEDF

*KEEP,BEAM.
      INTEGER          INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,IPAIR,NQUARK
      REAL*8           INPE,INPP
      COMMON /BEAM/    INPE,INPP,INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,
     &                 IPAIR,NQUARK
      save /beam/

*KEEP,KINVAR.
      REAL*8           GMUX,GMUY,GMUW,GMUNU
      COMMON /KINVAR/  GMUX,GMUY,GMUW,GMUNU
      save /kinvar/
*KEND.
      INTEGER NLIMAX
      PARAMETER (NLIMAX=13)
      REAL*4 PL(4,NLIMAX),PMXDA(4)
C
C   LUND COMMON <===================================================
      double precision kchg,pmas,parf,vckm
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      DOUBLE PRECISION P(4000,5),V(4000,5)
      INTEGER       N,K(4000,5), npad
      COMMON/PYJETS/N, npad, K, P, V
      save /pyjets/
C
      INTEGER NLINES
      REAL*8 PLAB(4,9)
      REAL*4 RANPHI,SINPHI,COSPHI,RANY,RANUDQ
C  INFORMATION FOR JETSET PACKAGE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      INTEGER njoin
      PARAMETER(NJOIN=2)
      INTEGER JLPSF(NJOIN),JLVAL(NJOIN),JLSEA1(NJOIN),JLSEA2(NJOIN)
      DATA    JLPSF/10,11/
      DATA    JLVAL/2,7/
      DATA    JLSEA1/4,9/
      DATA    JLSEA2/2,3/
      LOGICAL LSEA,LVAL
      save jlpsf,jlval,jlsea1,jlsea2
*KEEP,PARTIC.
      INTEGER I2STAT(13),I2PART(13),
     &        I2MO1(13),I2DA1(13),I2DA2(13)
      COMMON /PARTIC/  I2LINE,I2STAT,I2PART,I2MO1,I2DA1,I2DA2
*KEND.
C
      real *8 pairm
      common /gmulpm/ pairm
      save /gmulpm/
C
      REAL*4  I2MASS(NLIMAX)
      DATA    I2MASS/NLIMAX*-9999.9/
C
      DATA NCALL/0/
      DATA NFRAC3/0/
      DATA NTERM3/0/
      DATA NEXTW/1/
      save ncall,nfrac3,nterm3,nextw
      real ran2
      integer idum
      data idum/-1/
C
      NFRACS=0
C
      NLINES=9
      NFINAL=0
      NINIT =0
C
C LORENZ TRANFORMATION AND COMPUTING OF INTRMED PARTICLES <======
C PARTICLE 1 = "PROTON"    <==================
      PLAB(1,1)=0.0
      PLAB(2,1)=0.0
      PLAB(3,1)=GAMMA*P1 + BETGAM*E1
      PLAB(4,1)=GAMMA*E1 + BETGAM*P1
C PARTICLE 2 = ELEKTRON IN <================
      PLAB(1,2)=0.0
      PLAB(2,2)=0.0
      PLAB(3,2)= -GAMMA*P1 + BETGAM*E2
      PLAB(4,2)=  GAMMA*E2 - BETGAM*P1
C PARTICLE 9 = ELEKTRON OUT <===============
      PLAB(1,9)=P5*ST5*CP5
      PLAB(2,9)=P5*ST5*SP5
      PLAB(3,9)=GAMMA*CT5*P5 + BETGAM  *  E5
      PLAB(4,9)=GAMMA  *  E5 + BETGAM*CT5*P5
c      print *,'electron :', p5*p5-e5*e5,gamma,betgam
c      print *,PLAB(4,9)**2-PLAB(3,9)**2-PLAB(2,9)**2-PLAB(1,9)**2
C PARTICLE 4 = GAMMA_E   <==================
      PLAB(1,4)=PLAB(1,2)-PLAB(1,9)
      PLAB(2,4)=PLAB(2,2)-PLAB(2,9)
      PLAB(3,4)=PLAB(3,2)-PLAB(3,9)
      PLAB(4,4)=PLAB(4,2)-PLAB(4,9)
C PARTICLE 5 = QUARK OUT <==================
      PLAB(1,5)=P3*ST3*CP3
      PLAB(2,5)=P3*ST3*SP3
      PLAB(3,5)=GAMMA*CT3*P3 + BETGAM  *  E3
      PLAB(4,5)=GAMMA  *  E3 + BETGAM*CT3*P3
C PARTICLE 3 = GAMMA_P   <==================
      PLAB(1,3)=PLAB(1,1)-PLAB(1,5)
      PLAB(2,3)=PLAB(2,1)-PLAB(2,5)
      PLAB(3,3)=PLAB(3,1)-PLAB(3,5)
      PLAB(4,3)=PLAB(4,1)-PLAB(4,5)
C PARTICLE 6 = MUON1    <==================
      PLAB(1,6)=P6*ST6*CP6
      PLAB(2,6)=P6*ST6*SP6
      PLAB(3,6)=GAMMA*CT6*P6 + BETGAM  *  E6
      PLAB(4,6)=GAMMA  *  E6 + BETGAM*CT6*P6
C PARTICLE 7 = MUON1->2 <==================
      PLAB(1,7)=PLAB(1,3)-PLAB(1,6)
      PLAB(2,7)=PLAB(2,3)-PLAB(2,6)
      PLAB(3,7)=PLAB(3,3)-PLAB(3,6)
      PLAB(4,7)=PLAB(4,3)-PLAB(4,6)
C PARTICLE 8 = MUON2     <==================
      PLAB(1,8)=P7*ST7*CP7
      PLAB(2,8)=P7*ST7*SP7
      PLAB(3,8)=GAMMA*CT7*P7 + BETGAM  *  E7
      PLAB(4,8)=GAMMA  *  E7 + BETGAM*CT7*P7

c----> Lepton pair mass
      pairm=sqrt(
     &           ((PLAB(4,3)+PLAB(4,4))**2-(PLAB(3,3)+PLAB(3,4))**2)
     &          -((PLAB(1,3)+PLAB(1,4))**2+(PLAB(2,3)+PLAB(2,4))**2)
     &          )
c      print *,'Leptons pair mass =',pairm,'GeV'

C====> SET KINEMATIC VARIABLES FOR GKI   <============
      GMUX= -T2 /(EP*PLAB(4,4)-PP*PLAB(3,4))/2.0D0
      GMUY= (EP*PLAB(4,4)-PP*PLAB(3,4))/
     &     (EE*PLAB(4,4)+PE*PLAB(3,4))
      GMUW = (EP+PLAB(4,4))**2 - (PP+PLAB(3,4))**2
!     gmuw =((EP-PP)+(PLAB(4,4)-PLAB(3,4)))*
!     &                  ((EP+PP)+(PLAB(4,4)+PLAB(3,4)))
      IF (GMUW .GE. 0) THEN
         GMUW=SQRT(GMUW)
      ELSE
         WRITE(6,*) ' GMUFIL : NEGATIV W**2 COMPUTED : GENW**2 = ', GMUW
     &        ,' GENW IS SET TO 0 '
         print *,'plab(1)',(plab(j,1),j=1,4)
         print *,'plab(9)',(plab(j,9),j=1,4)
         print *,'gamma',gamma,betgam
         print *,'p,e,s,c',p5,e5,st5,cp5
         GMUW=0.0
      ENDIF
      GMUNU= GMUY*2.0*ULMASS(2212)/EP/EE
C===> RANDOM REFLECTION AT XZ-PLAIN <==================
      IF (ran2(idum)  .GE. 0.5) THEN
         RANY=-1.0
      ELSE
         RANY=1.0
      ENDIF
C===> RANDOM ROTATION AT Z-AXIS <=====================
      RANPHI=PI2*ran2(idum)
      SINPHI=SIN(RANPHI)
      COSPHI=COS(RANPHI)
C====> ROTATE, REFELECT AND TRANSFORM TO REAL*4 VALUES <=============
      DO 100 I=1,9
         PL(1,I) = SNGL(PLAB(1,I))*COSPHI + RANY*SNGL(PLAB(2,I))*SINPHI
         PL(2,I) =-SNGL(PLAB(1,I))*SINPHI + RANY*SNGL(PLAB(2,I))*COSPHI
         PL(3,I) = SNGL(PLAB(3,I))
         PL(4,I) = SNGL(PLAB(4,I))
c         print *,'Particle',I,'P=',(PL(j,I),j=1,4)
 100  CONTINUE
      
c      print *,'after rotation:'
c      print *,PL(4,9)**2-PL(3,9)**2-PL(2,9)**2-PL(1,9)**2
C===> RANDOM DISTRIBUTION OF LEPTON+ AND LEPTON- <===========
      IF (RAN2(idum) .LT. 0.5) THEN
         I2PART(6) = IPAIR
         I2PART(7) =-IPAIR
         I2PART(8) =-IPAIR
      ELSE
         I2PART(6) =-IPAIR
         I2PART(7) = IPAIR
         I2PART(8) = IPAIR
      ENDIF
C===> SELECTION OF HADRON MODE IN PARTON MODEL <================
      LVAL=.FALSE.
      LSEA=.FALSE.
      IF (PMOD .EQ. 101) THEN
         LVAL=.TRUE.
      ELSEIF (PMOD .EQ. 102) THEN
         LSEA=.TRUE.
      ELSEIF (PMOD .EQ. 103) THEN
         IF (ran2(idum) .GT. PSEA) THEN
            LVAL=.TRUE.
         ELSE
            LSEA=.TRUE.
         ENDIF
      ENDIF

c      print *,'LVAL=',LVAL,', LSEA=',LSEA
C===> add. Particles for Val.Quark scatering in Parton model <=====
      IF (LVAL) THEN
         NINIT =2
C===> ADD INITIAL PROTON <====
         I2STAT(10)=21
         I2PART(10)=2212
         I2MO1(10)=0
         I2DA1(10)=2
         I2DA2(10)=3
C====> CORRECT PART 1 <===
         I2MO1(1)=10
C====> ADD DIQUARK  <====
         I2STAT(11)=1
         I2MO1(11)=1
         I2DA1(11)=0
         I2DA2(11)=0
C===> RANDOM SELECTION OF U AND D QUARKS <==
         RANUDQ=ran2(idum)
         IF (RANUDQ .LT. PVALD) THEN
            I2PART(1)=1
            I2PART(11)=2203
            I2PART(5)=1
         ELSEIF (RANUDQ .LT. 0.5+0.5*PVALD) THEN
            I2PART(1)=2
            I2PART(11)=2101
            I2PART(5)=2
         ELSE
            I2PART(1)=2
            I2PART(11)=2103
            I2PART(5)=2
         ENDIF
         IUSEDF=I2PART(1)
C==> SET MASSES <=============
         I2MASS(1)=SNGL(MQ)
         I2MASS(11)=ULMASS(I2PART(11))
         I2MASS(10)=ULMASS(2212)
C==> SET MOMENTA <============
         PL(1,10)=0.0
         PL(2,10)=0.0
         PL(3,10)=SNGL(PP)
         PL(4,10)=SNGL(EP)
         PL(1,11)=0.0
         PL(2,11)=0.0
         PL(3,11)=SNGL(PP-PQ)
         PL(4,11)=PL(4,10)-PL(4,1)
C???> computed energy may be in conflict to mass !!!!  ????????
C
C===> add. Particles for SEA Quark scatering in Parton model <=====
      ELSEIF (LSEA) THEN
         NINIT =4
C===> ADD INITIAL PROTON <====
         I2STAT(10)=21
         I2PART(10)=2212
         I2MO1(10)=0
         I2DA1(10)=2
         I2DA2(10)=5
C====> ADD PART CODE FOR ADD PARTICLES <======
         DO 234 I=11,13
            I2DA1(I)=0
            I2DA2(I)=0
            I2MO1(I)=1
            I2STAT(I)=1
 234     CONTINUE
C====> CORRECT PART 1<===
         I2MO1(1)=1
C===> SET SCATTERED QUARK AND HIS ANTI QUARK <==
         IF (NQUARK .NE. 12) THEN
            I2PART(1)=NQUARK
         ELSE
            IF (ran2(idum) .LT. 0.2) THEN
               I2PART(1)=1
            ELSE
               I2PART(1)=2
            ENDIF
         ENDIF
         IUSEDF=I2PART(1)
         IF (ran2(idum) .LE. 0.5) THEN
            I2PART(12)=-I2PART(1)
         ELSE
            I2PART(12)=I2PART(1)
            I2PART(1)=-I2PART(1)
         ENDIF
         I2PART(5)=I2PART(1)
C===> ADD QUARK AND DIQUARK FROM P <==
         IF (I2PART(1) .LT. 0) THEN
            IPQ=13
            IPDQ=11
         ELSE
            IPQ=11
            IPDQ=13
         ENDIF
         RANUDQ=ran2(idum)
         IF (RANUDQ .LT. 1.0/3.0) THEN
            I2PART(IPQ)=1
            I2PART(IPDQ)=2203
         ELSEIF (RANUDQ .LT. 2.0/3.0) THEN
            I2PART(IPQ)=2
            I2PART(IPDQ)=2101
         ELSE
            I2PART(IPQ)=2
            I2PART(IPDQ)=2103
         ENDIF
C==> SET MASSES <=============
         I2MASS(1)=SNGL(MQ)
         I2MASS(5)=SNGL(MQ)
         I2MASS(12)=SNGL(MQ)
         I2MASS(13)=ULMASS(I2PART(13))
         I2MASS(11)=ULMASS(I2PART(11))
         I2MASS(10)=ULMASS(2212)
C==> SET MOMENTA <=============
         PL(1,10)=0.0
         PL(2,10)=0.0
         PL(3,10)=SNGL(PP)
         PL(4,10)=SNGL(EP)
         PL(1,12)=0.0
         PL(2,12)=0.0
clf         PL(3,12)=GMUSELX(-IABS(IUSEDF),QSCALE)*PP
         PL(4,12)=SQRT(PL(3,12)**2+I2MASS(12)**2)
         PL(1,IPQ)=0.0
         PL(2,IPQ)=0.0
clf         PL(3,IPQ)=GMUSELX(I2PART(IPQ),QSCALE)*PP
         PL(4,IPQ)=SQRT(PL(3,IPQ)**2+I2MASS(IPQ)**2)
         PL(1,IPDQ)=0.0
         PL(2,IPDQ)=0.0
         PL(3,IPDQ)=PL(3,10)-PL(3,1)-PL(3,IPQ)-PL(3,12)
         PL(4,IPDQ)=PL(4,10)-PL(4,1)-PL(4,IPQ)-PL(4,12)
      ENDIF
C
C    FOR INELASTIC MODE WITH STRUCTURE FUNCTIONS BUILD <========
C    HADRONIC SYSTEM USING LUND SHOWER MC.
C
      IF (PMOD.GE.10 .AND. PMOD.LE.99) THEN
         NFINAL=2
C====> INSERT THE MASS OF THE HADRONIC SYSTEM <==================
         I2MASS(5)=SNGL(MX)
C===> RANDOM SELECTION OF U , D AND DI QUARKS <===========
         RANUDQ=ran2(idum)
         IF (RANUDQ .LT. 1.0/9.0) THEN
            I2PART(10)=1
            I2PART(11)=2203
         ELSEIF (RANUDQ .LT. 5.0/9.0) THEN
            I2PART(10)=2
            I2PART(11)=2101
         ELSE
            I2PART(10)=2
            I2PART(11)=2103
         ENDIF
         ULMDQ=ULMASS(I2PART(11))
         ULMQ =ULMASS(I2PART(10))
C====> SET OF LUND CODES <====================================
         I2MO1(10)=5
         I2DA1(10)=0
         I2DA2(10)=0
         I2STAT(10)=1
         I2MO1(11)=5
         I2DA1(11)=0
         I2DA2(11)=0
         I2STAT(11)=1
C====> CHOOSE RANDOM DIRECTION IN MX FRAME <===================
         RANMXP=PI2*ran2(idum)
         RANMXT=ACOS(2.0*ran2(idum)-1.0)
C====> COMPUTE MOMENTUM OF DECAY PARTICLE FROM MX <=============
         PMXP=(MX**2-ULMDQ**2+ULMQ**2)**2/4.0/MX/MX - ULMQ**2
         if (pmxp.lt.0) return !FIXME FIXME FIXME FIXME !!!!!!!!
         pmxp=dsqrt(pmxp)
c         print *,ulmdq,ulmq,mx,pmxp,
c     +        (MX**2-ULMDQ**2+ULMQ**2)**2/4.0/MX/MX-ULMQ**2
C=====> BUILD 4-VECTORS AND BOOST DECAY PARTICLES <===============
         PMXDA(1)=SIN(RANMXT)*COS(RANMXP)*PMXP
         PMXDA(2)=SIN(RANMXT)*SIN(RANMXP)*PMXP
         PMXDA(3)=COS(RANMXT)*PMXP
         PMXDA(4)=SQRT(PMXP**2+ULMDQ**2)
c         PRINT *,' GMUFIL : PMXDA BEFORE LB:',(PMXDA(I),I=1,4)
         CALL LORENB(I2MASS(5),PL(1,5),PMXDA(1),PL(1,11))
c         PRINT *,' GMUFIL : PL(11) AFTER LB:',(PL(I,11),I=1,4)
         PMXDA(1)=-PMXDA(1)
         PMXDA(2)=-PMXDA(2)
         PMXDA(3)=-PMXDA(3)
         PMXDA(4)=SQRT(PMXP**2+ULMQ**2)
c         PRINT *,' GMUFIL : PMXDA BEFORE LB:',(PMXDA(I),I=1,4)
         CALL LORENB(I2MASS(5),PL(1,5),PMXDA(1),PL(1,10))
      ENDIF
C====> PREPARE THE LUND COMMON <================================
 10   CONTINUE
c      print *,'before lunset, NLINES=',NLINES,', NINIT=',NINIT,
c     +     ', NFINAL=',NFINAL
      CALL LUNSET(NLINES+NINIT+NFINAL)
C ====> FILLING THE LUND COMMON <================================
      DO 200 I=1+NLINES,NINIT+NLINES
C SET MOTHER/DAUGHTER VALUES, MARKING PARTICLES AS DECAYED <=======
         CALL LUKSET(I-NLINES,I2STAT(I),I2PART(I),
     &        I2MO1(I),I2DA1(I),I2DA2(I),0)
C SET PULS, ENERGY AND MASS OFF THE PARTICLES <==================
         CALL LUPSET(I-NLINES,PL(1,I),PL(2,I),PL(3,I),PL(4,I),I2MASS(I))
 200  CONTINUE
      DO 201 I=1,NLINES+NFINAL
C SET MOTHER/DAUGHTER VALUES, MARKING PARTICLES AS DECAYED <=======
         CALL LUKSET(I+NINIT,I2STAT(I),I2PART(I),
     &        I2MO1(I),I2DA1(I),I2DA2(I),NINIT)
C SET PULS, ENERGY AND MASS OF THE PARTICLES <==================
         CALL LUPSET(I+NINIT,PL(1,I),PL(2,I),PL(3,I),PL(4,I),I2MASS(I))
c         IF(I.EQ.9) THEN
c            PRINT *,I,I2MASS(I),(PL(J,I),J=1,4)
c         ENDIF
 201  CONTINUE
C PUTTING QUARK AND DIQUARK TO A COLOR SINGLET <======================
      IF (LVAL) CALL LUJOIN(NJOIN,JLVAL)
      IF (LSEA) THEN
         CALL LUJOIN(NJOIN,JLSEA1)
         CALL LUJOIN(NJOIN,JLSEA2)
      ENDIF
c      DO 1024, J=1,11
c         print *,'I=',J,'STATUS=',I2STAT(J)
c 1024 CONTINUE
c      print *,'Before LUJOIN================================'
c      CALL LULIST(2)
      IF (PMOD.GE.10 .AND. PMOD.LE.99) CALL LUJOIN(NJOIN,JLPSF)
C EXECUTE LUND FRAGMENTATION PROGRAM  <==============================
c      print *,'Before LUEXEC================================'
c      call LULIST(2)
      CALL LUEXEC
c      print *,'After  LUEXEC================================'
clf      CALL LUEXEC
C Check wether the Hadronic system is inelastic  <===================
      IF (PMOD.GE.10 .AND. PMOD.LE.99) THEN
         NPOUT=0
         DO 300 I=1,N
            IF (K(I,1) .EQ. 1) NPOUT=NPOUT+1 ! History code (KH=1 = mother particle)
 300     CONTINUE
         NFRAC3=NFRAC3+1
         NFRACS=NFRACS+1
         IF (NPOUT .EQ. 4 .AND. NFRACS .LE. 1000) GOTO 10
         IF (NFRACS .GT. 1000) NTERM3=NTERM3+1
      ENDIF
      NCALL=NCALL+1
c      CALL LULIST(1)
      IF (NCALL .GE. NEXTW) THEN
         IF (PMOD.GE.10 .AND. PMOD.LE.99) THEN
            WRITE(6,*) ' GMUFIL : NUMBER OF CALLS IS ',NCALL
     &           ,'  PMOD 10-99:  # FRAC TRY :',NFRAC3,'  # FRAC TERM :'
     &           ,NTERM3
         ELSE
c            WRITE(6,*) ' GMUFIL : NUMBER OF CALLS IS ',NCALL,' W',GMUW,
c     &           pairm
         ENDIF
         CALL LULIST(2)
         NEXTW=NEXTW*2
      ENDIF
!-      CALL LUHEPC(1)
      END
