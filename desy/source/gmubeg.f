*-- Author :    ZEUS Offline Group   18/08/94
      SUBROUTINE GMUBEG


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---->
C---- OUR MODIFICATION / ADDITION
C----
*KEEP,INPU.
      REAL*8       ME,MU,MP,MX,S,SQ,PE,PP,EE,EP,CONST,PI
      COMMON /INPU/ME,MU,MP,MX,S,SQ,PE,PP,EE,EP,CONST,PI

*KEEP,BEAM.
      INTEGER          INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,IPAIR,NQUARK
      REAL*8           INPE,INPP
      COMMON /BEAM/    INPE,INPP,INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,
     &                 IPAIR,NQUARK

*KEEP,VEGPAR.
      INTEGER          NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND,NGEN
      COMMON /VEGPAR/  NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND,NGEN

*KEEP,CUTS.
      INTEGER      MODCUT
      REAL*4       THMAX,THMIN,MXMN,MXMX,Q2MN,Q2MX
      REAL*8       COTTH1,COTTH2,ECUT,PTCUTMIN,PTCUTMAX,MXMIN2,MXMAX2,
     &             QP2MIN,QP2MAX
      COMMON /CUTS/COTTH1,COTTH2,ECUT,PTCUTMIN,PTCUTMAX,MXMIN2,MXMAX2,
     &             THMAX,THMIN,QP2MIN,QP2MAX,MODCUT,MXMN,MXMX,Q2MN,Q2MX

*KEEP,PARTIC.
      INTEGER I2STAT(13),I2PART(13),
     &        I2MO1(13),I2DA1(13),I2DA2(13)
      COMMON /PARTIC/  I2LINE,I2STAT,I2PART,I2MO1,I2DA1,I2DA2
*KEEP,COMGNA.
      INTEGER NGNA
      COMMON /COMMUP/  NGNA
*KEND.
C PARAMETER FOR PDFLIB <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      REAL*8           XMIN,XMAX,Q2MIN,Q2MAX
      COMMON /W50513/  XMIN,XMAX,Q2MIN,Q2MAX

C----- MGVH  27 APR 96 -------
      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VAL(20)
C-----------------------------

      EXTERNAL F
      double precision x(10),test

      REAL ULMASS

      DATA BCC/0.1D-03/

      INTEGER M2STAT(9)
      DATA    M2STAT/  21, 21, 21, 21,  1,  1, 21,  1,  1/
      INTEGER M2PART(9)
      DATA    M2PART/2212, 11, 22,22,2212,  0,  0,  0, 11/
      INTEGER M2MO1(9)
      DATA    M2MO1/    0,  0,  1,  2,  1,  3,  3,  7,  2/
      INTEGER M2DA1(9)
      DATA    M2DA1/    3,  4,  6,  0,  0,  0,  8,  0,  0/
      INTEGER M2DA2(9)
      DATA    M2DA2/    5,  9,  7,  0,  0,  0,  0,  0,  0/

c---- Lund common for the masses
      double precision kchg,pmas,parf,vckm
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      
C-----
      QP2MIN=-DBLE(Q2MN)
      QP2MAX=-DBLE(Q2MX)
C-----
      IF (PMOD .EQ. 1) THEN
         NDIM=7
         INTGP=1
         M2PART(1)=-11
         M2PART(5)=-11
      ELSEIF (PMOD .EQ. 2) THEN
         NDIM=7
         INTGP=2
      ELSEIF (PMOD .EQ. 11 .OR. PMOD .EQ. 12 .OR. PMOD .EQ. 13) THEN
         NDIM=8
         INTGP=PMOD-8
         MXMIN2=DBLE(MXMN)**2
         MXMAX2=DBLE(MXMX)**2
         M2STAT(5)=21
      ELSEIF (PMOD.GE.101 .AND. PMOD.LE.103) THEN
         NDIM=8
         INTGP=1
         MXMIN2=DBLE(MXMN)**2
         MXMAX2=DBLE(MXMX)**2
      ELSE
         WRITE(6,*) ' GMUBEG : WRONG PROTON MODE ; PMOD =',PMOD
         STOP
      ENDIF
C-----
      IF (EMOD .EQ. 1) THEN
         INTGE=1
      ELSEIF (EMOD .EQ. 2) THEN
         INTGE=2
         M2PART(2)=2212
         M2PART(9)=2212
      ELSE
         WRITE(6,*) ' GMUBEG : WRONG ELEKTRON MODE ; EMOD =',EMOD
         STOP
         ENDIF
C-----
      DO 100 I=1,9
         I2STAT(I)=M2STAT(I)
         I2PART(I)=M2PART(I)
         I2MO1(I)=M2MO1(I)
         I2DA1(I)=M2DA1(I)
         I2DA2(I)=M2DA2(I)
 100  CONTINUE
C     
C-----
      IF (.NOT.(IPAIR.EQ.11 .OR. IPAIR.EQ.13 .OR. IPAIR.EQ.15)) THEN
         WRITE(6,*) 'WRONG CODE FOR LEPTON PAIR,  PAIR=',IPAIR,
     &        '  VALID ARE 11, 13 AND 15 '
         STOP
      ENDIF
C-----
      NGNA=0
C-----
      PI    = DACOS(-1.D+00)
      MU    = ULMASS(IPAIR)
c      WRITE(6,*) '*** THE LEPTON PAIR CODE IS ',IPAIR,
c     &     '  THE MASS IS ',MU,' ****************'
      MP=ULMASS(I2PART(1))
      ME=ULMASS(I2PART(2))
c      print *,mu,mp,me
C-----
      PE    = INPE
      PP    = INPP
      EE    = DSQRT(PE*PE + ME*ME)
      EP    = DSQRT(PP*PP + MP*MP)
      S     = MP*MP  + ME*ME + 2. * (EE*EP + PE*PP)
      SQ    = DSQRT(S)
      PTOT  = PP - PE
      ETOT  = EP + EE
      CONST = (19.732D+03)**2
C COMPUTING OF COTAN FOR THETA CUTS <====================
      COTTH1= 1D0 / DTAN(DBLE(THMAX)*PI/180D0)
      COTTH2= 1D0 / DTAN(DBLE(THMIN)*PI/180D0)
c      print *,'GMUBEG: ETA in range [',
c     1     -dlog(dtan(THMAX*PI/180.d0/2.d0)),':',
c     2     -dlog(dtan(THMIN*PI/180.d0/2.d0)),']'
C SETTING OF THE PROTON STRUCTURE FUNCTION <====================
      IF (PMOD .GE. 100) THEN
C----- MGVH  27 APR 96 -------
Cxx       CALL PDFSET('MODE',MPDF)
         parm(1)= 'NPTYPE'
         VAL(1)  = 1
         PARM(2) = 'NGROUP'
         VAL(2)  = gpdf
         PARM(3) = 'NSET'
         VAL(3)  = spdf
clf         CALL PDFSET(PARM,VAL)
C-----------------------------
         WRITE(6,'(/A,/A,I8,A,I8,A,
     &  /A,G13.4,A,G13.4,A,/A,G13.4,A,G13.4,A,/A)')
     &'***************************************************************',
     &'* PDFLIB GROUP :',gpdf    ,'  SET :',spdf    ,'               *',
     &'*  XMIN = ', XMIN   ,'       XMAX = ', XMAX     ,'            *',
     &'* Q2MIN = ',Q2MIN   ,'      Q2MAX = ',Q2MAX     ,'            *',
     &'***************************************************************'
      ENDIF
      IF (Q2MAX .EQ. 0.D0) Q2MAX=1.D20
c=====
c      print *,'====================='
c      print *,'Masses',MP,ME,' Interaction',INTGP,INTGE
c      print *,'====================='
c=====
C-----
      IF (IBEG .EQ. 1) THEN
c         PRINT *,'GMUBEG : ===> VEGAS  IS OPERATIVE... '
         CALL  VEGAS(F,BCC,NDIM,NCVG,ITMX,NPRN,IGRAPH)
c         WRITE(6,*)
c     &        'SORRY SAVING OF VEGAS-PARAMETER IS ONLY POSSIBLE'//
c     &        ' ON THE IBM'
      ELSE
         WRITE(6,*)
     &        'SORRY READING OF VEGAS-PARAMETER IS ONLY POSSIBLE'//
     &        ' ON THE IBM'
         STOP
      ENDIF
C-----
      IF (IBEG .GT. 2) THEN
         PRINT *,'GMUBEG : ===> STORED SETGEN VARIABLES READ ...'
         CALL RESTR2(NDIM)
         CALL GMUPSG
      ELSEIF (IEND .LT. 2) THEN
c         WRITE(6,*)
c     &        ' GMUBEG : PROGRAM STOPS WITHOUT SETGEN FILE AND'//
c     &        ' EVENT GENERATION , IEND < 2 '
c         STOP
      ELSE
         PRINT *,'GMUBEG : ===> SETGEN IS OPERATIVE... '
         CALL SETGEN(F,NDIM,NPOIN,NPRIN,NTREAT)
c         WRITE(6,*)
c     &        'SORRY SAVING OF SETGEN-PARAMETER IS ONLY POSSIBLE'//
c     &        ' ON THE IBM'
      ENDIF
C-----
      IF (IEND .LT. 3) THEN
c         WRITE(6,*)
c     &        ' GMUBEG : PROGRAM STOPS WITHOUT EVENT GENERATION ,'//
c     &        ' IEND < 3 '
c         STOP
      ENDIF
C-----
      RETURN
      END
