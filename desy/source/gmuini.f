*-- Author :    ZEUS Offline Group   18/08/94
      SUBROUTINE GMUINI

      Implicit NONE

C ####################################################################
C #   SET DEFAULT PARAMETERS OF THE GENERATOR.                       #
C ####################################################################

      Real         PI
      PARAMETER   (PI=3.141592654)

*KEEP,VEGPAR.
      INTEGER          NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND,NGEN
      COMMON /VEGPAR/  NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND,NGEN
      save /vegpar/

*KEEP,BEAM.
      INTEGER          INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,IPAIR,NQUARK
      REAL*8           INPE,INPP
      COMMON /BEAM/    INPE,INPP,INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,
     &                 IPAIR,NQUARK
      save /beam/

*KEEP,CUTS.
      INTEGER      MODCUT
      REAL*4       THMAX,THMIN,MXMN,MXMX,Q2MN,Q2MX
      REAL*8       COTTH1,COTTH2,ECUT,PTCUTMIN,PTCUTMAX,MXMIN2,MXMAX2,
     &             QP2MIN,QP2MAX
      COMMON /CUTS/COTTH1,COTTH2,ECUT,PTCUTMIN,PTCUTMAX,MXMIN2,MXMAX2,
     &             THMAX,THMIN,QP2MIN,QP2MAX,MODCUT,MXMN,MXMX,Q2MN,Q2MX
      save /cuts/
*KEND.
*
*  PULS OF INCOMING PROTON AND ELECTRON
         INPE   = 26.7D0
         INPP   = 820.0D0
*
*  MODE OF INCOMING PROTON AND ELECTRON
         EMOD   = 1
         PMOD   = 2
         NQUARK = 12
*
*  CODE OF THE PRODUCED LEPTON PAIR - Electron pair default
         IPAIR  = 11
*
*  START AND STOP POINT OF THE GENERATOR PROGRAMM (1-3 = WHOLE RUN)
         IBEG   = 1
         IEND   = 3
         NGEN   = 100
*
*  NUMBER OF CALLS PER VEGAS ITERRATION
         NCVG   = 14000
*
*  NUMBER OF VEGAS ITERRATIONS
         ITMX   = 10
*
*  VEGAS PRINT PARAMETER
         NPRN   = 1
         IGRAPH = 0
*
*  NUMBER CALL PER BIN IN SETGEN (NR. OF BINS IS 3**8=6561 OR 3**7=2187)
         NPOIN  = 100
*
*  VEGAS PRINT FLAG AND STRATEGY NUMBER
         NPRIN  = 1
         NTREAT = 1
         NGEN = 10000
*
C MODCUT : MODE FOR CUT 0=NO; 1=VERMASEREN DET; 2=GIVEN PARAMETER <==
         MODCUT = 2
C THMIN,THMAX : MIN AND MAX THETA OF BOTH MUONS <====================
         THMIN  = 5.0
         THMAX  = 175.0
C ECUT : MIN. ENERGY OF BOTH MUONS <============================
         ECUT   = 1.0D0
C PTCUT : MIN. AND MAX. TRANSVERS MOMENT OF BOTH MUONS <==================
         PTCUTMIN = 0.5D0
         PTCUTMAX = 1.0D5
C GPDF, SPDF : DEFAULT PDF GRV LO <======
         GPDF = 5
         SPDF = 4
C MXMN, MXMX : GIVE THE LIMITS FOR THE MASS IM FINAL HADRONIC SYSTEM
         MXMN = 1.070
         MXMX = 320.0
C Q2MN, Q2MX : GIVE THE LIMITS FOR ABS(Q**2) AT THE PROTON SIDE
         Q2MN = 0.0
         Q2MX = 1E5
      RETURN
      END
