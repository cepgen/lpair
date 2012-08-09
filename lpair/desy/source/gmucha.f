*-- Author :    ZEUS Offline Group   18/08/94

      SUBROUTINE GMUCHA

***********************************************************************
*
*    SUBROUTINE GMUCHA
*
* PURPOSE: Interpret the data cards to change the default parameters of MUPAIR.
*
* INPUT: Bos text bank GMUP
*
* OUTPUT: updated parameters in MUPAIR commons.
*
* CALLED BY: GMUPA
*
* AUTHOR: OLAF DUENGER             CREATED AT: 91/12/12
*
* CHANGED BY:                              AT:
* REASON:
*
************************************************************************
*
        Implicit NONE
*
*KEEP,VEGPAR.
      INTEGER          NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND
      COMMON /VEGPAR/  NDIM,NCVG,ITMX,NPRN,IGRAPH,
     &                 NPOIN,NPRIN,NTREAT,IBEG,IEND

*KEEP,BEAM.
      INTEGER          INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,IPAIR,NQUARK
      REAL*8           INPE,INPP
      COMMON /BEAM/    INPE,INPP,INTGE,INTGP,GPDF,SPDF,PMOD,EMOD,
     &                 IPAIR,NQUARK

*KEEP,CUTS.
      INTEGER      MODCUT
      REAL*4       THMAX,THMIN,MXMN,MXMX,Q2MN,Q2MX
      REAL*8       COTTH1,COTTH2,ECUT,PTCUT,MXMIN2,MXMAX2,QP2MIN,QP2MAX
      COMMON /CUTS/COTTH1,COTTH2,ECUT,PTCUT,MXMIN2,MXMAX2,
     &             THMAX,THMIN,QP2MIN,QP2MAX,MODCUT,MXMN,MXMX,Q2MN,Q2MX

*KEND.
*
      Real         RINPP, RINPE, RECUT, RPTCUT
        Common/DumKL/RINPP, RINPE, RECUT, RPTCUT
*
C*  End of common
*
        Integer IErr
*
*--------  Read data cards and overwrite defaults:
*
*..  Reinitialise FFRead:
        Call KWFFRD(IErr)
*
*..  Define Key Words:
        Call FFKEY ('IBEG',   IBEG,    1,    'NONE')
        Call FFKEY ('IEND',   IEND,    1,    'NONE')
        Call FFKEY ('NTRT', NTReat,    1,    'NONE')
        Call FFKEY ('PRVG',  NPrin,    1,    'NONE')
        Call FFKEY ('NCVG',   NCVG,    1,    'NONE')
        Call FFKEY ('ITVG',   ITMX,    1,    'NONE')
        Call FFKEY ('NCSG',  NPOIN,    1,    'NONE')
        Call FFKEY ('INPP',  RINPP,    1,    'NONE')
        Call FFKEY ('PMOD',   PMOD,    1,    'NONE')
        Call FFKEY ('GPDF',   GPDF,    1,    'NONE')
        Call FFKEY ('SPDF',   SPDF,    1,    'NONE')
        Call FFKEY ('INPE',  RINPE,    1,    'NONE')
        Call FFKEY ('EMOD',   EMOD,    1,    'NONE')
        Call FFKEY ('PAIR',  IPAIR,    1,    'NONE')
        Call FFKEY ('QPDF', NQUARK,    1,    'NONE')
        Call FFKEY ('MCUT', MODCUT,    1,    'NONE')
        Call FFKEY ('THMX',  THMAX,    1,    'NONE')
        Call FFKEY ('THMN',  THMIN,    1,    'NONE')
        Call FFKEY ('ECUT',  RECUT,    1,    'NONE')
        Call FFKEY ('PTCT', RPTCUT,    1,    'NONE')
        Call FFKEY ('Q2MN',   Q2MN,    1,    'NONE')
        Call FFKEY ('Q2MX',   Q2MX,    1,    'NONE')
        Call FFKEY ('MXMN',   MXMN,    1,    'NONE')
        Call FFKEY ('MXMX',   MXMX,    1,    'NONE')
*
*..  Read cards and fill commons:
        Call KWFFGO('ZLPAIR', IErr)
*
*..  Fix up double precision entries:
         INPP = RINPP
         INPE = RINPE
         ECUT = RECUT
        PTCUT = RPTCUT

*
*--------  Return
*
        Return
        End
