*-- Author :    ZEUS Offline Group   18/08/94
      FUNCTION GMUSELX(IFLAV,DQ2)
* SELECT MOMENTUM FRACTION USING REJECTION METHOD AND PDFLIB
*
************************************************************
      REAL*4 H1RN
      REAL*8          XMIN,XMAX,Q2MIN,Q2MAX
      COMMON /W50513/ XMIN,XMAX,Q2MIN,Q2MAX
      REAL*8 TMAS
      INTEGER       PDFMOD,PDFNFL,PDFLO,PDFDUM
      COMMON/W50511/PDFMOD,PDFNFL,PDFLO,PDFDUM,TMAS
*
      REAL*8 X,QSCALE,DENS(-6:2),DQ2
      REAL*8 XDMAX(-6:2,10),QSDMAX(10),QSMAX,QSMIN,HXDMAX
      REAL*8 XX,DX
      DATA XDMAX/90*0.0D0/
      DATA NCALL,NPDFC,NWRITE/0,0,1/
*
      IF (IFLAV.GT.2 .OR. IFLAV.LT.PDFNFL) THEN
        WRITE(6,*) 'GMUSELX : WRONG FLAVOUR GIVEN; IFLAV =',IFLAV
        RETURN
      ENDIF
*
      IF (NCALL.LE.0) THEN
        QSMIN=DSQRT(Q2MIN)
        QSMAX=DSQRT(Q2MAX)
        DO 100 I=1,10
          QSCALE=QSMIN*(QSMAX/QSMIN)**((I-1)/9.0D0)
          QSDMAX(I)=QSCALE
          DO 200 II=1,1000
            CALL GMURANX(IHELP,XX,DX)
            CALL PDF2PDG(XX,QSCALE,DENS)
            DO 300 III=-6,2
              IF (DENS(III)/XX*DX.GT.XDMAX(III,I))
     &          XDMAX(III,I)=DENS(III)/XX*DX
300         CONTINUE
200       CONTINUE
100     CONTINUE
      ENDIF
      NCALL=NCALL+1
      IF (NCALL .GE. NWRITE) THEN
        NWRITE=NWRITE*2
        WRITE(6,*) 'GMUSELX : NCALL=',NCALL,' # OF PDF CALLS =',NPDFC
      ENDIF
      QSCALE=DSQRT(DQ2)
      HXDMAX=0.0D0
      DO 101 I=1,9
        IF (QSDMAX(I).LE.QSCALE .AND. QSDMAX(I+1).GT.QSCALE) THEN
          IF (QSDMAX(I) .GT. QSDMAX(I+1)) THEN
            IQS=I
            HXDMAX=QSDMAX(I)
          ELSE
            IQS=I+1
            HXDMAX=QSDMAX(I+1)
          ENDIF
        ENDIF
101   CONTINUE
      IF (HXDMAX .EQ. 0.0D0) THEN
        WRITE(6,*) 'GMUSELX : QSCALE OUT OF RANGE (Q2,Q2MIN,Q2MAX) =',
     &    DQ2,Q2MIN,Q2MAX
        GOTO 9999
      ENDIF
      DO 102 I=1,10000
        CALL GMURANX(IFLAV,X,DX)
        CALL PDF2PDG(X,QSCALE,DENS)
        NPDFC=NPDFC+1
        IF (DENS(IFLAV)*DX/X.GE.XDMAX(IFLAV,IQS)*H1RN(DUMMY)) THEN
          GMUSELX=X
          RETURN
        ENDIF
102   CONTINUE
      WRITE(6,*)
     &  'GMUSELX : # CYCLES REACHES 10000 : AS RETURN X=.001 IS TAKEN',
     &   ' FLAVOUR =',IFLAV,' Q2 =',DQ2
 9999 CONTINUE
      GMUSELX=.001
      END
*****************************************************
      SUBROUTINE GMURANX(IFLAV,X,DX)
*  MAPPING OF X
*****************************************************
      REAL*4 H1RN
      REAL*8 X,DX,Y
      REAL*8          XMIN,XMAX,Q2MIN,Q2MAX
      COMMON /W50513/ XMIN,XMAX,Q2MIN,Q2MAX
        Y=XMAX/XMIN
        X=XMIN*Y**H1RN(DUMMY)
        DX=X*DLOG(Y)
      END
***********************************************************************
      SUBROUTINE PDF2PDG(X,QSCALE,DENS)
*  CONVERT STRUCTF OUTPUT TO AN ARREAY USING PDG CODE AS INDEX
***********************************************************************
      REAL*8 DENS(-6:2)
      REAL*8 X,QSCALE
      CALL STRUCTF(X,QSCALE,DENS(2),DENS(1),DENS(-1),DENS(-3),
     &             DENS(-4),DENS(-5),DENS(-6),DENS(0))
      DENS(-2)=DENS(-1)
      END
