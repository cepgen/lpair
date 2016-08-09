*-- Author :    O. Duenger   17/12/91
      DOUBLE PRECISION FUNCTION TREAT(F,X,NDIM)
C
C  AUTHOR      : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      DIMENSION X(10),Z(10),XIN(10)
      COMMON/VGB2/NDO,IT,SI,SI2,SWGT,SCHI,XI(50,10),SCALLS
     + ,D(50,10),DI(50,10)
      save /vgb2/
      COMMON/TREATB/W,VALTREAT,XIN,Z
      COMMON/VGASIO/NINP,NOUTP,NSGOUT
C
C
      DATA NCALL/0/
      save ncall,r
C
      IF(NCALL.EQ.0)THEN
         NCALL=1
         R=NDO
         R=R**NDIM
      ENDIF
      W=R
      I=1
      DO 4 I=1,NDIM
         XX=X(I)*NDO
         J=XX
         JJ=J+1
         Y=XX-J
c         PRINT *,I,XX,J,JJ,Y
         IF(J.LE.0)THEN
            DD=XI(1,I)
         ELSE
            DD=XI(JJ,I)-XI(J,I)
         ENDIF
         Z(I)=XI(JJ,I)-DD*(1.-Y)
c         print *,DD,z(i),y,w
         XIN(I)=X(I)
         W=W*DD
4     CONTINUE
c      PRINT *,W,F(Z)
c      if (F(Z).lt.0) print *,F(Z)
      TREAT=W*F(Z)
      VALTREAT=TREAT
C
      RETURN
      END
