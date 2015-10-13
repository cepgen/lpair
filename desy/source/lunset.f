*-- Author :    O. Duenger   19/12/91
*
      SUBROUTINE LUNSET(LINE)
*=================================================================
      IMPLICIT NONE
      INTEGER LINE,I,II
*
C---JETSET and GENOUT common
      REAL P(4000,5),V(4000,5)
      INTEGER       N,K(4000,5)
      COMMON/LUJETS/N,K,P,V
*
      IF ((LINE .LT. 1) .OR. (LINE .GT. 4000)) THEN
         WRITE(6,*) ' LUNSET : WRONG LINE, LINE =',LINE
         RETURN
      ENDIF
*     
      N=LINE
*     
      DO 100 I=1,5
         DO 200 II=1,N
            V(II,I)=0.0
 200     CONTINUE
 100  CONTINUE
      END
