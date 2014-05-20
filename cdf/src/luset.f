*--   Author :    O. Duenger   14/11/91
*     
      SUBROUTINE LUKSET(LINE,STATUS,PART,MOTH,DAUG1,DAUG2,NOFF)
*=================================================================
      IMPLICIT NONE
      INTEGER LINE,STATUS,PART,MOTH,DAUG1,DAUG2,NOFF
*
C---JETSET and GENOUT common
      REAL*4        P(4000,5),V(4000,5)
      INTEGER       N,K(4000,5)
      COMMON/LUJETS/N,K,P,V
*     
      IF (LINE .GT. N) THEN
         WRITE(6,*)' LUKSET : LINE TOO BIG; LINE=',LINE,' N=',N
         RETURN
      ENDIF
*     
      IF (STATUS .NE. -9999) K(LINE,1)=STATUS
      IF (PART   .NE. -9999) K(LINE,2)=PART
      IF (MOTH. NE.-9999 .AND. MOTH .NE.0) K(LINE,3)=MOTH +NOFF
      IF (DAUG1.NE.-9999 .AND. DAUG1.NE.0) K(LINE,4)=DAUG1+NOFF
      IF (DAUG2.NE.-9999 .AND. DAUG2.NE.0) K(LINE,5)=DAUG2+NOFF
      END


*--   Author :    O. Duenger   19/12/91
*     
      SUBROUTINE LUNSET(LINE)
*=================================================================
      IMPLICIT NONE
      INTEGER LINE,I,II
*     
C---JETSET and GENOUT common
      REAL*4        P(4000,5),V(4000,5)
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

*--   Author :    O. Duenger   19/12/91
*     
      SUBROUTINE LUPSET(LINE,PX,PY,PZ,E,M)
*=================================================================
      IMPLICIT NONE
      INTEGER LINE
      REAL PX,PY,PZ,E,M,ULMASS
*
C---JETSET and GENOUT common
      REAL*4        P(4000,5),V(4000,5)
      INTEGER       N,K(4000,5)
      COMMON/LUJETS/N,K,P,V
*     
      IF (LINE .GT. N) THEN
         WRITE(6,*) ' LUPSET : TOO BIG LINE NUMBER, LINE =',LINE,
     +        ', N =',N
         RETURN
      ENDIF
*     
      P(LINE,1)=PX
      P(LINE,2)=PY
      P(LINE,3)=PZ
      P(LINE,4)=E
      IF (M .GE. -9998.0) THEN
         P(LINE,5)=M
      ELSE
         P(LINE,5)=ULMASS(K(LINE,2))
      ENDIF
      END
