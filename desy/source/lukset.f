*-- Author :    O. Duenger   14/11/91
*
      SUBROUTINE LUKSET(LINE,STATUS,PART,MOTH,DAUG1,DAUG2,NOFF)
*=================================================================
      IMPLICIT NONE
      INTEGER LINE,STATUS,PART,MOTH,DAUG1,DAUG2,NOFF
*
C---JETSET and GENOUT common
      REAL P(4000,5),V(4000,5)
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
