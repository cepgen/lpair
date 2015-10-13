*-- Author :    O. Duenger   19/12/91
*
      SUBROUTINE LUPSET(LINE,PX,PY,PZ,E,M)
*=================================================================
      IMPLICIT NONE
      INTEGER LINE
      REAL*4 PX,PY,PZ,E,M
      REAL ulmass
*
C---JETSET and GENOUT common
      REAL P(4000,5),V(4000,5)
      INTEGER          N,K(4000,5)
      COMMON/LUJETS/N,K,P,V
*     
      IF (LINE .GT. N) THEN
         WRITE(6,*) ' LUPSET : TOO BIG LINE NUMBER, LINE =',LINE,', N =',N
         RETURN
      ENDIF
*
      P(LINE,1)=PX
      P(LINE,2)=PY
      P(LINE,3)=PZ
      P(LINE,4)=E
      IF (M .GE. -9998.0) THEN
         P(LINE,5)=M
c     print *,'--> Setting the mass for the particle',LINE
c     print *,'    E**2-P**2=',E**2-PX**2-PY**2-PZ**2
c     print *,'    M**2=',M**2
      ELSE
         P(LINE,5)=ULMASS(K(LINE,2))
      ENDIF
      END
