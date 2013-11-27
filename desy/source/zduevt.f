*-- Author :
      SUBROUTINE ZDUEVT(IWANT)
*     ========================
*
*------------------------------------------------------------------------
*
*  ZDUEVT: Adminstrate event generation by LPair.
*  =======
*
*------------------------------------------------------------------------
*
      Implicit NONE
*     
      REAL*8          S1,S2,T1,T2
      COMMON /EXTRA/  S1,S2,T1,T2
*     
      Integer IWant
*     
*
*--------  Initialise
*
*..  Want the event is the default:
        IWant = 1
*
*--------  Generate the event:
*
*..  Generate the event with LPair:
        Call GMUGNA
*
*..  Fill /LUJETS/
        Call GMUFIL
        Call GMULHE
*
        End
