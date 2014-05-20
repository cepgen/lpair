      SUBROUTINE ACCEPT(n)
C
C     Save events
C
      implicit double precision (a-h,o-z)
*
      integer n
*
      CALL pawfil1              ! fill PAW ntuple/histograms
      CALL genzfil              ! fill GENZ output
*
      RETURN

      END
