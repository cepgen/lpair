      SUBROUTINE terminate
* Subroutine finishing the job, and closing all opened
* tasks...
* Author: L.Forthomme, Oct 2016
*
      IMPLICIT none

      INTEGER nout,ilhef
      COMMON/outp/nout,ilhef

      INTEGER ipar(20)
      DOUBLE PRECISION lpar(20)
      COMMON/datapar/ipar,lpar

      IF(ipar(2).eq.2) THEN
        CALL LHEEND
      ENDIF

      END

