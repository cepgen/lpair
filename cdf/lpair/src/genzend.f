       subroutine genzend
*
* --- LPAIR data common block
*
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*
       IF(IPAR(3).EQ.0) RETURN
*
* ---- end GENZ
       CALL GNZEND
*
       return
       end
