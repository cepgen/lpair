       subroutine pawend
*
* --- LPAIR data common block
*
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*
       INTEGER ICYCLE
*
       IF(IPAR(2).EQ.0) RETURN
*
* ---- close histograms
*
       CALL HCDIR('//HISTO',' ')
       CALL HROUT(0,ICYCLE,' ')
       CALL HREND('HISTO')
       CLOSE(UNIT=31)
*
       return
       end
