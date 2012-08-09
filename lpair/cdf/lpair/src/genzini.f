       subroutine genzini
*
*  assume 11 words/particle * max particles + 1k overhead
*
       INTEGER MEMGNZ
       PARAMETER (MEMGNZ=268000)
       REAL GJUNK
       COMMON/GNCSTO/GJUNK(MEMGNZ)
       INTEGER IE,IRET
*
* --- LPAIR data common block
*
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*
       IF(IPAR(3).EQ.0) RETURN
*
* ---- BEGIN GENZ --- (MAIN)
*
*   set GENZ parameters: run number and generator name
       CALL GNZPARI(1,'RUNN',33)
       CALL GNZPARC(1,'CFZO','XOL')
       CALL GNZPARC(1,'GENE','LPR')
*   print current values of GENZ parameters
       CALL GNZPPAR
*
*   initialize GENZ. Call with negative argument to prevent
*   re-initialization of ZEBRA
*
       CALL GNZINIT(-MEMGNZ)
*
*   output to file on stream 28
*
       CALL GNZOPEN(28,'OUTP','mygenz.output')
*
* ---- END GENZ ---  (MAIN)
*
       return
       end
