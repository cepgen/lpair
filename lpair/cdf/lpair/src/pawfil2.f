      subroutine pawfil2(id,xmin,xmax,weight)
*
* --- LPAIR data common block
*
       INTEGER IPAR(20)
       REAL*8 LPAR(20)
       COMMON/DATAPAR/IPAR,LPAR
*
      integer id
      real xvalue,xmin,xmax,weight
*
      IF(IPAR(2).EQ.0) RETURN
*
      xvalue = (xmin+xmax)/2.
*
      CALL HCDIR('//HISTO',' ')
      CALL HFILL(ID,REAL(XVALUE),0.,REAL(WEIGHT))
*
      return
      end
