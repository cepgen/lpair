*-- Author :    ZEUS Offline Group   18/08/94
*******************************************************************
*
      FUNCTION VGRAN(DUMMY)
*
*   USING H1RN FOR VEGAS, SETGEN AND GENERA
*
* JC changed to use Geant generator for easy random number series control 
*
*******************************************************************
      integer buf_size,buf_size1
      parameter (buf_size = 10000)
      parameter (buf_size1 = buf_size+1)
      REAL*4 H1RN,BUFF(buf_size)
      REAL*8 VGRAN
      integer i,j
      data j / 0 /
      if(j.eq.0) then
        call grndm(buff,buf_size)
        j = 1
      endif
      if(j.lt.buf_size1) then
        vgran = buff(j)
        j = mod(j+1,buf_size1)
      endif  
*      VGRAN=RLU(DUMMY)
      END
