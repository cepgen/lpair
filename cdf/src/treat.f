      DOUBLE PRECISION FUNCTION treat(f,x,ndim)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      EXTERNAL f
      DOUBLE PRECISION x(10),z(10)
      PARAMETER(MAXDIV=250)
      COMMON/bveg2/ndo,it,si,si2,swgt,schi,xi(MAXDIV,10),scalls
     + ,d(MAXDIV,10),di(MAXDIV,10),nxi(MAXDIV,10)
      SAVE ncall,r
      DATA ncall/0/
      IF ( ncall .EQ. 0 ) THEN
        ncall = 1
        r = ndo
        r = r**ndim
      ENDIF
      w = r
      do 4 i = 1,ndim
        xx = x(i)*ndo
        j  = xx
        jj = j+1
        y  = xx-j
        IF ( j .GT. 0 ) THEN
           dd = xi(jj,i)-xi(j,i)
        ELSE
           dd = xi(1,i)
        ENDIF
        z(i) = xi(jj,i)-dd*(1-y)
        w = w*dd
4     CONTINUE
      treat = w*f(z)
      RETURN
      END
